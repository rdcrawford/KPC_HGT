# ------------------------------------------------------------------------------
# Get nucleotide cluster alignment data
# 2022/02/25
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# Alignments based on the cd-hit clusters were divided into multiple
# clusters based on nucleotide identity. Partition each alignment and
# corresponding distance matrix to match the nucleotide clusters.
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( cognac )

# ---- Import the data ---------------------------------------------------------

# Load the paths to the gene alignments: "algnPaths"
load( "../data/2022_02_04_make_nt_aligns.rData" )

# Load the data on the nt gene clusters: "clustList", "plGeneMat", and "geneEnv"
load( "../data/2022_02_04_parse_plasmid_gene_algns.rData" )

# "kpc", "kpcPlGeneList", "kpcContigMat"
load( "../data/2022_02_06_find_kpc_contig_genes.rData" )

# ----- Parse the gene alignments into the corresponding nt clusters -----------

# Initialize the lists to output
kpcPlGenes = as.integer( colnames( kpcContigMat ) )
nGenes     = length( kpcPlGenes )
algns      = vector( "list", nGenes )
geneDists  = vector( "list", nGenes )
kpcClList  = vector( "list", nGenes )
aIdx       = 1
algn       = ParseFasta( algnPaths[ aIdx ] )
distMat    = CreateAlgnDistMat( algnPaths[ aIdx ], "raw" )

# For each nt cluster, get the alignment sequences
for ( i in seq( nGenes ) )
{
  gIdx = kpcPlGenes[i]
  isCl = names( algn ) %in% clustList[[ gIdx ]]

  # If all sequences in the alignment have been classified, advance to the
  # next gene alignment
  while ( !TRUE %in% isCl )
  {
    aIdx    = aIdx + 1
    algn    = ParseFasta( algnPaths[ aIdx ] )
    distMat = CreateAlgnDistMat( algnPaths[ aIdx ], "raw" )
    isCl    = names( algn ) %in% clustList[[ gIdx ]]
  }
  cat( i, ": ", gIdx, " --> ", aIdx, '\n', sep = '' )

  if ( length( clustList[[ gIdx ]] ) > 1 )
  {
    clDistMat      = distMat[ isCl, isCl ]
    seqIdOrd       = order( row.names( clDistMat ) )
    clDistMat      = clDistMat[ seqIdOrd, seqIdOrd ]
    algns[[i]]     = algn[ isCl ][ seqIdOrd ]
    geneDists[[i]] = clDistMat
    kpcClList[[i]] = clustList[[ gIdx ]][ order( clustList[[ gIdx ]] ) ]
  }
}

isCorrect = sapply( seq( nGenes ), function(i)
{
  if ( length( kpcClList[[i]] ) == 1  ) return( TRUE )
  aa = identical( kpcClList[[i]], row.names( geneDists[[i]] ))
  bb = identical( kpcClList[[i]], names( algns[[i]] ) )
  return( aa && bb )
})

if ( FALSE %in% isCorrect )
  stop( "Number incorrect = ", sum( !isCorrect ), '\n' )

names( algns )     = kpcPlGenes
names( geneDists ) = kpcPlGenes
names( kpcClList ) = kpcPlGenes

# ---- Save the data -----------------------------------------------------------

print( pryr::mem_used() )
rData = "../data/2022_02_25_get_nt_cluster_algn_data.rData"
save( file = rData, list = c( "algns", "geneDists", "kpcClList" ) )

# ------------------------------------------------------------------------------