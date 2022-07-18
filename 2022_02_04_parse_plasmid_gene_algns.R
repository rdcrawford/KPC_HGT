# ------------------------------------------------------------------------------
#
#
# Ryan D. Crawford
# ------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( cognac )
library( Rcpp )

# ---- Constant declarations ---------------------------------------------------

IDENT_THRESH = 0.01
SUB_THRESH   = 2

# ---- Load the data -----------------------------------------------------------

# Load the data with the path to the neucleotide alignments
# load( "../data/2022_02_04_make_nt_aligns.rData" )
load(  "../data/2022_02_04_2022_02_04_parse_cd_hit_data.rData" )
algnDir   = "../analysis/2022_02_04_eip_cognac_analysis/nt_gene_algns/"
algnExt   = "_all_seqs_mafft_alignnt.fasta"
nGenes    = 72145
algnPaths = paste0( algnDir, "gene_cluster_", seq( nGenes ), algnExt )

for ( i in 1:length( algnPaths ) )
{
  if ( !file.exists( algnPaths[i] ) )
    algnPaths[i] = paste0( algnDir, "gene_cluster_", i, ".fasta" )
}

# ---- Function definitions ----------------------------------------------------

source( "R/VecToList.R" )

FindGeneClusters = function( algnPath )
{
  distMat   = CreateAlgnDistMat( algnPath, "shared" )
  nSubMat   = CreateAlgnDistMat( algnPath, "raw" )
  nSeqs     = nrow( distMat )
  clReps    = vector( "integer", 0 )
  geneClIds = vector( "integer", nSeqs )

  # Iterate over rows in the distance matrix and make the clusters
  # based off of the distance to the cluster representitive.
  for ( i in seq( nSeqs ) )
  {
    bestCl       = GetGeneMatch( i, distMat, nSubMat, clReps )
    geneClIds[i] = bestCl
    if ( bestCl > length( clReps ) ) clReps[ length( clReps ) + 1 ] = i
  }

  # If all sequences are in the same cluster, return the seq ids
  if ( length( clReps ) == 1 ) return( list( row.names( distMat ) ) )

  # Reevaluate the cluster assignments and pick the best representative
  # for each sequence
  repSubs = nSubMat[ , clReps ]
  geneClIds = sapply( seq( nSeqs ), function(i) which.min( repSubs[ i, ] ) )
  names( geneClIds ) = row.names( distMat )

  # Return the list of gene ids for
  return( VecToList( geneClIds ) )
}

GetGeneMatch = function( i, distMat, nSubMat, clReps )
{
  if ( length( clReps ) == 0 ) return( 1 )
  bestIdx = 0
  bestMatch = 1
  for ( j in 1:length( clReps ) )
  {
    if (
      distMat[ i, clReps[j] ] <= IDENT_THRESH ||
      nSubMat[ i, clReps[j] ] <= SUB_THRESH
      )
    {
      if ( bestMatch > distMat[ i, clReps[j] ] )
      {
        bestMatch = distMat[ i, clReps[j] ]
        bestIdx   = j
      }
    }
  }
  if ( bestIdx == 0 ) return( length( clReps ) + 1 )
  return( bestIdx )
}

# Add the cd-hit cluster name to the gff file
AddClustIdxToGff = function( geneEnv, gIdx, clustList, gIdList )
{
  G_ID = 1
  # Initialize a vector to store the cluster id for each gene
  clIdxs = vector( "integer", nrow( geneEnv$gfList[[ gIdx ]] ) )

  # For each cd-hit cluster...
  for ( i in  seq( nGenes ) )
  {
    # Find if the current genome has this gene
    isThisGenome = gIdList[[i]] == geneEnv$genomeNames[ gIdx ]

    # Match the gene id to the corresponding row in the gff file, and
    # mark the corresponding cluster as this gene index
    for ( j in which( isThisGenome ) )
    {
      isClGene =
        geneEnv$gfList[[ gIdx ]][ , G_ID ] == clustList[[i]][j]
      clIdxs[ isClGene ] = i
    }
  }

  # Add cluster indexes to the data frame
  geneEnv$gfList[[gIdx]] = cbind.data.frame(
    geneEnv$gfList[[ gIdx ]],
    "ntClIds" = clIdxs,
    stringsAsFactors = FALSE
    )
}

# ---- Parse the alignments and call orthologous genes -------------------------

clustList = vector( "list", 0 )
for ( algnPath in algnPaths )
{
  geneCls = FindGeneClusters( algnPath )
  sPos    = length( clustList ) + 1
  ePos    = sPos + length( geneCls ) - 1
  clustList[ sPos:ePos ] = geneCls
}

# ---- Make a gene presence absence matrix -------------------------------------

gIdList = lapply( 1:length( clustList ),
  function(i) sapply( clustList[[i]], GetGenomeId )
  )

nGenomes  = length( geneEnv$genomeNames )
nGenes    = length( clustList )
plGeneMat = matrix( 0, nGenomes, nGenes )

for ( j in seq( nGenes ) )
{
  isGenome = geneEnv$genomeNames %in% gIdList[[j]]
  for ( i in which( isGenome ) ) plGeneMat[ i, j ] = 1
}
row.names( plGeneMat ) = geneEnv$genomeNames
colnames( plGeneMat )  = seq( ncol( plGeneMat ) )
save( file = "../data/plGeneMat.rData", list = "plGeneMat" )

# ------------------------------------------------------------------------------

# Add the nt gene id as a column in the gff file
for ( i in seq( nGenomes ) ) AddClustIdxToGff( geneEnv, i, clustList, gIdList )

# ---- Save the data -----------------------------------------------------------

print( pryr::mem_used( ) )
rData = "../data/2022_02_04_parse_plasmid_gene_algns.rData"
save( file = rData, list = c( "clustList", "plGeneMat", "geneEnv" ) )

# ------------------------------------------------------------------------------
