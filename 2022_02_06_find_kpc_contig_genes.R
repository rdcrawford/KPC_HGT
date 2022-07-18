# ------------------------------------------------------------------------------
# Find KPC contig genes
# 2022/02/06
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# Find the KPC plasmid genes using the contigs of complete genomes
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( ape )
library( Rcpp )
library( maps )
library( cognac )
library( future )
library( cluster )
library( phytools )
library( future.apply )

# ---- Constant declarations ---------------------------------------------------

# Column indices of the gff files
G_ID       = 1
G_ANNOT    = 2
CONTIG     = 3
G_START    = 4
G_END      = 5
G_CLUST    = 7
G_NT_CLUST = 8

# ---- Load the data -----------------------------------------------------------

# Read in the relevant meta-data on the analysis genomes
isolateMetaData = read.table(
  "../data/isolateMetaData.tsv",
  sep              = '\t',
  header           = TRUE,
  stringsAsFactors = FALSE
  )

# Load the "geneEnv" with the parsed gff files and cd-hit results
load( "../data/2022_02_04_parse_plasmid_gene_algns.rData" )

# ---- Function definitions ----------------------------------------------------

GetIdFromAnnot = function( annot, geneEnv )
{
  annotIdxs = lapply( geneEnv$gfList,
    function(x) grep( annot, x[ , G_ANNOT ] )
    )
  gIds = unlist( sapply( 1:length( annotIdxs ),
    function(i) geneEnv$gfList[[i]][ annotIdxs[[i]], G_NT_CLUST ]
    ))
  gCounts = table( gIds )
  gCounts = gCounts[ order( gCounts, decreasing = TRUE ) ]
  print( gCounts )

  return( as.integer( names( gCounts )[1] ) )
}

# ---- Find the genes on each KPC plasmid --------------------------------------

# Find which genomes are complete
compGenomes =
  isolateMetaData$genomeNames[ isolateMetaData$genomeStatus == "complete" ]
compIdxs = which( geneEnv$genomeNames %in% compGenomes )

# Find the gene id assigned to KPC
kpc = GetIdFromAnnot( "KPC", geneEnv )

# Make a list of genes on each KPC+ contig
kpcPlGeneList = list()
for ( i in compIdxs )
{
  kpcIdxs    = which( geneEnv$gfList[[i]][ , G_NT_CLUST ] == kpc )
  kpcContigs = unique( geneEnv$gfList[[i]][ kpcIdxs, CONTIG ] )

  if ( length( kpcContigs ) )
  {
    contGenes = lapply( kpcContigs, function(x)
    {
      isContig = geneEnv$gfList[[i]][ , CONTIG ] == x
      return( geneEnv$gfList[[i]][ isContig, G_NT_CLUST ] )
    })

    names( contGenes ) = paste0( geneEnv$genomeNames[i], '_', kpcContigs )
    kpcPlGeneList = append( kpcPlGeneList, contGenes )
  }
}

# Remove any instances where KPC was on the chromosome
nContigGenes  = sapply( kpcPlGeneList, length )
kpcPlGeneList = kpcPlGeneList[ nContigGenes < 2000 ]

# ---- Make a presence absence matrix for each KPC contig ----------------------

# Make a vector with the gene Ids of the kpc plasmid genes
genes = unique( unlist( kpcPlGeneList ) )
genes = genes[ genes != 0 ]
genes = genes[ order( genes ) ]

# Make a plasmid x gene matrix with the presence of each plasmid gene on the
# plasmid
kpcContigMat = sapply( genes,
  function(x) as.integer( sapply( kpcPlGeneList, function(y) x %in% y ) )
  )
row.names( kpcContigMat ) = names( kpcPlGeneList )
colnames( kpcContigMat )  = genes

kpcPlGeneMat = plGeneMat[ , genes ]
save( file = "../data/kpcPlGeneMat.rData", list = "kpcPlGeneMat" )

# ---- Save the data -----------------------------------------------------------

print( pryr::mem_used() )
rData = "../data/2022_02_06_find_kpc_contig_genes.rData"
save( file = rData, list = c( "kpc", "kpcPlGeneList", "kpcContigMat" ) )

# ------------------------------------------------------------------------------