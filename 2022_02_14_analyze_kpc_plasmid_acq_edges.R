# ------------------------------------------------------------------------------
#
# 2022/02/14
# Ryan D. Crawford
# ------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( ape )
library( Rcpp )
library( maps )
library( cognac )
library( future )
library( cluster )
library( phytools )
library( pheatmap )
library( future.apply )

# ---- Constant declarations ---------------------------------------------------

# Column indices of the cluster stat data
G_ID       = 1
C_ID       = 2
CL_CONTIG  = 3
N_CL_GENES = 4
N_TP       = 5
N_FP       = 6
PPV        = 7

# ---- Import the data ---------------------------------------------------------

# Load the data on relation between clusters and the plasmids: "clStats",
# "chrContigs", "kpcContigs", and "kpcPlClstrMat"
load( "../data/2022_02_11_find_cluster_precision.rData" )

# Load the data on kpc contigs
load( "../data/2022_02_06_find_kpc_contig_genes.rData" )

# Load the data on the gene clusters: "geneEdgeMat", "clTrData", "clMat",
# "mbrList", "clAcqMat", "clAcqEdgeList", and "clDelEdgeList"
load( "../data/2022_02_08_get_cl_gene_counts_by_edge.rData" )

# ---- Function definitions ----------------------------------------------------


# ----

nCls = length( mbrList )
nPls = length( )
dim( clAcqMat )
dim( kpcPlClstrMat )

plGIds = sapply( row.names( kpcContigMat ),
  function(x) strsplit( x, '_' )[[1]][1], USE.NAMES = FALSE
  )

gIdxs = sapply( plGIds, function(x) which( row.names( clAcqMat ) == x ) )

plAcqEdges = lapply( 1:nrow( kpcPlClstrMat ), function(i)
{
  plCls = which( kpcPlClstrMat[ i, ] == 1 )
  if ( length( plCls ) == 0 ) return( NULL )
  return( clAcqMat[ gIdxs[i], plCls ] )
})

nPlAcqEdges = sapply( plAcqEdges, length )
isMultiCl =  which( nPlAcqEdges > 1 )
edgeCountMat = sapply( isMultiCl, function(i)
{
  edgeVals = unique( plAcqEdges[[i]] )
  if ( length( edgeVals )  == 1 ) return( c( length( plAcqEdges[[i]] ), 0 ) )
  nEdge    = sapply( edgeVals, function(x) sum(  plAcqEdges[[i]] == x ) )
  edOrd    = order( nEdge, decreasing = TRUE )
  nEdge    = nEdge[ edOrd ]
  return( c( nEdge[1], sum( nEdge[ 2:length( nEdge ) ] ) ) )
})

nAcqEdges     = sapply( plAcqEdges, length )
acqEdgeCounts = lapply( plAcqEdges, table )


# ---- Save the data -----------------------------------------------------------

print( pryr::mem_used() )
rData = "../data/.rData"
save( file = rData, list = c( "clStats", "chrContigs", "kpcContigs" ) )

# ---- Make the figures --------------------------------------------------------

figDir = "../figures/"

nn = sapply( seq( 0, max( nPlAcqEdges ) ), function(x) sum( nPlAcqEdges == x ) )
names( nn ) = seq( 0, max( nPlAcqEdges ) )
cols = ifelse( nn == 0, "white", "#395983" )
barplot(
  nn,
  col    = cols,
  border = cols
  )
sum( edgeCountMat[ 2, ] == 0 ) / 239
plEdgeCounts = c(
  "On majority plasmid\ntrasition edge" = sum( edgeCountMat[ 1, ] ),
  "Not on majority plasmid\ntrasition edge" = sum( edgeCountMat[ 2, ] )
  )

barplot( plEdgeCounts, col = "#395983", border = "#395983" )
# ------------------------------------------------------------------------------