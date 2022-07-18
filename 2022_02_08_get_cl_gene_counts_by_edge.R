# ------------------------------------------------------------------------------
# Get Cluster gene counts by edge
# 2022/02/08
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# For each edge, find the number of genes in each cluster which are present
# of the edge. Then use these data to find the edge on which a cluster was
# acquired in each genome.
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( ape )
library( maps )
library( future )
library( cluster )
library( pheatmap )
library( phytools )
library( future.apply )

# ---- Constant declarations ---------------------------------------------------

AQUISTION = 1  # Integer indicating acquisitions in "trMat"
DELETION  = 2  # Integer indicating deletions in "trMat"
CL_ABSENT = 0  # Value specifying that a gene isn't present on any edge

# ---- Import the data ---------------------------------------------------------

# Load the MPR data, kpc plasmids data, and the data on the presence/absence
# of clusters across the analysis genomes
load( "../data/2022_02_07_find_best_hclust_height_threshold.rData" )
load( "../data/2022_02_07_run_mpr_kpc_plasmid_genes.rData" )
load( "../data/2022_02_06_find_kpc_contig_genes.rData" )

if ( nrow( cgTree$edge ) != nrow( trMat ) )
  stop( "Number of edges in the tree and transition matrix are not the same!" )

# ---- Function definitions ----------------------------------------------------

source( "R/TreeParsingFunctions.R" )

source( "R/GetTreeOrder.R" )

FindAcqEdge = function( edgeVals, geneCounts, minMbrs )
{
  if ( geneCounts[1] < minMbrs ) return( CL_ABSENT )
  idx = 1
  while( geneCounts[ idx + 1 ] >= minMbrs && idx < length( geneCounts ) )
    idx = idx + 1
  return( edgeVals[ idx ] )
}

GetDeletionEdges = function( tree, acqEdges, clEdgeCounts, minMbrs )
{
  delEdges = lapply( acqEdges,
    function(x) FindDelEdge( tree, x, clEdgeCounts, minMbrs )
    )
  return( unique( unlist( delEdges ) ) )
}

FindDelEdge = function( tree, edgeVal, clEdgeCounts, minMbrs )
{
  # If the feature was lost on this edge of the tree, return null value
  if ( clEdgeCounts[ edgeVal ] < minMbrs ) return( edgeVal )

  # If this corresponds to a tip, return the tip label
  if ( tree$edge[ edgeVal, 2 ] <= Ntip( tree ) ) return( NULL )

  # Get the edges extending from the RHS node
  nextEdges = LookUpNextEdges( tree, edgeVal )

  # Return the deletions on the subsequent edges
  delEdges = c(
    FindDelEdge( tree, nextEdges[1], clEdgeCounts, minMbrs ),
    FindDelEdge( tree, nextEdges[2], clEdgeCounts, minMbrs )
    )
  return( delEdges )
}

RepClId = function( lst )
{
  return( lapply( 1:length( lst ), function(i) rep( i, length( lst[[i]] ) ) ) )
}

# ---- Find how many genes in each cluster are present on each edge ------------

# Make an edge x gene matrix to store the presence of each gene on each gene
geneEdgeMat = matrix( 0, nrow = nrow( trMat ), ncol = ncol( trMat ) )
for ( j in 1:ncol( geneEdgeMat ) )
{
  # Find all of the edges that a gene is present on
  geneEdges = GetFeatureTransitionEdges( cgTree, trMat[ , j ] )

  # For each edge where a gene is predicted to be, indicate it as present
  # in the matrix
  for ( i in geneEdges ) geneEdgeMat[ i, j ] = 1
}
dimnames( geneEdgeMat ) = dimnames( trMat )

# Make an edge x gene cluster matrix with the number of cluster genes
# present on each edge
clGeneEdgeMat = sapply( mbrList, function(x) rowSums( geneEdgeMat[ , x ] ) )

# ---- Find Where each genome acquired each cluster ----------------------------

# Look up the index of the edges connected to each tip
tipEdges = sapply( cgTree$tip.label, function(x) GetTipEdge( x, cgTree ) )

# Get all of the edges down to the root
ancEdges = lapply( tipEdges, function(x) GetAncestorEdges( cgTree, x, NULL ) )

# Make a genome x cluster matrix to store the acquisition edge for each
# genome.
nClMbrs   = sapply( mbrList, length )
minClMbrs = ceiling( nClMbrs * 0.51 )

# Make a genome x cluster matrix to store the edge of acquisition
clAcqMat = sapply( 1:ncol( clGeneEdgeMat ), function(j)
{
  acqEdge = sapply( ancEdges,
    function(x) FindAcqEdge( x, clGeneEdgeMat[ x, j ], minClMbrs[j] )
    )
  return( acqEdge )
})
row.names( clAcqMat ) = cgTree$tip.label

# isEip = cgTree$tip.label %in% eipData$Sample_ID
clAcqEdgeList = lapply( 1:ncol( clAcqMat), function(j)
{
  edgeVals = unique( clAcqMat[ , j ] )
  return( edgeVals[ edgeVals != CL_ABSENT ] )
})
cat( "clAcqEdgeList:", length( clAcqEdgeList ), '\n' )

# Find the cluster deletion edges
clDelEdgeList = lapply( 1:length( clAcqEdgeList ),
  function(j)
    GetDeletionEdges(
      cgTree, clAcqEdgeList[[j]], clGeneEdgeMat[ , j ], minClMbrs[j]
      )
  )

nClAcqs = sapply( clAcqEdgeList, length )
nClDels = sapply( clDelEdgeList, length )

clMat =
  sapply( 1:ncol( clAcqMat ), function(j) as.integer( clAcqMat[, j ] > 0 ) )
row.names( clMat ) = cgTree$tip.label

clIdxs = unlist( c( RepClId( clAcqEdgeList ), RepClId( clDelEdgeList ) ) )
clIds = sapply( clIdxs, function(x) isGood[x])

# Make a data-frame to store the relevant information on the clusters
clTrData = data.frame(
  clIds,
  unlist( c( clAcqEdgeList , clDelEdgeList ) ),
  c( rep( '+', sum( nClAcqs ) ), rep( '-', sum( nClDels ) ) )
  )
colnames( clTrData ) = c( "clId", "edgeIdx", "trDir" )


# ---- Save the data -----------------------------------------------------------

print( pryr::mem_used() )
rData = "../data/2022_02_08_get_cl_gene_counts_by_edge.rData"
save(
  file = rData,
  list = c(
    "geneEdgeMat", "clTrData", "clMat", "mbrList", "clAcqMat",
    "clAcqEdgeList", "clDelEdgeList"
    )
  )

# ---- Make the figures --------------------------------------------------------

figDir = "../figures/"

# Read in the cognac fastTree
cgTree = midpoint.root( read.tree( "../data/cgTree.tre" ) )
nCls      = length( mbrList )
nTotalTrs = sapply( seq( nCls ), function(i) sum( nClDels[i], nClAcqs[i] ) )
plOrd     = order( nTotalTrs )

par( mfrow = c( 2, 1 ) )
yLim = c( 0, 200 )
png(
  filename = paste0( figDir, "2022_02_08_tr_frequency_by_cluster.png" ),
  units    = "in",
  width    = 5,
  height   = 5,
  res      = 300
  )
barplot(
  nClAcqs[ plOrd ], col = "#395983", border = "#395983", ylim = yLim,
  main = "Number of acquisitions by cluster"
  )
barplot(
  nClDels[ plOrd ], col = "#395983", border = "#395983", ylim = yLim,
  main = "Number of deletions by cluster"
  )
dev.off()

load( "../data/2022_02_07_make_annotations_and_colors.rData" )
load( "../data/eipData.rData" )

hmCols  = colorRampPalette( c( "#395983", "#e7e7e7", "#952827" ) )( 100 )
treeOrd = GetTreeOrder( cgTree, row.names( clMat ) )
treeOrd = GetTreeOrder( midpoint.root( cgTree ), row.names( clMat ) )
clMat   = clMat[ treeOrd, ]
pheatmap(
  clMat[ , isCl ],
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  color             =  c( "#e7e7e7", "#2E4769" ),
  annotation_row    = annots,
  annotation_colors = annotCols,
  cluster_rows      = FALSE,
  cluster_cols      = TRUE,
  filename          = paste0( figDir, "2022_02_08_cl_mat.png" )
  )

isEip    = row.names( clMat ) %in% eipData$Sample_ID
isCl     = colSums( clMat[ isEip, ] ) > 1
eipClMat = clMat[ isEip, isCl ]


pheatmap(
  eipClMat,
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  color             =  c( "#e7e7e7", "#2E4769" ),
  annotation_row    = annots,
  annotation_colors = annotCols,
  cluster_rows      = FALSE,
  cluster_cols      = TRUE,
  filename          = paste0( figDir, "2022_02_08_eip_cl_mat.png" )
  )

# st307Isos = row.names( annots )[ annots$ST == "K. pneumoniae ST307" ]
# isSt307   = row.names( eipClMat ) %in% st307Isos
#
# is307Cl = colSums( eipClMat[ isSt307, ] ) > 1
# isClIso = rowSums( eipClMat[ , is307Cl ] ) > 0
#
# pheatmap(
#   eipClMat[ isClIso, is307Cl ],
#   show_rownames     = FALSE,
#   show_colnames     = FALSE,
#   color             =  c( "#e7e7e7", "#2E4769" ),
#   annotation_row    = annots,
#   annotation_colors = annotCols,
#   cluster_rows      = FALSE,
#   cluster_cols      = TRUE
#   )
source( "R/MakeDistMatFns.R" )
acqJacMat = MakeJaccardMat( clAcqEdgeList, names( mbrList ) )
hCl = agnes( acqJacMat )
pheatmap(
  acqJacMat[ hCl$order, hCl$order ],
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  color             = hmCols,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  filename          = paste0( figDir, "2022_02_08_cl_edge_jac_mat.png" )
  )

nClTrEdge = lapply( 1:ncol( clGeneEdgeMat ),
  function(j) clGeneEdgeMat[ clAcqEdgeList[[j]], j ]
  )

fracClGenes = lapply( 1:length( nClTrEdge ), function(j)
  nClTrEdge[[j]] / nClMbrs[j]
  )

nClGs  = unlist( nClTrEdge )
frVals = unlist( fracClGenes )
names( nClGs ) = NULL
names( frVals ) = NULL
frOrd  = order( frVals )

cat( "Fraction of cluster genes at transition edge:" )
summary( frVals )

par( mfrow = c( 2, 1 ) )
png(
  filename = paste0( figDir, "2022_02_08_n_cl_genes_frac_cl_genes.png" ),
  units    = "in",
  width    = 5,
  height   = 5,
  res      = 300
  )
barplot(
  nClGs[ frOrd ],
  col    = "#395983",
  border = "#395983",
  xlab   = "Cluster acquisition edge",
  ylab   = "Number of genes"
  )
barplot(
  frVals[ frOrd ],
  col    = "#395983",
  border = "#395983",
  xlab   = "Cluster acquisition edge",
  ylab   = "Fraction of cluster genes"
  )
dev.off()



# ------------------------------------------------------------------------------
