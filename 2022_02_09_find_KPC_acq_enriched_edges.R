# ------------------------------------------------------------------------------
# Find cluster KPC acquisition enrichment
# 2022/02/09
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# Use a permutation test to determine clusters which are statistically
# enriched for KPC acquisition
# ------------------------------------------------------------------------------

set.seed( 144 )

# ---- Load libraries ----------------------------------------------------------

library( ape )
library( future )
library( pheatmap )
library( future.apply )

# ---- Constant declarations ---------------------------------------------------

N_PERMUTES  = 1e6 # Number of permutations to preform
MIN_EDGES   = 5   # Min number of edges required for statistical analysis
IS_KPC      = 1   # Transition edge data column index indicating if edge is KPC+
KPC_DIR     = 2   # Transition edge data column index indicating KPC transition
EDGE_ID     = 2   # Transition edge data column index index of the edge
TR_DIR      = 3   # Transition edge data column index direction of transition
NO_KPC_EDGE = 1   # Transition occurs on an edge without KPC
KPC_EDGE    = 2   # Transition occurs on an edge with KPC
KPC_CO_TR   = 3   # Simultaneous transition with KPC in the same direction

# ---- Import the data ---------------------------------------------------------

# Load the data on the kpc plasmids: "kpc", "kpcPlGeneList", and
# "kpcContigMat"
load( "../data/2022_02_06_find_kpc_contig_genes.rData" )

# Load the data on the cluster transition edges: "geneEdgeMat", "clTrData",
# "clMat", "mbrList", "clAcqMat", "clAcqEdgeList", and "clDelEdgeList"
load( "../data/2022_02_08_get_cl_gene_counts_by_edge.rData" )

# Load the KPC plasmid gene MPR results
load( "../data/2022_02_07_run_mpr_kpc_plasmid_genes.rData" )

# ---- Function definitions ----------------------------------------------------

GetEdgeStatus = function( edgeIdx, trDir )
{
  if ( trDir == '+' )
  {
    if ( edgeIdx %in% kpcAqEdges ) return( KPC_CO_TR )
    if ( edgeIdx %in% kpcEdges ) return( KPC_EDGE )
    return( NO_KPC_EDGE )
  }

  if ( edgeIdx %in% kpcDelEdges ) return( KPC_CO_TR )
  if ( edgeIdx %in% kpcEdges ) return( KPC_EDGE )
  return( NO_KPC_EDGE )
}

SampleKpcEdgeStatus = function( nEdges, fracAcqs )
{
  # Sample transition directions adjusted for the probability of acquisition
  # and deletion
  sampleTrDir = sample(
    x       = c( "-", "+" ),
    size    = nEdges,
    prob    = fracAcqs,
    replace = TRUE
    )

  # Sample edges adjusted for the length of the branches
  sampleEdges = sample(
    x       = edgeIds,
    size    = nEdges,
    prob    = cgTree$edge.length,
    replace = TRUE
    )

  # Set the status of the ransomized data
  edgeStatus = sapply( seq( nEdges ),
    function(i) GetEdgeStatus( sampleEdges[i], sampleTrDir[i] )
    )

  # Return the counts of each transition type
  return( CountEdgeTypes( edgeStatus ) )
}

CountEdgeTypes = function( edgeStatus )
{
  return( sapply( seq( 3 ), function(x) sum( edgeStatus == x ) ) )
}

CountEdgeClasses = function( countList )
{
  clCounts = sapply( seq( 3 ),
    function(i) sum( sapply( countList, function(x) x[i] ) )
    )
  return( clCounts )
}

# ---- Preform the permutation test --------------------------------------------

kpcIdx      = which( colnames( trMat ) == kpc )
kpcAqEdges  = which( trMat[ , kpcIdx ] == ACQUISTION )
kpcDelEdges = which( trMat[ , kpcIdx ] == DELETION )
kpcEdges    = which( geneEdgeMat[ , kpcIdx ] == 1 )
clIds       = unique( clTrData$clId )
clIdxs      = lapply( clIds, function(x) which( clTrData$clId == x ) )
nClEdges    = sapply( clIdxs, length )

# Find all of the gene clusters with at least n transitions
idxs = which( nClEdges >= MIN_EDGES )

# Permute the transition edges
edgeIds = seq( nrow( cgTree$edge ) )
expKpcEdgeStats = future_lapply( idxs, function(i)
{
  # Get the fraction of these clusters that occur on acquisition and deletion
  # edges
  fracAcqs = sapply( c( "-", "+" ),
    function(x) sum( clTrData[ clIdxs[[i]], TR_DIR ] == x ) / nClEdges[i]
    )

  # Get the status of the randomized edges
  rndData = lapply( seq( N_PERMUTES ),
    function(j) SampleKpcEdgeStatus( nClEdges[i], fracAcqs )
    )
  return( rndData )
}, future.seed = NULL )

# Get the status of each edge
kpcStatus = lapply( clIdxs[ idxs ], function(x)
{
  clEdgeStats = sapply( x,
    function(i) GetEdgeStatus( clTrData[ i, EDGE_ID ], clTrData[ i, TR_DIR ] )
    )
  return( clEdgeStats )
})

# Find the KPC status of the cluster edges
obsKpcEdgeStats = lapply( kpcStatus, CountEdgeTypes )

# Calculate the p-value for each edge
pValsTrEdge = sapply( 1:length( expKpcEdgeStats ), function(i)
{
  nObs  = obsKpcEdgeStats[[i]][ KPC_CO_TR ]
  nExp  = sapply( expKpcEdgeStats[[i]], function(x) x[ KPC_CO_TR ] )
  nGrtr = sum( nExp > nObs )
  return( nGrtr / N_PERMUTES )
})

# Get the bonferroni corrected significance threshold
sigThresh = 0.05 / length( idxs )

cat( "Corrected significant hits:", sum( pValsTrEdge <= sigThresh ), '\n' )
cat( "Uncorrected significant hits:", sum( pValsTrEdge <= 0.05 ), '\n' )
cat( "P-value summary:\n")
summary( pValsTrEdge )


obsCounts = CountEdgeClasses( obsKpcEdgeStats )

nCls = length( idxs )
expCounts = sapply( seq( N_PERMUTES ), function(i)
{
  classCounts = lapply( seq( nCls ), function(j) expKpcEdgeStats[[j]][[i]] )
  return( CountEdgeClasses( classCounts ) )
})

print( obsCounts )
summary( expCounts[ 1, ] )
summary( expCounts[ 2, ] )
summary( expCounts[ 3, ] )

# No KPC edges: 2051 observed
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   827.0   865.5   880.5   880.2   897.5   936.0
# KPC+ edges: 2045 observed
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    3196    3233    3250    3252    3270    3305
# KPC/Gene cluster co-transitions: 233 observed
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    6.00   12.75   15.00   15.26   18.00   23.00

# ---- Save the data -----------------------------------------------------------

print( pryr::mem_used() )
rData = "../data/2022_02_09_find_cluster_kpc_acquisition_enrichment.rData"
save( file = rData, list = ls() )

rData = "../data/2022_02_09_find_cluster_kpc_acquisition_enrichment_sig.rData"
save( file = rData, list = c( "pValsTrEdge", "sigThresh", "idxs" ) )

# ---- Make the figures --------------------------------------------------------

figDir = "../figures/2022_02_09_find_KPC_acq_enriched_edges/"

clKpcStat = sapply( obsKpcEdgeStats, function(x)
{
  isKpc = x != NO_KPC_EDGE
  if ( !TRUE %in% isKpc ) return( 1 )
  if ( !FALSE %in% isKpc ) return( 3 )
  return( 2 )
})
kpcClStatus = sapply( seq( 3 ), function(x) sum( clKpcStat == x ) )
names( kpcClStatus ) =
  paste( c( "Never", "Sometimes", "Always" ), "present\non KPC edge" )
par( mar = c( 12, 4.1, 4.1, 2.1 ) )
png(
  filename = paste0( figDir, "2022_02_08_edge_kpc_prevalance.png" ),
  units    = "in",
  width    = 5,
  height   = 5,
  res      = 300
  )
# par( mfrow = c( 1, 1) )
barplot( kpcClStatus, col = "#2E4769", border = "#2E4769", las = 2 )
dev.off()

isAcq = sapply( kpcStatus, function(x) KPC_CO_TR %in% x )
nAcq  = c( sum( !isAcq ), sum( isAcq ) )
names( nAcq ) =
  c( "No KPC acquisitions\nobserved", "KPC acquisitions\nobserved" )

png(
  filename = paste0( figDir, "2022_02_09_n_cls_with_kpc_acq.png" ),
  units    = "in",
  width    = 5,
  height   = 5,
  res      = 300
  )
par( mar = c( 10, 4.1, 4.1, 2.1 ) )
barplot( nAcq, col = "#2E4769", "border" = "#2E4769", las = 2 )
dev.off()

nVal = 1 / N_PERMUTES
ltPvs = -log10( pValsTrEdge + nVal )
ltPvs = ltPvs[ order( ltPvs ) ]
ltSigTrh = -log10( sigThresh + nVal )
png(
  filename = paste0( figDir, "2022_02_09_kpc_enrichment_pvalues.png" ),
  units    = "in",
  width    = 5,
  height   = 5,
  res      = 300
  )
par( mar = c( 5.1, 4.1, 4.1, 2.1 ) )
cols = ifelse( ltPvs >= ltSigTrh, "#952827", "#2E4769" )
barplot(
  ltPvs,
  col    = cols,
  border = cols,
  xlab   = "Gene cluster",
  ylab   = "-log10 p-value"
  )
legend(
  "topleft",
  col = c( "#952827", "#2E4769" ),
  pch = 15,
  legend = c( "Significant", "Not significant" )
  )
dev.off()

expCoTrCounts = expCounts[ 3, ]
obsCoTr       = obsCounts[3]



BIN_VAL = 5

maxVal = ceiling( obsCoTr / 25 ) * 25
obsBin = obsCoTr / BIN_VAL
binVals   = seq( 0, maxVal, by = BIN_VAL )
binCounts = sapply( binVals,
  function(x) sum( x <= expCoTrCounts & expCoTrCounts < x + BIN_VAL  )
  )

binCols = ifelse( binCounts > 0, "#395983", "white" )


barplot(
  binCounts,
  col    = binCols,
  border = binCols
  )sharted
abline(
  v = obsBin,
  col = "#952827",
  lwd = 2,
  lty = 2
  )
points(
  x   = obsBin,
  y   = max( binCounts ) * 0.975,
  col = "#952827",
  pch = 8,
  cex = 3
  )
axisBins = seq( 0, maxVal /  25 )
atVals   = axisBins * BIN_VAL
labVals   = axisBins * 25
axis(
  side   = 1,
  at     = atVals,
  labels = labVals,
  las    = 2,
  outer  = FALSE
  )

# ------------------------------------------------------------------------------
# kpcIsos = isolateMetaData$genomeNames[ isolateMetaData$kpc == 1]
# fracKpc = sapply( 1:ncol( clMat ), function(j)
# {
#   isCl = clMat[ , j ] == 1
#   nKpc = sum( row.names( clMat )[ isCl ] %in% kpcIsos )
#   return( nKpc / sum( isCl ) )
# })
#
# cgTree = read.tree( "../data/cgTree.tre" )
# # Load the annotations and colors: "annots" and "annotCols"
# load( "../data/2022_02_07_make_annotations_and_colors.rData" )
#
#
# HighlightEdges = function( featureEdges )
# {
#   edgeCols = sapply( 1:nrow( cgTree$edge ), function(x)
#   {
#     if ( x %in% featureEdges ) return( "#952827" )
#     return( "black" )
#   })
#
#   plot(
#     cgTree,
#     type            = "p",
#     label.offset    = 100,
#     cex             = 0.001,
#     no.margin       = TRUE,
#     show.tip.label  = TRUE,
#     use.edge.length = FALSE,
#     edge.color      = edgeCols,
#     edge.width      = 2
#     )
#
#   # # The "adjVal" variable sets the position of the tip labels. The "adjInc"
#   # # sets how much to move the tip labels to the right
#   # adjVal = 100
#   # adjInc = 75
#   # rowIdxs = GetAnnotRowIdxs( cgTree, annots )
#   # for ( j in c( 1) )
#   # {
#   #   PlotTipLabels( rowIdxs, annots, j, annotCols, adjVal )
#   #   adjVal = adjVal + adjInc
#   # }
# }
#
# source( "R/PlotTree.R" )
# source( "R/TreeParsingFunctions.R" )
#
#
# isTip = cgTree$edge[ kpcEdges, 2 ] <= length( cgTree$tip.label )
# xx = cgTree$tip.label[ cgTree$edge[ kpcEdges[ isTip ], 2 ] ]
# isXx = kpcIsos %in% xx
# sum( !isXx )
# xx[ !xx %in% kpcIsos ]
# cgTree$edge[ kpcEdges, ]
#
#
# clEdges = GetFeatureEdges( cgTree, clAcqEdgeList[[1]], clDelEdgeList[[1]] )
# kpcEdges = GetFeatureEdges( cgTree, kpcAqEdges, kpcDelEdges )
# HighlightEdges( kpcEdges )
# HighlightEdges( clEdges )
# ------------------------------------------------------------------------------
