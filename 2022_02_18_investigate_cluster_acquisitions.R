# ---- Import the data ---------------------------------------------------------

# Load the annotations and colors: "annots" and "annotCols"
load( "../data/2022_02_07_make_annotations_and_colors.rData" )

# ---- Function definitions ----------------------------------------------------

PlotClEdgeFracts = function( edgeFracs )
{
  edgePerc = round( edgeFracs * 100 )

  edgeColVals = sapply( edgePerc, function(x)
  {
    if ( x == 0 ) return( "black" )
    return( edgeCols[x] )
  })

  plot(
    cgTree,
    type            = "p",
    label.offset    = 500,
    cex             = 0.001,
    no.margin       = TRUE,
    show.tip.label  = TRUE,
    use.edge.length = FALSE,
    edge.color      = edgeColVals,
    edge.width      = 2
    )

  # The "adjVal" variable sets the position of the tip labels. The "adjInc"
  # sets how much to move the tip labels to the right
  adjVal = 50
  adjInc = 75
  rowIdxs = GetAnnotRowIdxs( cgTree, annots )
  for ( j in 1:ncol( annots ) )
  {
    PlotTipLabels( rowIdxs, annots, j, annotCols, adjVal )
    adjVal = adjVal + adjInc
  }

  legVals = seq( 0, 100, 10 )
  legend(
    "topleft",
    legend = paste0( legVals, '%' ),
    col    = c( "black", edgeCols[ legVals[ 2:11 ] ] ),
    pch    = 15
    )
}

# ----

acqOrd = order( nClAcqs, decreasing = TRUE )

mVals = sapply( acqOrd[ 1:10 ], function(j)w
{
  xx = c(
    median( clGeneEdgeMat[ clAcqEdgeList[[j]], j ] ),
    median( clGeneEdgeMat[ clDelEdgeList[[j]], j ] )
    )
  return( xx / nClMbrs[j] )
})

mVals = sapply( acqOrd, function(j)
{
  edgeVals = clAcqEdgeList[[j]]
  edgeVals = edgeVals[ edgeVals != 1 ]
  pEdges   = sapply( edgeVals, function(x) LookUpPrevEdge( cgTree, x ) )
  eCounts  = clGeneEdgeMat[ edgeVals, j ]
  pCounts  = clGeneEdgeMat[ pEdges, j ]
  geneDiff = sapply( 1:length( edgeVals ), function(i) eCounts[i] - pCounts[i] )
  return( geneDiff )
})
sapply( mVals, function(x) sum( x == 1 ) )
xx = sapply( mVals,
  function(x) sapply( seq(5), function(y) sum( sum( x <= y ) ) )
  )
dim( xx )
yy = sapply( 1:ncol( xx), function(j) xx[ , j ] / nClMbrs[acqOrd[j]] )
summary( yy[,1 ])

sapply( mVals[ 1:10 ], summary )

fracClMat = sapply( 1:ncol( clGeneEdgeMat  ),
  function(j) clGeneEdgeMat[ , j ] / nClMbrs[j]
  )

cols = c( "#ba3231", "#f37a48", "#eaa221", "#00bc48", "#1697B7", "#7f4485")
edgeCols = colorRampPalette( rev( cols ) )( 100 )
barplot( rep( 10, 100 ), col = edgeCols, border = edgeCols )

# ----

for ( j in c( 79, 43 ) )
{
  edgeFracs = fracClMat[ , j ]
  PlotClEdgeFracts( edgeFracs )
}
nClMbrs[  acqOrd[1] ]
clGeneMat = geneEdgeMat[ , mbrList[[ acqOrd[1] ]] ]
clGeneMat = clGeneMat[ rowSums( clGeneMat ) > 1, ]
pheatmap(
  clGeneMat,
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  color             = c( "#e7e7e7", "#395983" ),
  cluster_rows      = TRUE,
  cluster_cols      = TRUE
  )
