# ------------------------------------------------------------------------------
#
# 2022/03/09
# Ryan D. Crawford
# ------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------

set.seed( 144 )

# ---- Load libraries ----------------------------------------------------------

library( vegan )
library( scales )
library( pheatmap )

# ---- Constant declarations ---------------------------------------------------

HGT_EVENT   = 1
IMPORTATION = 2

CL_ID  = 1
ST_ID  = 2
STR_ID = 3
N_EDGE = 4
N_FACS = 5
HGT    = 6

# ---- Import the data ---------------------------------------------------------

# Load the data on the hgt clusters
load( "../data/2022_03_09_analyze_hgt_cluster_distribution.rData" )

# Load the meta-data from the EIP isolates
load( "../data/eipData.rData" )

# ---- Function definitions ----------------------------------------------------

MakeFacIds = function( eipData )
{
  facilityId = vector( "character", nrow( eipData ) )
  for ( st in unique( eipData$State ) )
  {
    facVals = unique( eipData$Fac_ID[ eipData$State == st ] )
    for ( j in 1:length( facVals ) )
    {
      facId = paste0( st, j )
      for ( i in which( eipData$Fac_ID == facVals[j] ) ) facilityId[i] = facId
    }
  }
  names( facilityId ) = eipData$Sample_ID
  return( facilityId )
}

MakeIdxVec = function( idxList )
{
  idxVals = vector( "integer", nrow( annots ) )
  for ( j in 1:length( idxList ) )
    for ( i in idxList[[j]] ) idxVals[i] = j
  return( idxVals )
}

CalcShannonDiv = function( clEdgeData )
{
  isState = lapply( seq( nSts ), function(x)  clEdgeData[ , ST_ID ] == x )
  isStr   = lapply( seq( nStrs ), function(x)  clEdgeData[ , STR_ID ] == x )
  stStrMat = sapply( isState,
    function(x) sapply( isStr, function(y) sum( x & y ) )
    )
  colnames( stStrMat )  = stVals
  row.names( stStrMat ) = strainVals
  return( diversity( t( stStrMat ), index = "shannon" ) )
}

PlotNumProps = function( clEdgeData, colIdx, mTxt )
{
  isState = lapply( seq( nSts ), function(x)  clEdgeData[ , ST_ID ] == x )
  isHgt    = clEdgeData[ , HGT ] == HGT_EVENT
  isImp    = clEdgeData[ , HGT ] == IMPORTATION
  nProps = c(
   sapply( isState, function(x)
     sum( clEdgeData[ x & isHgt, colIdx ] > 1 ) / sum( x & isHgt ) ),
   sapply( isState, function(x)
     sum( clEdgeData[ x & isImp, colIdx ] > 1 ) / sum( x & isImp ) )
   )
  names( nProps ) = c( paste( stIds, "HGT" ), paste( stIds, "Importation" ) )

  plOrd  = order( names( nProps ) )
  bpCols = c( stCols, stCols )[ plOrd ]
  nProps = nProps[ plOrd ]

  barplot(
    nProps,
    col    = bpCols,
    border = bpCols,
    las    = 2,
    main   = mTxt,
    ylab   = "Number of propagation events"
    )
}

PlotCounts = function( clEdgeData, colIdx, yInc = 25, mTxt )
{
  # Get the positions on the Y-axis: the counts of the imput column
  yPos =  unlist( c(
    lapply( isState, function(x) clEdgeData[ x & isHgt, colIdx ] ),
    lapply( isState, function(x) clEdgeData[ x & isImp, colIdx ] )
    ))
  maxVal = max( yPos )
  yPos   = rescale( yPos, c( 0.05, 0.95 ) )

  # Get the positions on the X-axis: the state and importation status
  nStHgt = sapply( isState, function(x) sum( x & isHgt ) )
  nStImp = sapply( isState, function(x) sum( x & isImp ) )
  xPos = unlist(c(
    lapply( seq(3), function(i) rep( i, nStHgt[i] ) ),
    lapply( seq(3), function(i) rep( i  + 0.5, nStImp[i] ) )
    ))
  xPos = rescale( xPos, c( 0.05, 0.95 ) )

  # Get the colors
  cols = unlist(c(
    lapply( seq(3), function(i) rep( stCols[i], nStHgt[i] ) ),
    lapply( seq(3), function(i) rep( stCols[i], nStImp[i] ) )
    ))

  plot.new()
  points(
    x    = jitter( xPos ),
    y    = yPos, # jitter( yPos ),
    col  = cols,
    pch  = 16,
    xlab = "State",
    ylab = "Isolate count"
    )
  mtext( mTxt )
  axis(
    side   = 1,
    at     = seq( 0.05, 0.95, length.out = 6 ),
    labels = stratVals[ order( stratVals ) ],
    las    = 2
    )
  nTicks = ceiling( maxVal / yInc )
  axis(
    side   = 2,
    at     = seq( 0.05, 0.95, length.out = nTicks + 1 ),
    labels = seq( 0, yInc * nTicks, yInc ),
    las    = 2
    )
}

GetStrCols = function( clEdgeData )
{
  return( sapply( clEdgeData[ , STR_ID ], function(x) strainCols[x] ) )
}

PlotLogCounts = function( clEdgeData, colIdx, yInc = 25, mTxt, baseVal = 2 )
{
  plRange = c( 0.05, 0.95 )

  # Get the positions on the Y-axis: the counts of the imput column
  countVals =  unlist( c(
    lapply( isState, function(x) clEdgeData[ x & isHgt, colIdx ] ),
    lapply( isState, function(x) clEdgeData[ x & isImp, colIdx ] )
    ))
  isCount   = countVals != 0
  countVals = countVals[ isCount ]

  nTicks = ceiling( max( countVals ) / yInc )
  maxVal = yInc * nTicks
  yPos   = log( c( 1, maxVal, countVals ), base = exp( baseVal ) )
  yPos   = rescale( yPos, plRange )
  yPos   = yPos[ 3:length( yPos ) ]


  axisVals = seq( 0, maxVal, yInc )
  axisVals[1] = 1
  yAxisPos  = rescale( log( axisVals, base = exp( baseVal ) ), plRange )
  xAxisPos  = seq( plRange[1], plRange[2], length.out = 6 )

  # Get the positions on the X-axis: the state and importation status
  nStHgt = sapply( isState, function(x) sum( x & isHgt ) )
  nStImp = sapply( isState, function(x) sum( x & isImp ) )
  xPos = unlist(c(
    lapply( seq(3), function(i) rep( i, nStHgt[i] ) ),
    lapply( seq(3), function(i) rep( i  + 0.5, nStImp[i] ) )
    ))
  xPos = xPos[ isCount ]
  xPos = rescale( xPos, plRange )
  # Get the colors
  # cols = unlist(c(
  #   lapply( seq(3), function(i) rep( stCols[i], nStHgt[i] ) ),
  #   lapply( seq(3), function(i) rep( stCols[i], nStImp[i] ) )
  #   ))
  cols = unlist( c(
    lapply( isState, function(x) GetStrCols( clEdgeData[ x & isHgt, ] ) ),
    lapply( isState, function(x) GetStrCols( clEdgeData[ x & isImp, ] ) )
    ))
  cols = cols[ isCount ]

  plot.new()
  points(
    x    = jitter( xPos, 1.5 ),
    y    = jitter( yPos, 1.5 ), # yPos
    col  = cols,
    pch  = 16,
    xlab = "State",
    ylab = "Isolate count"
    )
  title( mTxt )
  axis(
    side   = 1,
    at     = xAxisPos,
    labels = stratVals[ order( stratVals ) ],
    las    = 2
    )
  axis(
    side   = 2,
    at     = yAxisPos,
    labels = axisVals,
    las    = 2
    )
  xIncYal = ( xAxisPos[2] - xAxisPos[1] ) / 2
  segments(
     x0  = xAxisPos[1:5] + xIncYal,
     y0  = rep( plRange[1], 5 ),
     x1  = xAxisPos[1:5] + xIncYal,
     y1  = rep( plRange[2], 5 ),
     col = "#a1a1a1",
     lwd = 2
     )
  # legend(
  #   "topleft",
  #   col    = strainCols,
  #   legend = strainVals,
  #   pch    = 15,
  #   cex    = 0.75
  #   )
}

# ---- Get the data on the facilities ------------------------------------------

# Get the facility of isolation and sort to match the other meta-data
facIds = MakeFacIds( eipData )
eipOrd =
  sapply( row.names( annots ), function(x) which( eipData$Sample_ID == x ) )
facIds = facIds[ eipOrd ]
patIds = eipData$Patient_ID[ eipOrd ]

# Make sure everything is sorted
identical( names( facIds ), row.names( annots ) )
identical( eipData$Sample_ID[ eipOrd ], row.names( annots ) )
identical( row.names( eipAcqEdge), row.names( annots ) )
identical( row.names( stHgtStatus ), row.names( eipAcqEdge ) )
identical( row.names( annots ), row.names( eipAcqEdge ) )

# Get the data on the strains in the study
strainId   = GetStrainId( annots )
strainVals = unique( strainId )
nStrs      = length( strainVals )
strainVals = strainVals[ order( strainVals ) ]
strIdxList = lapply( strainVals, function(x) which( strainId == x ) )
colVals    = c( annotCols$ST, annotCols$species[ 1:5 ] )
strainCols = colVals[ order( names( colVals ) ) ]

# Get the data on the states in the study
stVals    = names( annotCols$sources )[ 1:3 ]
nSts      = length( stVals )
stCols    = annotCols$sources[ 1:3 ]
stIdxList = lapply( stVals, function(x) which( annots$sources == x ) )
nStIsos   = sapply( stIdxList, length )

# Convert the list into a vector
strIdxs = MakeIdxVec( strIdxList )
stIdxs  = MakeIdxVec( stIdxList )

# ---- Make a data-frame with the relevant information on each edge ------------

stRow = lapply( c( "EIP_MN", "EIP_TN", "EIP_CT" ),
  function(x) annots$sources == x
  )


# Parse out the data for each cluter acquisition edge
clEdgeData = data.frame()
for ( j in 1:ncol( eipAcqEdge ) )
{
  for ( st in stRow )
  {
    edgeVals    = unique( eipAcqEdge[ st, j ] )
    edgeVals    = edgeVals[ edgeVals != 0 ]
    if ( length( edgeVals ) > 0 )
    {
      edgeIdxs    =
        lapply( edgeVals, function(x) which( eipAcqEdge[ , j ] == x & st ) )
      edgeIdxs    =
        lapply( edgeIdxs, function(x) x[ !duplicated( patIds[x] ) ] )
      edgeReps    = sapply( edgeIdxs, function(x) x[1] )
      edgeStVals  = stIdxs[ edgeReps ]
      edgeStrVals = strIdxs[ edgeReps ]
      hgtStats    = stHgtStatus[ edgeReps, j ]
      nEdge       = sapply( edgeIdxs, length )
      edgeFacs    =
        sapply( edgeIdxs, function(x) unique( facIds[x][ facIds[x] != '' ] ) )
      nFacs       = sapply( edgeFacs, length )

      stEdgeData = cbind.data.frame(
        "clId"     = rep( j, length( edgeVals ) ),
        "stId"     = edgeStVals,
        "strId"    = edgeStrVals,
        "nEdge"    = nEdge,
        "nFacs"    = nFacs,
        "hgt"      = hgtStats,
        "edgeVals" = edgeVals
        )
      clEdgeData = rbind.data.frame( clEdgeData, stEdgeData )
    }
  }
}

# ---- Analyze the HGT cluster distribution ------------------------------------

isState  = lapply( seq( nSts ), function(x)  clEdgeData[ , ST_ID ] == x )
nStEdges = sapply( isState, sum )

# Find which edges correspond to HGT, and importation. And which correspond to
# the respective states and strains
isHgt    = clEdgeData[ , HGT ] == HGT_EVENT
isImp    = clEdgeData[ , HGT ] == IMPORTATION
evCounts = c( "HGT" = sum( isHgt ), "Importation" = sum( isImp ) )

# Calculate the Shannon diversity
shDiv = c(
  CalcShannonDiv( clEdgeData[ isHgt, ] ), CalcShannonDiv( clEdgeData[ isImp, ] )
  )
stIds = c( "CT", "MN", "TN" )
stratVals = c( paste( stIds, "HGT" ), paste( stIds, "Importation" ) )
names( shDiv ) = stratVals

#
wilResFacHgt = wilcox.test( nFacs ~ hgt, data = clEdgeData, exact = FALSE )
wilResIsoHgt = wilcox.test( nEdge ~ hgt, data = clEdgeData, exact = FALSE )
wilResFacHgt$p.value


stPvals = sapply( isState, function(x)
{
  wilResFacHgt =
    wilcox.test( nFacs ~ hgt, data = clEdgeData[ x, ], exact = FALSE )
  wilResIsoHgt =
    wilcox.test( nEdge ~ hgt, data = clEdgeData[ x, ], exact = FALSE )
  pVals =
    c( "facHgt" = wilResFacHgt$p.value, "isoHgt" = wilResIsoHgt$p.value )
  return( pVals )
})
colnames( stPvals )  = stVals
row.names( stPvals ) = c( "Facility", "Number of isolates" )
stPvals

GetStrData = function( clEdgeData )
{
  strCounts = sapply( seq( length( strIdxList ) ),
    function(x) sum( clEdgeData[ , STR_ID ] == x )
    )
  names( strCounts ) = strainVals
  strFracs = sapply( strCounts, function(x) x / nrow( clEdgeData ) )
  return( cbind.data.frame( strCounts, strFracs ) )
}

hgtStrData = GetStrData( clEdgeData[ isHgt, ] )
impStrData = GetStrData( clEdgeData[ isImp, ] )

ki = grep( "K. pneumoniae", row.names( hgtStrData ) )

sum( hgtStrData[ ki, 1 ] ) / sum( hgtStrData[ , 1 ] )
sum( impStrData[ ki, 1 ] ) / sum( impStrData[ , 1 ] )
# [1] 0.7459677
# > sum( impStrData[ ki, 1 ] ) / sum( impStrData[ , 1 ] )
# [1] 0.8962264

dim( annots )
sum( annots$species == "K. pneumoniae" )


clStrMat = sapply( unique( clEdgeData[ , CL_ID ] ), function(y)
{
  isCl = clEdgeData[ ,CL_ID ] == y
  strClCounts = sapply( seq( length( strIdxList ) ),
    function(x) sum( clEdgeData[ isCl, STR_ID ] == x )
    )
  return( strClCounts )
})
row.names( clStrMat ) = strainVals

ki = grep( "K. pneumoniae", strainVals )
xx = sapply( 1:nrow( clStrMat ), function(i) sum( clStrMat[ i, ] > 0 ) ) /
  ncol( clStrMat )
names( xx ) = row.names( clStrMat )
pheatmap( clStrMat )


testCols = colorRampPalette(
  colors      = c( "#395983", "#e7e7e7", "#952827" ),
  interpolate = "linear",
  bias = 2.5
  )( 100 )
testCols = rev( testCols )
barplot(
  seq( 100 ),
  col    = testCols,
  border = testCols
  )

# ---- Save the data -----------------------------------------------------------

print( pryr::mem_used() )
rData = "../data/2022_03_09_find_forward_hgt_cluster_transmission.rData"
save( file = rData, list = ls() )

# ---- Make the figures --------------------------------------------------------

figDir = "../figures/2022_03_09_find_forward_hgt_cluster_transmission/"
figPaths = paste0( figDir, "2022_03_09_hgt_imp_counts.png" )

clCounts = c(
  "HGT"         = sum( clEdgeData[ isHgt, N_EDGE ] ),
  "Importation" = sum( clEdgeData[ isImp, N_EDGE ] )
  )
png(
  filename = figPaths[1],
  units    = "in",
  width    = 6.5,
  height   = 8,
  res      = 300
  )
par( mar = c( 5.1, 4.1, 4.1, 2.1 ), mfrow = c( 1, 2 ) )
barplot(
  evCounts,
  col    = "#395983",
  border = "#395983",
  ylim   = c( 0, 800 ),
  main   = "Counts of cluster\nacquisition events"
  )
barplot(
  clCounts,
  col    = "#395983",
  border = "#395983",
  ylim   = c( 0, 800 ),
  main   = "Counts of clusters\ncarried"
  )
dev.off()


png(
  filename = figPaths[1],
  units    = "in",
  width    = 6.5,
  height   = 8,
  res      = 300
  )
par( mar = c( 5.1, 4.1, 4.1, 2.1 ) )
barplot(
  evCounts,
  col    = "#395983",
  border = "#395983",
  ylim   = c( 0, 250 )
  )
dev.off()

figPaths[2] = paste0( figDir, "2022_03_09_shannon_diversity.png" )
png(
  filename = figPaths[2],
  units    = "in",
  width    = 6.5,
  height   = 8,
  res      = 300
  )
par( mar = c( 7.5, 4.1, 4.1, 2.1 ) )
shOrd  = order( names( shDiv ) )
bpCols = c( stCols, stCols )[ shOrd ]
barplot(
  shDiv[ shOrd ],
  col    = bpCols,
  border = bpCols,
  ylim   = c( 0, 2 ),
  las    = 2,
  main   = "Shannon diversity by state"
  )
dev.off()

figPaths[3] = paste0( figDir, "2022_03_09_isolates_by_acquisition.png" )
png(
  filename = figPaths[3],
  units    = "in",
  width    = 6.5,
  height   = 8,
  res      = 300
  )
par( mar = c( 7.5, 4.1, 4.1, 2.1 ) )
PlotCounts( clEdgeData, N_EDGE, 25, "Number of isolates per acquisition" )
dev.off()

figPaths[4] = paste0( figDir, "2022_03_09_facilities_by_acquisition.png" )
png(
  filename = figPaths[4],
  units    = "in",
  width    = 6.5,
  height   = 8,
  res      = 300
  )
par( mar = c( 7.5, 4.1, 4.1, 2.1 ) )
PlotCounts( clEdgeData, N_FACS, 10, "Number of isolates per acquisition" )
dev.off()

figPaths[5] = paste0( figDir, "2022_03_09_fac_propagation_per_acq.png" )
png(
  filename = figPaths[5],
  units    = "in",
  width    = 6.5,
  height   = 8,
  res      = 300
  )
par( mar = c( 7.5, 4.1, 4.1, 2.1 ) )
plTitle = "Fraction of isolate propagation events per acquisition"
PlotNumProps( clEdgeData, N_EDGE, plTitle )
dev.off()


figPaths[6] = paste0( figDir, "2022_03_09_fac_propagation_per_acq.png" )
png(
  filename = figPaths[6],
  units    = "in",
  width    = 6.5,
  height   = 8,
  res      = 300
  )
par( mar = c( 7.5, 4.1, 4.1, 2.1 ) )
plTitle = "Fraction of facility propagation events per acquisition"
PlotNumProps( clEdgeData, N_FACS, plTitle )
dev.off()


figPaths[7] = paste0( figDir, "2022_03_09_iso_counts_per_acq_log.png" )
png(
  filename = figPaths[7],
  units    = "in",
  width    = 6.5,
  height   = 8,
  res      = 300
  )
par( mar = c( 7.5, 4.1, 4.1, 2.1 ) )
plTitle = "Isolates per acquisition"
PlotLogCounts( clEdgeData, N_EDGE, yInc = 25, plTitle, baseVal = 2 )
dev.off()

figPaths[8] = paste0( figDir, "2022_03_09_fac_counts_per_acq_log.png" )
png(
  filename = figPaths[8],
  units    = "in",
  width    = 6.5,
  height   = 8,
  res      = 300
  )
par( mar = c( 7.5, 4.1, 4.1, 2.1 ) )
plTitle = "Facilities per acquisition"
PlotLogCounts( clEdgeData, N_FACS, yInc = 10, plTitle, baseVal = 2 )
dev.off()



figPaths[9] = paste0( figDir, "2022_03_09_fac_counts_per_acq_log.png" )
hmCols = colorRampPalette( c( "#395983", "#e7e7e7", "#952827" ) )( 100 )
pheatmap(
  stPvals,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  color             = hmCols,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  display_numbers   = TRUE,
  number_format     = "%.3e", #"%.3f",
  number_color      = "white",
  filename          = figPaths[9]
  )

# ---- Generate the Rmarkdown --------------------------------------------------

sink( "2022_03_09_find_forward_hgt_cluster_transmission.Rmd" )

cat(
"---
title: \"\"
author: \"Ryan D. Crawford\"
date: \"03/09/2022\"
output: pdf_document
---\n\n\n"
)

cat( paste0( "![](", figPaths, ")\n\n", collapse = '' ) )

sink()

# ------------------------------------------------------------------------------