# ------------------------------------------------------------------------------
#
# 2022/02/25
# Ryan D. Crawford
# ------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( ape )

# ---- Constant declarations ---------------------------------------------------

YEAR  = 1
MONTH = 2
DAY   = 3

HGT_EVENT   = 1
IMPORTATION = 2


STRAIN   = 1
STATE    = 2
FACILITY = 3

# ---- Import the data ---------------------------------------------------------

# Load the annotations and colors: "annots" and "annotCols"
load( "../data/2022_02_07_make_annotations_and_colors.rData" )

# Load the data on the plasmid clusters: "eipClMat", "eipAcqEdge", "plMbrList",
# and "plMat"
load( "../data/2022_02_16_analyze_hgt_cluster_distribution.rData" )

# Load the meta-data from the EIP isolates
load( "../data/eipData.rData" )

# Load the data on whether cluster were acquired prior to importation to state
# or facility: "stHgtStatus", and "facHgtStatus"
load( "../data/2022_02_21_find_cluster_hgt_status.rData" )

# Load the mpr results on the edge in which an isolate was included in a region.
# includes: "arFeatureMats", "edgeColMat", "featEdgeMat", "edgeFeatList",
# and "tipTrEdgeMat"
load( "../data/2022_02_19_run_mpr_for_location.rData" )

# "geneEdgeMat", "clTrData", "clMat", "mbrList", "clAcqMat",
# "clAcqEdgeList", and "clDelEdgeList"
load( "../data/2022_02_08_get_cl_gene_counts_by_edge.rData" )

# Read in the core-gene tree
cgTree = read.tree( "../data/cgTree.tre" )

# ---- Function definitions ----------------------------------------------------

source( "R/PlotTree.R" )

source( "R/PlotTrace.R" )

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
  return( facilityId )
}

GetStrainId = function( annots )
{
  strainId = sapply( 1:nrow( annots ), function(i)
  {
    if ( annots$ST[i] != "Other" ) return( annots$ST[i] )
    return( annots$Species[i] )
  })
  return( strainId )
}

FixDate = function( inDate )
{
  if ( is.na( inDate ) ) return( NA )
  if ( grepl( '-', inDate, fixed = TRUE ) )
  {
    dateVals = strsplit( inDate, '-' )[[1]]
  } else {
    dateVals = strsplit( inDate, '/' )[[1]]
    # if ( length( dateVals ) != 3 ) return( NA )
    if ( length( dateVals ) == 3 ) dateVals = dateVals[ c( 3, 1, 2 ) ]
  }
  dateVals = as.integer( dateVals )
  if ( length( dateVals ) == 1 ) return( c( dateVals, NA, NA ) )
  if ( length( dateVals ) == 2 ) stop( inDate )
  return( dateVals )
}

GetDateRange = function( dateVals, yearInc, monthInc, dayInc )
{
  if ( !missing( yearInc ) )
  {
    yearRange =
      seq( min( dateVals[ , YEAR ] ), max( dateVals[ , YEAR ] ), yearInc )
    return( yearRange )
  }

  yearRange = seq( min( dateVals[ , YEAR ] ), max( dateVals[ , YEAR ] ))

  if ( !missing( monthInc ) )
  {
    dateRng = unlist( lapply( yearRange,
      function(x) paste0( x, '/', seq( 1, 12, monthInc ) )
      ))
    return( dateRng )
  }
  daysPerMonth =
    c( 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30 , 31 ) - dayInc + 1
  dateRng = unlist( lapply( yearRange,
    function(x) sapply( seq( 1, 12 ),
      function(y) paste0( x, '/', y, '/', seq( 1, daysPerMonth[y] , dayInc ) )
      )
    ))

  return( dateRng )
}

RoundDateMatrix = function( dateMat, yearInc, monthInc, dayInc )
{
  # Round to the nearest year
  if ( yearInc > 1 )
  {
    yearVals =
      RoundDateVals( dateMat[ , YEAR ], yearInc, min( dateMat[ , YEAR ] ) )
    return( as.matrix( yearVals, ncol = 1 ) )
  }
  yearVals = dateMat[ , YEAR ]


  if ( missing( monthInc ) ) return( as.character( yearVals ) )

  # Round to the yearest month
  if ( monthInc > 1 )
  {
    monthVals = RoundDateVals( dateMat[ , MONTH ], monthInc, 1 )

    dateVals = t( sapply( 1:nrow( dateMat ),
      function(i) c( yearVals[i], monthVals[i] )
      ))
    return( dateVals )
  }
  monthVals = dateMat[ , MONTH ]

  # Round to the nearest day
  dayVals = RoundDateVals( dateMat[ , DAY ], dayInc, 1 )

  for ( i in 1:length( dateVals ) )
  {
    if ( dateVals[i] == 29 )
    {
      dateVals[i] = 1
      if ( monthVals[i] == 12 )
      {
        monthVals[i] = 1
      } else {
        monthVals[i] = monthVals[i] + 1
      }
    }
  }

  dateVals = t(sapply( 1:nrow( dateMat ),
    function(i) c( yearVals[i], monthVals[i], dayVals[i]  )
    ))
  return( dateVals )
}

RoundDateVals = function( dateVals, inc, minVal )
{
  return( minVal + floor( ( dateVals - minVal ) / inc ) * inc )
}

MakeDateStrs = function( dateMat )
{
  dateStrVals = sapply( 1:nrow( dateMat ), function(i)
  {
    if ( TRUE %in% is.na( dateMat[ i, ] ) ) return( NA )
    return( paste0( dateMat[ i, ], collapse = '/' ) )
  })
  return( dateStrVals )
}

ConStatusToStr = function( facHgtStatus )
{
  hgtStats = sapply( 1:ncol( facHgtStatus ), function(j)
  {
    sapply( facHgtStatus[ , j ], function(x)
    {
      if ( is.na(x) ) return( NA )
      if ( x == HGT_EVENT ) return( "HGT" )
      if ( x == IMPORTATION ) return( "Importation" )
      return( NA )
    })
  })

  return( hgtStats )
}

# ---- Merge the data on the plasmids with the meta-data -----------------------

# Subset to the isolates who have a facility and sort the data.
eipData   = eipData[ !is.na( eipData$Fac_ID ), ]
eipData   = eipData[ eipData$Sample_ID %in% row.names( annots ), ]
eipData   = eipData[ order( eipData$Sample_ID ), ]
eipAnnots = annots[ row.names( annots ) %in% eipData$Sample_ID, ]
eipAnnots = eipAnnots[ order( row.names( eipAnnots ) ), ]

# Parse the meta-data on facility, strain, and date to create the meta-data
# to plot
facIds   = MakeFacIds( eipData )
strainId = GetStrainId( eipAnnots )
dateMat  = t( sapply( eipData$Cult_Date, FixDate ) )

# Make a data-frame with the parsed meta-data
plData   = cbind.data.frame( facIds, strainId )
row.names( plData ) = row.names( eipAnnots )
colnames( plData )  = c( "Facility", "Strain" )

# Sort the data on the cluster HGT status
facHgtStatus =
  facHgtStatus[ row.names( facHgtStatus ) %in% row.names( plData ), ]
facHgtStatus = facHgtStatus[ order( row.names( facHgtStatus ) ), ]

stHgtStatus =
  stHgtStatus[ row.names( stHgtStatus ) %in% row.names( plData ), ]
stHgtStatus = stHgtStatus[ order( row.names( stHgtStatus ) ), ]

# Sort the data on the presence of the gene clusters in
eipClMat = eipClMat[ row.names( eipClMat ) %in% row.names( plData ), ]
eipClMat = eipClMat[ order( row.names( eipClMat ) ), ]
eipClMat = eipClMat[ , colnames( eipClMat ) %in% colnames( facHgtStatus ) ]

# Make a list of the strain colors
colVals  = c( annotCols$ST, annotCols$Species[ 1:5 ] )
# Yellow: "#f0bc5e", "#eaa221"; Dark gray: "#212e41"; Light gray: "#dde1ea"
hgtCol = "#eaa221"
impCol = "#dde1ea"
absCol = "#212e41"

hgtColList = list( "hgtStatus" = c( "HGT" = hgtCol, "Importation" = impCol ) )

plotCols = c(
  list( "Strain" = colVals[ order( names( colVals ) )] ), hgtColList
  )


hgtStatus = ConStatusToStr( facHgtStatus )

# Make an edge x gene cluster matrix with the number of cluster genes
# present on each edge
clEdgeMat = sapply( names( plMbrList ), function(y)
{
  isCl = sapply( mbrList[ as.integer( strsplit( y, '/' )[[1]] ) ],
    function(x) rowSums( geneEdgeMat[ , x ] ) >= ceiling( length(x) * 0.51 )
    )
  return( sapply( 1:nrow( isCl ), function(i)  TRUE %in% isCl[ i, ] ) )
})

clEdgeMat = sapply( plMbrList,
  function(x) rowSums( geneEdgeMat[ , x ] ) >= ceiling( length(x) * 0.51 )
  )

# ---- Make the figures --------------------------------------------------------

# Determine which plots to make
clCounts   = colSums( eipClMat )
countOrder = order( clCounts, decreasing = TRUE )
plotIdxs   = countOrder[ 1:32 ]


figDir = "../figures/2022_02_25_plot_plasmid_trace/"

# Make the target plots
targetPlPaths = sapply( plotIdxs, function(j)
{
  isCl = eipClMat[ , j ] == 1

  clDateMat = RoundDateMatrix( dateMat[ isCl, ], 1, 3 )
  dateRange = GetDateRange( clDateMat, monthInc = 3 )
  clPlData  = cbind.data.frame(
    MakeDateStrs( clDateMat ), plData[ isCl, ], hgtStatus[ isCl, j ]
    )
  clPlData  = clPlData[ complete.cases( clPlData ), ]

  if (
    length( unique( clPlData[ , 1 ] ) ) == 1 &
    length( unique( clPlData[ , 2 ] ) ) == 1
    )
  {
    stop( j, ": This one is not good..." )
  }

  GetStHgtStats( clPlData, "MN" )
  GetStHgtStats( clPlData, "TN" )


  figPath = paste0( figDir, "2022_02_25_", j, "_target_plot.png" )
  png(
    filename = figPath,
    units    = "in",
    width    = 8,
    height   = 9,
    res      = 300
    )
  par( mar = c( 5.1, 4.1, 4.1, 2.1 ) )
  TargetPlot(
    plData     = clPlData,
    xColIdx    = 1,
    xAxisVals  = dateRange,
    yColIdx    = 2, # yAxisLabs  = yAxisLabs,
    catColIdxs = c( 3, 4 ),
    plotCols   = plotCols,
    logScale   = TRUE,
    logBase    = exp(2), # baseCexVal = 1.5,
    cexAdjVal  = 0.5,
    baseCexVal = 2
    )
  # mtext( paste0( "Gene cluster ", j, "\n" ) )
  dev.off()

  return( figPath )
})


GetStHgtStats = function( clPlData, stId )
{
  isTn    = grepl( stId, row.names( clPlData ) )
  isHgt   = clPlData$`hgtStatus[isCl, j]` == "HGT"
  isTnHgt = isTn & isHgt

  nImpFacs =
    sum( clPlData$Facility[ isTnHgt ] %in% clPlData$Facility[ isTn & !isHgt ])

  nHgt = sum( isTnHgt )
  fracFacOvlp = nImpFacs/ nHgt

  cat(
    "Number of hgt events: ", nHgt, '\n',
    "Fraction with facility overlap: ", fracFacOvlp, '\n',
    sep = ''
    )

}

# Create a vector to rainbow colors to assign to the facilities
cols = c( "#ba3231", "#f37a48", "#eaa221", "#00bc48", "#1697B7", "#7f4485")

# Find the edges connected to the tips
isTip = cgTree$edge[ , 2 ] <= length( cgTree$tip.label )

cat( "cgTreeFacPaths\n" )

# Find the corresponding row in the HGT hgt data
annotHgtIdxs = sapply( row.names( annots ), function(x)
{
  isG = which( row.names( facHgtStatus ) == x )
  if ( length( isG ) ) return( isG )
  return( NA )
})

identical(
  row.names( facHgtStatus )[ annotHgtIdxs[ !is.na( annotHgtIdxs ) ] ],
  row.names( annots )[ !is.na( annotHgtIdxs ) ]
  )

hgtColList = list( "hgtStatus" = c(
  "HGT" = hgtCol, "Importation" = impCol, "Absent" = "#212e41"
  ))
annotCols = c( annotCols, hgtColList )

cgTreeFacPaths = sapply( plotIdxs, function(j)
{
  #
  isClTipEdge = clEdgeMat[ , j ] == 1 & isTip

  # Find the facilities with this cluster
  facVals = unique( featEdgeMat[ isClTipEdge, FACILITY ] )
  facVals = facVals[ facVals != 0 ]
  facCols = colorRampPalette( cols )( length( facVals ) )

  clHgtStatus = sapply( annotHgtIdxs, function(i)
  {
    if ( is.na(i) ) return( "Absent" )
    if ( is.na( hgtStatus[ i, j ] ) ) return( "Absent" )
    return( hgtStatus[ i, j ] )
  })
  clPlData = cbind.data.frame( annots, clHgtStatus )
  # clPlData = clPlData[ complete.cases( clPlData ), ]

  edgeCols = sapply( 1:nrow( cgTree$edge ), function(i)
  {
    isFac = featEdgeMat[ i, FACILITY ] != 0
    isCl  = clEdgeMat[ i, j ] == 1
    if ( isFac & isCl )
      return( facCols[ facVals == featEdgeMat[ i, FACILITY ] ] )
    if ( isCl ) return( "#767676" )
    return( "#a1a1a1" )
  })

  edgeWidth = ifelse( edgeCols ==  "#a1a1a1", 0.75, 2 )

  par( mar = c( 5.1, 4.1, 4.1, 2.1 ) )
  figPath = paste0( figDir, "2022_02_25_", j, "_tree_plot.png" )
  png(
    filename = figPath,
    units    = "in",
    width    = 6.5,
    height   = 8,
    res      = 300
    )
  PlotTree( cgTree, clPlData, annotCols, edgeCols, edgeWidth )
  mtext( paste0( "Gene cluster ", j, "\n" ) )
  dev.off()

  return( figPath )
})
cat( " -- Done\n" )


# cat( "clCgTreeFacPaths\n" )
# clCgTreeFacPaths = sapply( plotIdxs, function(j)
# {
#   #
#   isClTipEdge = clEdgeMat[ , j ] == 1 & isTip
#
#   # Find the facilities with this cluster
#   facVals = unique( featEdgeMat[ isClTipEdge, FACILITY ] )
#   facVals = facVals[ facVals != 0 ]
#   facCols = colorRampPalette( cols )( length( facVals ) )
#
#
#   clPlData = cbind.data.frame( eipAnnots, facHgtStatus[ , j ]  )
#   clPlData = clPlData[ complete.cases( clPlData ), ]
#
#   edgeCols = sapply( 1:nrow( cgTree$edge ), function(i)
#   {
#     isFac = featEdgeMat[ i, FACILITY ] != 0
#     isCl  = clEdgeMat[ i, j ] == 1
#     if ( isFac & isCl )
#       return( facCols[ facVals == featEdgeMat[ i, FACILITY ] ] )
#     if ( isCl ) return( "#767676" )
#     return( "#a1a1a1" )
#   })
#
#   edgeWidth = ifelse( edgeCols ==  "#a1a1a1", 0.75, 2 )
#
#
#
#   par( mar = c( 5.1, 4.1, 4.1, 2.1 ) )
#   figPath = paste0( figDir, "2022_02_25_", j, "_tree_plot.png" )
#   png(
#     filename = figPath,
#     units    = "in",
#     width    = 6.5,
#     height   = 8,
#     res      = 300
#     )
#   PlotTree( cgTree, clPlData, annotCols, edgeCols, edgeWidth )
#   mtext( paste0( "Gene cluster ", j, "\n" ) )
#   dev.off()
#
#   return( figPath )
# })
# cat( " -- Done\n" )

row.names( cgTree$edge ) = seq( nrow( cgTree$edge ) )

facIdx = sapply( 1:nrow( arFeatureMats[[3]] ), function(i)
{
  ii = which( arFeatureMats[[3]][ i, ] == 1 )
  if ( length( ii ) == 0 ) return( NA )
  return( ii )
})
names( facIdx ) = row.names( arFeatureMats[[3]] )
cat( "clCgTreeFacPaths\n" )
clCgTreeFacPaths = sapply( plotIdxs, function(j)
{
  #
  isClTipEdge = clEdgeMat[ , j ] == 1 & isTip

  # Find the facilities with this cluster
  facVals = unique( featEdgeMat[ isClTipEdge, FACILITY ] )
  facVals = facVals[ facVals != 0 ]

  facCols  = colorRampPalette( cols )( length( facVals ) )
  clPlData = cbind.data.frame( eipAnnots,  hgtStatus[ , j ] )
  clPlData = clPlData[ complete.cases( clPlData ), ]


  clTree = keep.tip( cgTree, row.names( clPlData ) )

  edgeCols = sapply( 1:nrow( clTree$edge ), function(i)
  {
    ii    = as.integer( row.names( clTree$edge ) )[i]
    isFac = featEdgeMat[ ii, FACILITY ] != 0
    isCl  = clEdgeMat[ ii, j ] == 1
    if ( isFac & isCl )
      return( facCols[ facVals == featEdgeMat[ ii, FACILITY ] ] )
    if ( isCl ) return( "#767676" )
    if ( clTree$edge[ i, 2 ] <= length( clTree$tip.label ) )
    {
      isIso = names( facIdx ) == clTree$tip.label[ clTree$edge[ i, 2 ]  ]
      isFac = facVals == facIdx[ isIso ]
      if ( TRUE %in% isFac) return( facCols[ isFac ] )
      return( "#767676" )
    }
    return( "#a1a1a1" )
  })


  edgeWidth = ifelse( edgeCols ==  "#a1a1a1", 0.75, 2 )

  par( mar = c( 5.1, 4.1, 4.1, 2.1 ) )
  figPath = paste0( figDir, "2022_02_25_", j, "_cl_tree_plot.png" )
  png(
    filename = figPath,
    units    = "in",
    width    = 6.5,
    height   = 8,
    res      = 300
    )
  PlotTree( clTree, clPlData, annotCols, edgeCols, edgeWidth )
  mtext( paste0( "Gene cluster ", j, "\n" ) )
  dev.off()
  return( figPath )
})
cat( " -- Done\n" )

# ---- Generate the Rmarkdown --------------------------------------------------

sink( "2022_02_25_plot_plasmid_trace.Rmd" )

cat(
"---
title: \"Plot plasmid trace\"
author: \"Ryan D. Crawford\"
date: \"02/25/2022\"
output: pdf_document
---\n\n\n"
)

for ( i in 1:length( plotIdxs ) )
{
  cat( "## Gene cluster", plotIdxs[i], "\n\n" )
  cat( "![](", targetPlPaths[i], ")\n\n", sep = '' )
  cat( "![](", cgTreeFacPaths[i], ")\n\n", sep = '' )
  cat( "![](", clCgTreeFacPaths[i], ")\n\n", sep = '' )
}
sink()

# ------------------------------------------------------------------------------
# dateMat[ isCl, ]
# dateRange = GetDateRange( dateMat[ isCl, ], yearInc = 1  )
#
# clPlData = cbind.data.frame( dateMat[ isCl, YEAR ], plData[ isCl, ] )
#
# TargetPlot(
#   plData     = clPlData,
#   xColIdx    = 1,
#   xAxisVals  = dateRange,
#   yColIdx    = 2, # yAxisLabs  = yAxisLabs,
#   catColIdxs = 3,
#   plotCols   = plotCols,
#   logScale   = TRUE,
#   logBase    = exp(1), # baseCexVal = 1.5,
#   cexAdjVal  = 0.75
#   )
# plData     = clPlData;
# xColIdx    = 1;
# xAxisVals  = dateRange;
# yColIdx    = 2; # yAxisLabs  = yAxisLabs;
# catColIdxs = c( 3, 4 );
# plotCols   = plotCols;
# logScale   = TRUE;
# logBase    = exp(1); # baseCexVal = 1.5;
# cexAdjVal  = 0.5;
# baseCexVal = 2;
# x <- stats::runif(12); y <- stats::rnorm(12)
# i <- order(x, y); x <- x[i]; y <- y[i]
# plot(x,y, main = "arrows(.) and segments(.)")
# ## draw arrows from point to point :
# s <- seq(length(x)-1)  # one shorter than data
# arrows(x[s], y[s], x[s+1], y[s+1], col = 1:3)
# s <- s[-length(s)]
# segments(x[s], y[s], x[s+2], y[s+2], col = "pink")