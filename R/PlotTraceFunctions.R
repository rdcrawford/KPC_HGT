PlotGrid = function( xPoints, yPoints )
{
  clipVal = 0.01
  clip( 
    xPoints[1] - clipVal,
    xPoints[ length( xPoints ) ] + clipVal,
    yPoints[1] - clipVal,
    yPoints[ length( yPoints ) ] + clipVal
    )
  
  for ( x in yPoints ) abline( h = x, col = "grey" )
  for ( x in xPoints ) abline( v = x, col = "grey" )
}

GetAxisPoints = function( n )
{
  return( seq( n ) / n )
}

PlotAxis = function( pointVals, labelVals, axisVal, offSetVal )
{
  if ( axisVal == 'x' )
  {
    sideVal = 1
  } else if ( axisVal == 'y' ) {
    sideVal = 2
  } else {
    stop( "axis must be x or y\n" )
  }
  
  if ( missing( offSetVal ) ) offSetVal = 0
  
  # Add a legend for the number of plasmids
  axis(
    side   = sideVal, 
    at     = pointVals, 
    labels = labelVals,
    las    = 2,
    pos    = ,
    outer  = FALSE
    )
}

PlotTrace = function( 
  isolateMetaData,
  dateColIdx, 
  yColIdx,
  catColIdxs,
  plotCols, 
  logScale, 
  xOffSet, 
  yOffSet 
  )
{
  # Set the values to shift the x and y axis
  if ( missing( xOffSet ) )    xOffSet    = 0
  if ( missing( yOffSet ) )    yOffSet    = 0.05
  if ( missing( catColIdxs ) ) catColIdxs = NULL
  if ( missing( logScale ) )   logScale   = FALSE
  

  nDates     = length( plDateVals )
  
  # Get the points for the y-axis
  yAxisVals  = unique( isolateMetaData[ , yColIdx ] )
  nYaxisVals = length( yAxisVals )
  
  # Get the points to plot along the x and y axis
  xPoints = GetAxisPoints( nDates ) - xOffSet
  yPoints = GetAxisPoints( nYaxisVals ) - yOffSet
  
  catVals = unique( isolateMetaData[ , catColIdxs ] )
  catVals = catVals[ order( catVals, decreasing = TRUE ) ]
  catVals = names( plotCols )
  
  xIncVal = xPoints[2] - xPoints[1]
  # Make a list with the coordinates for each point on the plot 
  plotCoords = GetPointCoords( 
    isolateMetaData, dateColIdx, yColIdx, catColIdxs, plotCols, logScale, 
    xPoints, yPoints, plDateVals, yAxisVals, catVals, xIncVal
    )
  
  # Make a list with the coordinates for each point on the plot 
  isEip = grepl( "EIP", isolateMetaData$sources )
  plotSoCoords = GetPointCoords( 
    isolateMetaData[ isEip, ], 
    dateColIdx, yColIdx, 4, c( "red4",  "darkblue", "purple" ) , logScale, 
    xPoints, yPoints, plDateVals, yAxisVals, 
    c( "EIP_CT", "EIP_MN", "EIP_TN" ), xIncVal / 3
    )
  
  
  # Get the size of the dots to plot for the legend
  maxVal = max( sapply( plotCoords, function(x) x$cexVal ) )
  cexVals = GetCexVals( maxVal, logScale )
  
  yAxisLabs = sapply( yAxisVals, function(x) 
  {
    isVal = isolateMetaData[ , yColIdx ] == x
    idx = 
      which( eipData$Sample_ID %in% isolateMetaData$genomeNames[ isVal ] )[1]
    if ( is.na( idx ) ) return( "" )
    if ( eipData$ST[ idx ] != "Other species" )  return( eipData$ST[ idx ] )
    return( isolateMetaData$species[ which( isVal )[1] ] )
  })
  
  plot.new()                                        # Create a new plot
  PlotAxis( xPoints, plDateVals, 'x', yOffSet )     # Plot the x-axis
  PlotAxis( yPoints, yAxisLabs, 'y', xOffSet )      # Plot the y-axis
  PlotLegend( cexVals, xPoints, yOffSet, plotCols ) # Plot the legend points
  PlotPoints( plotCoords )                          # Plot the points 
  # PlotPoints( plotSoCoords )                        # Plot the source points
  PlotGrid( xPoints, yPoints )                      # Plot the xy grid
  PlotPoints( plotCoords )                          # Re-plot the points
  # PlotPoints( plotSoCoords )                        # Re-plot the source points
  PlotLegend( cexVals, xPoints, yOffSet, plotCols ) # Re-plot the legend points
}

GetCexVals = function( valCounts, logScale )
{
  
  maxVal = max( valCounts )
  if ( logScale )
  {
    cexVals = seq( ceiling( maxVal ) )
    names( cexVals ) = sapply( cexVals, function(x) round( exp( x ), 0 ) )
    
  } else {
    cexVals = seq( maxVal )
    names( cexVals ) = cexVals
  }
  
  return( cexVals )
}

PlotLegend = function( cexVals, xPoints, yOffSet, plotCols )
{
  names( plotCols ) = gsub( "United States", "US",  names( plotCols ))
  legVals   = c( names( cexVals ), gsub( ' ', '\n', names( plotCols ) ) )
  legCols   = c( rep( "black", length( cexVals ) ), plotCols )
  colCexVal = cexVals[ ceiling( length( cexVals ) / 2 ) ]
  cexVals   = c( cexVals, rep( colCexVal, length( plotCols ) ) )
  
  points(
    x   = xPoints[ 1:length( cexVals ) ], 
    y   = rep( 1, length( cexVals ) ), 
    cex = cexVals,
    col = legCols,
    pch = 16
    )
  text(
    legVals,
    x   = xPoints[ 1:length( cexVals ) ], 
    y   = 1 - yOffSet, 
    col = "black", 
    cex = 0.75
    )
  
}

GetPointCoords = function( 
  isolateMetaData, dateColIdx, yColIdx, catColIdxs, plotCols, logScale, 
  xPoints, yPoints, plDateVals, yAxisVals, catVals, xIncVal
  )
{
  # Create a vector of all of the 
  xyVals = sapply( 1:nrow( isolateMetaData ), 
    function(i) 
      paste( 
        isolateMetaData[ i, c( yColIdx, dateColIdx ) ], 
        collapse = ' ' 
        )
    )
  
  # 
  plVals = sapply( 1:nrow( isolateMetaData ), 
    function(i) paste( xyVals[i], isolateMetaData[ i, catColIdxs ] )
    )
  
  # 
  xyCats = unique( unique( xyVals ) )
  
  # 
  subPoints = sapply( xyCats, function(x)
  {
    isXy = xyVals == x
    plCats = unique( plVals[ isXy ] )
    if ( length( plCats ) == 1 ) return( list( plCats ) ) 
    return( plCats = plCats[ order( plCats, decreasing = TRUE ) ] )
  })
  nSubPoints = sapply( subPoints, length )
  

  valCounts = table( plVals ) 
  if ( logScale ) valCounts = log( valCounts + 1 )  

  # 
  
  
  # 
  plotCoords = vector( "list", 0 )
  
  # 
  for ( i in 1:length( xyCats ) )
  {
    rIdx = which( xyVals == xyCats[i] )[1]
    xIdx = which( plDateVals == isolateMetaData[ rIdx, dateColIdx ] )
    yIdx = which( yAxisVals == isolateMetaData[ rIdx, yColIdx ] )
    
    xPointSep = xIncVal / nSubPoints[i]
    for ( j in seq( nSubPoints[i] ) )
    {
      rIdx = which( plVals == subPoints[[i]][j] )[1]
      
      pointData = list(
        "xCoord" = xPoints[ xIdx ] + ( ( j - 1 ) * xPointSep ),
        "yCoord" = yPoints[ yIdx ],
        "cexVal" = valCounts[ names( valCounts ) == subPoints[[i]][j] ],
        "colVal" = plotCols[ catVals == isolateMetaData[ rIdx, catColIdxs ] ]
        )
      plotCoords[[ length( plotCoords ) + 1 ]] = pointData
    }
  }
  return( plotCoords ) 
}

PlotPoints = function( plotCoords )
{
  for ( i in seq( length( plotCoords ) ) )
  {
    points(
      x   = plotCoords[[i]]$xCoord,
      y   = plotCoords[[i]]$yCoord,
      cex = plotCoords[[i]]$cexVal,
      col = plotCoords[[i]]$colVal,
      pch = 16
      )
  }
}

# GetDateCount = function( dateOne, dateTwo )
# {
#   nYears  = dateTwo[1] - dateOne[1]
#   nMonths = dateTwo[2] - dateOne[2]
#   return( nYears * 12 + nMonths + 1 )
# }
# 
# GetDateComp = function( inDate )
# {
#   dateComps = strsplit( inDate, '-' )[[1]]
#   return( as.numeric( dateComps ) )
# }
# GetDateVals = function( annots )
# {
#   # Get the date values input in the input data-frame
#   dateVals  = unique( annots$dateVals )
#   dateVals  = dateVals[ order( dateVals ) ]
#   
#   # Get the number of months seperating the first and last date
#   firstDate = GetDateComp(  dateVals[1] )
#   lastDate  = GetDateComp( dateVals[ length( dateVals ) ] )
#   nDates    = SubtractDates( firstDate, lastDate )
#   
#   # If the number of dates is equal to the number of unique dates in the
#   # date frame, there are no dates missing, return the date values
#   if ( nDates == length( dateVals ) ) return( dateVals )
#   
#   # Make a vector to s
#   plDates    = vector( "character", nDates )
#   plDates[1] = dateVals[1]
#   yearVal    = firstDate[1]
#   monthVal   = firstDate[2]
#   for ( i in seq( 2, nDates ) )
#   {
#     if ( monthVal == 12 )
#     {
#       monthVal = 1
#       yearVal  = yearVal + 1
#     } else {
#       monthVal = monthVal + 1
#     }
#     
#     if ( monthVal < 10 )
#     {
#       plDates[i] = paste0( yearVal, "-0", monthVal )
#     } else {
#       plDates[i] = paste0( yearVal, '-', monthVal )
#     }
#     
#     
#   }
#   return( plDates )
# }
