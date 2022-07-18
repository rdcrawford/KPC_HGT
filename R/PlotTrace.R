# ------------------------------------------------------------------------------
# Target Plot
# 2021/07/24
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# Make a scatter plot where points can be nested within each other recursively
# to illustrate the prevalence of multiple variables.
# ------------------------------------------------------------------------------

# plData     = plData
# xColIdx    = xColIdx
# xAxisVals  = xAxisVals
# yColIdx    = yColIdx
# yAxisLabs  = yAxisLabs
# catColIdxs = c( regColIdxs, kpcColIdx )
# plotCols   = plotCols[ 1:2 ]
# logScale   = TRUE
# logBase    = 10
# baseCexVal = 1
# cexAdjVal  = 0.75

# Make a scatter plot where the size of the points corresponds to the
# number of observations with the intersecting XY values
TargetPlot = function(
  plData,     # Data-frame with the relevant values to plot
  xColIdx,    # Column index of the dates to plot on the x-axis
  xAxisVals,  # Vector with the variables to plot on the x-axis
  xAxisLabs,  # Labels to plot on the x-axis
  yColIdx,    # Column index of the dates to plot on the y-axis
  yAxisVals,  # Vector with the variables to plot on the y-axis
  yAxisLabs,  # Labels to plot on the x-axis
  catColIdxs, # Column indexes of the categories to plot, in order
  plotCols,   # List of colors to be used for the categories to plot
  logScale,   # Logical to take the log of the count values
  logBase,    # Base value for the logarithm. Defaults to "exp(1)"
  baseCexVal, # The size of the a point with a value equal to one
  cexAdjVal   # The fraction to remove from the inner circles
  )
{
  # Set the default variable values
  # if ( missing( xOffSet ) )    xOffSet    = 0
  # if ( missing( yOffSet ) )    yOffSet    = 0
  if ( missing( catColIdxs ) ) catColIdxs = NULL
  if ( missing( logScale ) )   logScale   = FALSE
  if ( missing( logBase ) )    logBase    = exp(1)
  if ( missing( plotCols ) )   plotCols   = list( "black" )
  if ( missing( baseCexVal ) ) baseCexVal = 1
  if ( missing( cexAdjVal ) )  cexAdjVal  = 0.75

  # Make the function that determines the size of the plot points
  GetCexVal = MakeGetCexValFn( logScale, logBase, baseCexVal, cexAdjVal )

  # Get the points for the x and y axis
  if ( missing( xAxisVals ) )
    xAxisVals = GetPointCatVals( plData, xColIdx )
  if ( missing( yAxisVals ) )
    yAxisVals = GetPointCatVals( plData, yColIdx )

  # Transform the xAxis values into a list
  if ( !is.list( xAxisVals ) ) xAxisVals = list( xAxisVals )

  # Set the values for the axis
  if ( missing( xAxisLabs ) ) xAxisLabs = xAxisVals[[1]]
  if ( missing( yAxisLabs ) ) yAxisLabs = yAxisVals[[1]]

  # Make a list of the values to plot
  catValList = GetPointCatVals( plData, catColIdxs )
  valList    = c( xAxisVals, yAxisVals, catValList )

  # Sort the colors to match the values in the list
  plotCols = SortColors( plotCols, catValList, plData, catColIdxs )

  # Get the points to plot along the x and y axis\
  nXPts     = length( xAxisLabs )
  xPoints   = GetAxisPoints( nXPts + 6 )
  legPoints = xPoints[ ( nXPts + 1 ):length( xPoints ) ]
  xPoints   = xPoints[ 1:nXPts ]
  yPoints   = GetAxisPoints( length( yAxisLabs ) )

  # Make a list with the coordinates for each point on the plot
  plotCoords = GetPointCoords(
    plData, valList, xColIdx, yColIdx, xPoints, yPoints, catColIdxs, plotCols,
    GetCexVal
    )

  # Get the size of the dots to plot for the legend
  maxVal  = max( sapply( plotCoords, function(x) x$cexVal ) )
  cexVals = GetLegendCexVals( maxVal, logScale, logBase, baseCexVal, GetCexVal )

  # Create the plot
  plot.new()
  PlotAxis( xPoints, xAxisLabs, 'x' )
  PlotAxis( yPoints, yAxisLabs, 'y' )
  PlotLegend( cexVals, legPoints, yPoints, catValList, plotCols )
  PlotXyPoints( plotCoords )
  PlotGrid( xPoints, yPoints )
  PlotXyPoints( plotCoords )
  PlotLegend( cexVals, legPoints, yPoints, catValList, plotCols )
}

# Plot a grid of the xy values to be the background of the plot
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

# Get a vector of the position on the axis for 'm' points to plot
GetAxisPoints = function( n )
{
  return( seq( n ) / n )
}

# Plot the axis labels on the 'x' or 'y' axis specified by the 'axisVal'
# The labels specified by 'labelVals' are plotted at the points along the
# axis specified by 'pointVals'
PlotAxis = function( pointVals, labelVals, axisVal, offSetVal )
{
  # Find the side of the plot to place the axis
  if ( axisVal == 'x' )
  {
    sideVal = 1
  } else if ( axisVal == 'y' ) {
    sideVal = 2
  } else {
    stop( "axis must be x or y\n" )
  }

  if ( missing( offSetVal ) ) offSetVal = 0

  # Plot the axis and labels
  axis(
    side   = sideVal,
    at     = pointVals,
    labels = labelVals,
    las    = 2,
    pos    = ,
    outer  = FALSE
    )
}

# Return a function to generate the the size of the point
MakeGetCexValFn = function( logScale, logBase, baseCexVal, cexAdjVal )
{
  if ( logScale )
  {
    GetCexVal = function( nVal, ptIdx )
    {
      cexVal = log( nVal, base = logBase ) + baseCexVal
      if ( missing( ptIdx ) ) return( cexVal )
      # Get the size of the point
      adjVal = 1 * ( cexAdjVal**( ptIdx - 1 ) )
      return( cexVal * adjVal )
    }
  } else {
    GetCexVal = function( nVal, ptIdx )
    {
      cexVal = nVal + baseCexVal
      if ( missing( ptIdx ) ) return( cexVal )
      # Get the size of the point
      adjVal = 1 * ( cexAdjVal**( ptIdx - 1 ) )
      return( cexVal * adjVal )
    }
  }

  return( GetCexVal )
}

# Plot the values corresponding to the size of the points and the
# categories which were plotted
PlotLegend = function( cexVals, legPoints, yPoints, catValList, plotCols )
{
  # cexVals    = rev( cexVals )
  colVals    = unlist( plotCols, use.names = FALSE )
  colNames   = unlist( catValList, use.names = FALSE )
  legVals    = c( colNames, names( cexVals ) )
  legCols    = c( colVals, rep( "black", length( cexVals ) ) )
  colCexVal  = 2 #cexVals[ ceiling( length( cexVals ) / 2 ) ]
  legCexVals = c( rep( colCexVal, length( colVals ) ), cexVals )
  # yTextPos   = 1 - ( yPoints[2] - yPoints[1]  )
  points(
    x   = rep( legPoints[1], length( legVals ) ),
    y   = yPoints[ 1:length( legVals ) ],
    cex = legCexVals,
    col = legCols,
    pch = 16
    )
  txtPos = legPoints[1] + ( legPoints[2] - legPoints[1] ) / 2
  text(
    legVals,
    x   = rep( txtPos, length( legVals ) ),
    y   = yPoints[ 1:length( legVals ) ],
    col = "black",
    cex = 0.5,
    pos = 4
    )
}

# Generate a range of values spanning the range of the size of the points
GetLegendCexVals = function( maxVal, logScale, logBase, baseCexVal, GetCexVal )
{
  # If this is the log of the maximium
  if ( logScale )
  {
    maxPtVal = round( logBase**( maxVal - baseCexVal ), 0 )
  } else {
    maxPtVal = maxVal - baseCexVal
  }

  if ( maxPtVal <= 10 )
  {
    nVals = maxPtVal / 2
  } else {
    nVals = 5
  }

  legCexVals = round( seq( 1, maxPtVal, length.out = nVals ) )
  cexVals = sapply( legCexVals, GetCexVal )
  names( cexVals ) = legCexVals
  return( cexVals )
}

# Find the unique values for each input category excluding NA values
GetPointCatVals = function( plData, catColIdxs )
{
  catValList = sapply( catColIdxs, function(j)
  {
    vals = unique( plData[ !is.na( plData[ , j ] ), j ] )
    return( list( vals[ order( vals ) ] ) )
  })
  return( catValList )
}

# For each category, find the the poiintersection of all of the points
LookUpPoints = function(
  idx, colIdxs, plData, valList, valIdxs = NULL, rIdxs = NULL
  )
{
  # If all of the columns to plot have been considered, return a list
  # containing a list of values, and a list of row indexes which have the
  # corresponding values
  if ( idx > length( colIdxs ) ) return( list( list( valIdxs, rIdxs ) ) )

  # Get the current column and the values contained in that column
  cIdx   = colIdxs[ idx ]
  ptVals = valList[[ idx ]]

  # Find which values are contained by these observations. Generate a list of
  #
  if ( is.null( rIdxs ) )
  {
    curIdxs = which( ptVals %in% plData[ , cIdx ] )
    idxList =
      lapply( curIdxs, function(x) which( plData[ , cIdx ] == ptVals[x] ) )

  } else {

    # Find which of the point values are in this column
    curIdxs = which( ptVals %in% plData[ rIdxs, cIdx ] )
    idxList =
      lapply( curIdxs, function(x) rIdxs[ plData[ rIdxs, cIdx ] == ptVals[x] ] )
  }

  nxtIdx = idx + 1
  ptIdxList = list( )
  for ( i in 1:length( idxList ) )
  {
    ptList = LookUpPoints(
      nxtIdx, colIdxs, plData, valList, c( valIdxs, curIdxs[i] ), idxList[[i]]
      )
    ptIdxList = append( ptIdxList, ptList )
  }

  return( ptIdxList )
}

# Create a list of points for the scatter plot
GetPointCoords = function(
  plData, valList, xColIdx, yColIdx, xPoints, yPoints, catColIdxs, plotCols,
  GetCexVal
  )
{
  # Get the distance between points on the xAxis
  xIncVal = xPoints[2] - xPoints[1]

  # Find the points with intersceting features
  colIdxs   = c( xColIdx, yColIdx, catColIdxs )
  plIdxList = LookUpPoints( 1, colIdxs, plData, valList )

  # Get the x and y values from
  xIdxVals = sapply( plIdxList, function(x) x[[1]][1] )
  yIdxVals = sapply( plIdxList, function(y) y[[1]][2] )

  # Initialize the list of plot coordinates
  plotCoords = vector( "list", 0 )

  for ( x in xIdxVals )
  {
    for ( y in yIdxVals )
    {
      # Find the spints with shared x and y value
      isXyVal = xIdxVals == x & yIdxVals == y

      if ( TRUE %in% isXyVal )
      {
        # At the position on the x-y grid, generate the points to plot
        xpPoints = CreatePointsAtXyPos(
          plIdxList[ isXyVal ], xPoints, yPoints, plotCols, xIncVal, GetCexVal
          )

        # Append the points to the output
        plotCoords = append( plotCoords, xpPoints )
      }
    }
  }

  return( plotCoords )
}

# For the current category Create a list of points for the plot
CreatePointsAtXyPos = function(
  xyValList, xPoints, yPoints, plotCols, xIncVal, GetCexVal
  )
{
  #
  nPts   = length( xyValList )
  adjVal = xIncVal / ( nPts + 1 )
  nVal   = sapply( xyValList, function(x) length( x[[2]] ) )
  xVal   = xPoints[ xyValList[[1]][[1]][1] ]
  yVal   = yPoints[ xyValList[[1]][[1]][2] ]

  # Get the values for this point
  xyPoints = sapply( seq( nPts ), function(i)
  {
    ptValIdxs = xyValList[[i]][[1]]
    xVal = xVal + adjVal * ( i - 1 )
    point =
      CreatePointVal( ptValIdxs, nVal[i], xVal, yVal, plotCols, GetCexVal )
    return( point )
  })

  # Return the points to plot sorted by size
  return( SortPoints( xyPoints ) )
}

# Generate a "point" to be plotted at the specified x and y position with
# the specified color. Additionally, if there are any "sub points" generate
# them
CreatePointVal = function(
  ptValIdxs, nVal, xVal, yVal, plotCols, GetCexVal, sPtIdx = 1
  )
{
  if ( sPtIdx > length( plotCols ) ) return( list() )
  point           = list()
  point$xCoord    = xVal
  point$yCoord    = yVal
  point$colVal    = plotCols[[ sPtIdx ]][ ptValIdxs[ sPtIdx + 2 ] ]
  point$cexVal    = GetCexVal( nVal, sPtIdx )
  point$subPoints = CreatePointVal(
    ptValIdxs, nVal, xVal, yVal, plotCols, GetCexVal, sPtIdx + 1
    )
  return( list( point ) )
}

# Sort a list of points by the "cexVal" attribute in decresing order
SortPoints = function( xyPoints )
{
  n = sapply( xyPoints, function(x) x$cexVal )
  return( xyPoints[ order( n, decreasing = TRUE ) ] )
}

# Plot a point with the specicified attributes
PlotPoint = function( point )
{
  # Plot the point using the specified coordinates and attributes from the
  # input list
  points(
    x   = point$xCoord,
    y   = point$yCoord,
    cex = point$cexVal,
    col = point$colVal,
    pch = 16
    )

  # If this point has any associated sub-points plot them
  if ( length( point$subPoints ) ) PlotPoint( point$subPoints[[1]] )
}

# Input the list of points and make the scatter plot
PlotXyPoints = function( plotCoords )
{
  for ( i in seq( length( plotCoords ) ) ) PlotPoint( plotCoords[[i]] )
}

# Sort the input colors to match the categories in the list
SortColors = function( plotCols, catValList, plData, catColIdxs )
{
  isIn = names( plotCols ) %in% colnames( plData )[ catColIdxs ]
  if ( !FALSE %in% isIn && !is.null( names( plotCols ) ) )
  {
    idxs = sapply( names( plotCols ),
      function(x) which( colnames( plData )[ catColIdxs ] == x )
      )
    plotCols = plotCols[ order( idxs ) ]
  }

  for ( i in 1:length( plotCols ) )
  {
     plotCols[[i]] =
       plotCols[[i]][ names( plotCols[[i]] ) %in% catValList[[i]] ]
    idxs = sapply( names( plotCols[[i]] ),
      function(x) which( catValList[[i]] == x )
      )
    plotCols[[i]] = plotCols[[i]][ order( idxs ) ]
  }
  return( plotCols )
}

# ------------------------------------------------------------------------------
