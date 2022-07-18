MakeHmCols = function( inMat )
{
  minVal = min( inMat )
  if ( minVal < 0 )
  {
    maxVal = max( inMat )
    nNeg   = which( seq( minVal, maxVal, length.out = 100 ) >= 0 )[1]
    nPos   = 100 - nNeg + 1
    mVal   = max( c( nPos, nNeg ) )
    negRng = ( mVal- nNeg + 1 ):mVal

    # Make a vector of heatmap colors
    hmCols = c(
      colorRampPalette( c( "#395983", "#e7e7e7" ) )( mVal )[ negRng ],
      colorRampPalette( c( "#e7e7e7", "#952827" ) )( mVal )[ 2:nPos ]
      )
    return( hmCols )
  }
  return( colorRampPalette( c( "#395983", "#e7e7e7", "#952827" ) )( 100 ) )
}
