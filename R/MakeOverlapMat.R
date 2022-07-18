# ------------------------------------------------------------------------------
# MakeOverlapMat
# 2021/09/01
# Ryan D. Crawford
# ------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------

RemoveDupGenes = function( genes )
{
  nZeros = sum( genes == 0 ) 
  isDuplicated = duplicated( genes )
  genes = genes[ !isDuplicated ]
  if ( nZeros <= 1 ) return( genes ) 
  return( c( genes, rep( 0, nZeros - 1 ) ) )
}

# Calculate the overlap coefficient for a pair of input vectors
# rhs = contigData[[i]]; lhs = contigData[[j]]
GetOverlapCoef = function( lhs, rhs )
{
  # First find the length of the smaller vector. Remove any zeros from 
  # one vector (zeros represent singletons). calculate the intersection and 
  # return the similarity
  lhs       = RemoveDupGenes( lhs )
  rhs       = RemoveDupGenes( rhs )
  minLen    = min( length( rhs ), length( lhs ) )
  rhs       = rhs[ rhs != 0 ]
  intersect = sum( lhs %in% rhs )
  return( intersect / minLen )
}

# For each pair of plasmids calculate the overlap coefficient. Return
# a square symmetric matrix with these similarity values
MakeOverlapMat = function( contigData )
{
  # Initialize a square matrix to store the overpal coef values
  nPlasmids               = length( contigData )
  overlapMat              = matrix( 1, nPlasmids, nPlasmids )
  row.names( overlapMat ) = names( contigData )
  colnames( overlapMat )  = names( contigData )

  # For each pair of plasmids, calculate the overlap
  for ( i in seq( nPlasmids - 1 ) )
  {
    for ( j in seq( i + 1, nPlasmids ) )
    {
      overlapVal         = GetOverlapCoef( contigData[[i]], contigData[[j]] )
      overlapMat[ i, j ] = overlapVal
      overlapMat[ j, i ] = overlapVal
    }
  }

  return( overlapMat )
}

# ------------------------------------------------------------------------------