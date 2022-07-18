
CalcUnion = function( lhs, rhs ) return( sum( lhs %in% rhs ) )

CalcUnionMat = function( inList, geneIds )
{
  n = length( inList )
  unMat = matrix( 1, n, n )
  row.names( unMat ) = geneIds
  colnames( unMat )  = geneIds

  for ( i in seq( n  - 1 ) )
  {
    unMat[ i, i ] = length( inList[[i]] )
    for ( j in seq( i + 1, n ) )
    {
      un = CalcUnion( inList[[i]], inList[[j]] )
      unMat[ i, j ] = un
      unMat[ j, i ] = un
    }
  }
  return( unMat )
}

# Sørensen-dice index
CalcSorensenDiceIdx = function( unMat, nGenes )
{
  n      = ncol( unMat )
  sdMat  = matrix( 0, n, n )
  row.names( sdMat ) = colnames( unMat )
  colnames( sdMat )  = colnames( unMat )

  for ( i in seq( n  - 1 ) )
  {
    for ( j in seq( i + 1, n ) )
    {
      sdIdx = ( 2 * unMat[ i, j ] ) / ( nGenes[i] + nGenes[j] )
      if ( sdIdx > 1 ) stop( i,  ', ', j )
      sdMat[ i, j ] = sdIdx
      sdMat[ j, i ] = sdIdx
    }
  }
  return( sdMat )
}

# Jaccard index
CalcJacIdx = function( unMat, nGenes )
{
  n      = ncol( unMat )
  jacMat = matrix( 1, n, n )
  row.names( jacMat ) = colnames( unMat )
  colnames( jacMat )  = colnames( unMat )

  for ( i in seq( n  - 1 ) )
  {
    for ( j in seq( i + 1, n ) )
    {
      jIdx = unMat[ i, j ] / ( nGenes[i] + nGenes[j] - unMat[ i, j ] )
      jacMat[ i, j ] = jIdx
      jacMat[ j, i ] = jIdx
    }
  }
  return( jacMat )
}

MakeJaccardMat = function( inList, listNames )
{
  unMat = CalcUnionMat( inList, listNames )
  nElms = sapply( inList, length )
  return( CalcJacIdx( unMat, nElms ) )
}