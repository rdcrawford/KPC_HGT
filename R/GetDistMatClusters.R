# ------------------------------------------------------------------------------
# Shared plasmid analysis
# 2021/07/07
# Ryan D. Crawford
# ------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------

# Look up the indexes in the list which share overlapping elements. 
FindOverlapIdxs = function( isShared, idxs )
{
  isClade = sapply( seq( ncol( isShared ) ),
    function(j) TRUE %in% isShared[ idxs, j ]
    )
  sIdxs = which( isClade )
  if ( !FALSE %in% ( sIdxs %in% idxs ) ) return( idxs )
  return( FindOverlapIdxs( isShared, sIdxs ) )
}

# Create a list of isolates which all share the same nearest neighbors
GetDistMatClusters = function( inMat, threshVal )
{
  idx      = 1
  nIsos    = nrow( inMat )
  clIds    = vector( "integer", nIsos )
  isShared = sapply( seq( nIsos ), function(j) inMat[ , j ] <= threshVal )

  while ( 0 %in% clIds )
  {
    clIdxs = FindOverlapIdxs( isShared, which( clIds == 0 )[1] )

    for ( i in clIdxs ) clIds[i] = idx
    idx = idx + 1
  }

  names( clIds ) = row.names( inMat )
  return( clIds )
}

# ------------------------------------------------------------------------------