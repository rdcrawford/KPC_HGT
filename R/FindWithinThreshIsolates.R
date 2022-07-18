# ------------------------------------------------------------------------------
# Make cognac output files
# 2021/05/09
# Ryan D. Crawford
# ------------------------------------------------------------------------------

FindWithinThreshIsolates = function( distMat, varThresh )
{
  # Create a logical vector to indicate which genomes are classigied
  isIn = vector( "logical", nrow( distMat ) )
  
  # Thi variable keeps track of the curent genome to test
  idx = 1
  
  # Initialize a list with the representitives within each cluster
  clustList = vector( "list", 0 )
  
  # While there are still observations in th matirx to test..
  while ( !is.na( idx ) ) 
  {
    isIn[ idx ] = TRUE
    
    isQuery = which( !isIn )
    isInCl  = distMat[ idx, isQuery ] <= varThresh
    listIdx = length( clustList ) + 1
    
    if ( TRUE %in% isInCl )
    {
      for ( i in isQuery[ isInCl ] ) isIn[i] = TRUE
      clustReps = 
        c( row.names( distMat )[ idx ], row.names( distMat )[ isQuery ][ isInCl ] )
    } else {
      clustReps = vector( "character", 0 )
    }
    
    clustList[ listIdx ] = list( clustReps )
    names( clustList )[ listIdx ] = row.names( distMat )[ idx ]
    
    # if ( length( clustReps ) )
    # {
    #   isIso  = row.names( distMat ) %in% clustReps
    #   clDist =  distMat[ idx, isIso ] 
    #   if ( min( clDist ) > varThresh )
    #   {
    #     cat( "The minimium variant threshold at", idx, "is", min( clDist ), '\n' )
    #     stop( )
    #   }
    # }
    
    idx =  isQuery[ !isInCl ][1]
    
  }
  
  return( clustList )
}

# ------------------------------------------------------------------------------

