# ------------------------------------------------------------------------------
# FindGeneClusters
# 2021/09/01
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# 
# ------------------------------------------------------------------------------

FindGeneClusters = function( distMat, distThreshold )
{
  clReps    = vector( "integer", 0 )
  geneClIds = vector( "integer", nrow( distMat ) )

  # Iterate over rows in the distance matrix
  for ( i in 1:nrow( distMat ) )
  {
    bestCl       = GetGeneMatch( i, distMat, clReps, distThreshold )
    geneClIds[i] = bestCl
    if ( bestCl > length( clReps ) ) clReps[ length( clReps ) + 1 ] = i
  }
  return( geneClIds )
}

# FindGeneClusters = function( distMat, distThreshold )
# {
#   clReps    = vector( "list", 0 )
#   geneClIds = vector( "integer", nrow( distMat ) )
# 
#   # Iterate over rows in the distance matrix 
#   for ( i in 1:nrow( distMat ) )  
#   {
#     bestCl       = GetGeneMatch( i, distMat, clReps, distThreshold )
#     geneClIds[i] = bestCl
#     if ( bestCl > length( clReps ) )
#     {
#       clReps[[ length( clReps ) + 1 ]] = i
#     } else {
#       clReps[[ bestCl ]][ length( clReps[[ bestCl ]] ) + 1 ] = i
#     }
#   }
#   
#   names( geneClIds ) = row.names( distMat )
#   return( geneClIds )
# }

GetGeneMatch = function( i, distMat, clReps, distThreshold )
{
  if ( length( clReps ) == 0 ) return( 1 )
  bestIdx = 0
  bestMatch = 1
  for ( j in 1:length( clReps ) )
  {
    minVal = min( distMat[ i, clReps[j] ])
    if ( minVal <= distThreshold )
    {
      if ( bestMatch > minVal )
      {
        bestMatch = minVal
        bestIdx   = j
      }
    }
  }

  if ( bestIdx == 0 ) return( length( clReps ) + 1 )
  return( bestIdx )
}

# ------------------------------------------------------------------------------

#   
#   names( geneClIds ) = row.names( distMat )
#   return( geneClIds )
# }
# 
# GetGeneMatch = function( i, distMat, clReps, distThreshold )
# {
#   cat( "Finding best match:", i, '\n' )
#   if ( length( clReps ) == 0 ) return( 1 )
#   bestIdx = 0
#   bestMatch = 1
#   for ( j in 1:length( clReps ) )
#   {
#     if ( distMat[ i, clReps[j] ] <= distThreshold )
#     {
#       if ( bestMatch > distMat[ i, clReps[j] ] )
#       {
#         
#         cat(  
#           "  -- New match:", bestIdx, '(', bestMatch, '>',
#           distMat[ i, clReps[j] ], ')', '\n' 
#           )
#         bestMatch = distMat[ i, clReps[j] ]
#         bestIdx   = j
#       }
#     }
#   }
#   cat( "  -- returning:", bestIdx, '\n' )
#   if ( bestIdx == 0 ) return( length( clReps ) + 1 )
#   return( bestIdx )
# }