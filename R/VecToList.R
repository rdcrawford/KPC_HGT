# ------------------------------------------------------------------------------
# VecToList
# 2021/08/23
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# Transform a vector to a list by each category in th list. The list 
# elements contain the names of the vector for each category in the vector
# ------------------------------------------------------------------------------

VecToList = function( clusts )
{
  return( tapply( names(clusts), clusts, function(x) return(x) ) )
}

# ------------------------------------------------------------------------------
