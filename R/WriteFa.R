
# ------------------------------------------------------------------------------
# WriteFa
# 2021/09/03
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# 
# ------------------------------------------------------------------------------

WriteFa = function( algnPath, algn, seqIds )
{
  if ( missing( seqIds ) )
  {
    if ( is.null( names( algn ) ) )
    {
      seqIds = seq( length( algn ) )
    } else {
      seqIds = names( algn )
    }
  }
  
  sink( algnPath )
  for ( i in 1:length( algn ) )
    cat( '>', seqIds[i], '\n', algn[i], '\n', sep = '' )
  sink()
}

# ------------------------------------------------------------------------------
