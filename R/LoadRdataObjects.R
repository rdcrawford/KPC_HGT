# ------------------------------------------------------------------------------
# LoadRdataObjects
# 2021/08/11
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# Assigen variables in an R-data object to the global environment
# ------------------------------------------------------------------------------

LoadRdataObjects = function( rDataPath, varNames )
{
  # Load the object from the path
  load( rDataPath )
  
  # Get a list of the variables in the function envrionment
  vars    = ls()
  vars    = vars[ !vars %in% c( "rDataPath", "varNames" ) ]
  isDescr = vars == "description"
  
  # If this contains a description, print it
  if ( TRUE %in% isDescr ) cat( description, '\n' )
  
  # If only selected variables are to be loaded, load each of them
  if ( missing( varNames ) ) varNames = vars[ !isDescr ]
  
  # Assign each variabe to the global environment
  for ( x in varNames ) assign( x, get( x ), envir = .GlobalEnv )
}

# -----------------------------------------------------------------------------

