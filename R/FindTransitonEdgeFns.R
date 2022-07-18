# ------------------------------------------------------------------------------
# Get Transition Edges
# 2021/12/02
# Ryan D. Crawford
# ------------------------------------------------------------------------------

# ---- Constant declarations ---------------------------------------------------

ACQUISTION = 1 # Integer indicating an acquisition occurred on an edge
DELETION   = 2 # Integer indicating an deletion occurred on an edge
NO_CHANGE  = 0 # Integer indication that there was no change between nodes
LHS_NODE   = 1 # Column index of the left node in the edge matrix
RHS_NODE   = 2 # Column index of the right node in the edge matrix
LOWER      = 1 # Column index of lower predicted ancestral state from "MPR"
UPPER      = 2 # Column index of upper predicted ancestral state from "MPR"

# ---- Function definitions ----------------------------------------------------

# Subset the data to the relevant features which are present at least twice
# and are not present in every genome in the tree. Return the relevant
# entries sorted by the tip labels in the tree.
SortFeatureMat = function( tree, featureMat )
{
  # Subset to the isolate in the tree
  featureMat = featureMat[ row.names( featureMat ) %in% tree$tip.label, ]
  # counts     = sapply( 1:ncol( featureMat ),
  #   function(j) sum( featureMat[ , j ] > 0 )
  #   )
  # isThere    = sapply( counts, function(x) 1 < x && x < nrow( featureMat ) )

  # Sort the input matrix so that the order is the same as the tip labels on
  # the tree and subset to the relevant plasmids
  treeOrder = sapply( tree$tip.label,
    function(x) which( row.names( featureMat ) == x )
    )
  return( featureMat[ treeOrder ,  ] ) #isThere
}

# For an input matrix of discrete character states, identify the edges on the
# input tree on which the each of the features are gained or lost. Return an
# edge by feature matrix for each relevant column in the input matrix.
GetTransitionEdges = function( tree, featureMat, method, outGroup, zeroRoot )
{
  # Parse the input arguments
  if ( missing( method ) || method == "ace" )
  {
    # If using the ace method the tree is required to be rooted. If the input
    # tree is not rooted, root at the midpoint
    if ( !is.rooted( tree ) ) tree = midpoint.root( tree )

    # Set the function pointer to ace: Ancestral Character Estimation
    arFn = GetAce

    # Initialize the outgroup to null to have a placeholder variable
    if ( missing( outGroup ) ) outGroup = NULL

  } else if ( method == "mpr" ) {

    #
    arFn = GetMpr

    #
    if ( missing( outGroup ) )
      stop( "When using mpr an outgroup must be specified" )

    #
    if ( is.rooted( tree ) ) tree = unroot( tree )
  } else {
    stop( "Method must be \"ace\" or \"mpr\"" )
  }

  # Set up multi-threadding via future
  # plan = future::tweak( future::multicore, workers = threadVal )
  future::plan( future::multicore )
  options( future.globals.maxSize = Inf )

  if ( missing( zeroRoot ) ) zeroRoot = FALSE

  # If requested, zero the out group for all of the features
  if ( zeroRoot ) featureMat = ZeroOutGroup( featureMat, outGroup )

  # Sort the input matrix to match the order of the tip labels and subset to
  # the relevant features for this tree
  featureMat = SortFeatureMat( tree, featureMat )

  # Make node labels for the tree, so you can be sure the nodes are sorted
  # correctly later
  tree = makeNodeLabel( tree )

  # Find where acquisitions or deletions occur on each edge. Generate the
  # results as an edge by feature matrix
  trMat = future.apply::future_sapply( 1:ncol( featureMat ), function(j)
  {
    # Call the ancestral reconstruction for this feature and get the ancestral
    # character states for each node
    nodeChrSts = arFn( featureMat[ , j ], tree, outGroup )
    if ( TRUE %in% is.na( nodeChrSts ) )
    {
      save( file = "../data/GetTransitionEdges.rData", list = ls() )
      stop( "Error at: ", j, '\n' )
    }
    # Find the transition edges using the character states of the nodes
    return( FindEdgeTrStatus( tree, nodeChrSts, featureMat[ , j ] ) )
  }, future.seed = TRUE )

  # Set the dimnames to the names of the genes and the index of the edges
  colnames( trMat )  = colnames( featureMat )
  row.names( trMat ) = seq( nrow( trMat ) )

  return( trMat )
}

# For each edge, use the ancestral states predicted at each node to determine
# if there was a transition at each node.
FindEdgeTrStatus = function( tree, nodeChrSts, features )
{
  # Make a vector with the indexes of the node values
  nodeIds = seq( tree$Nnode ) + length( tree$tip.label )

  # For each edge in the matrix find the transition status
  edgeTrStatus = sapply( 1:nrow( tree$edge ), function(i)
  {
    # Look up the status of the input feature on the left node in the tree
    lhs = nodeChrSts[ nodeIds == tree$edge[ i, LHS_NODE ] ]

    # Look up the status of the right node or tip
    if ( tree$edge[ i, RHS_NODE ] <= length( tree$tip.label ) )
    {
      rhs = features[ tree$edge[ i, RHS_NODE ] ]
    } else {
      rhs = nodeChrSts[ nodeIds == tree$edge[ i, RHS_NODE ] ]
    }
    return( GetEdgeTrDir( lhs, rhs ) )
  })

  return( edgeTrStatus )
}

# Assigned to "arFn" to perform maximum parsimony ancestral character estimation
GetMpr = function( features, tree, outGroup )
{
  mpr    = MPR( features, tree, outGroup )
  ndIdxs = as.integer( gsub( "Node", '', row.names( mpr ) ) )
  mpr    = mpr[ order( ndIdxs ),  ]
  return( mpr[ , LOWER ] )
}

# Assigned to "arFn" to perform Ancestral character estimation
GetAce = function( features, tree, outGroup )
{
  arResults = ace( features, tree, type = "discrete", model = "ARD" )

  isBad = sapply( c( NA, NaN, Inf ), function(x) x %in% arResults$lik.anc )
  nFeatures = sum( features > 0 )
  fVals = unique( features )
  # check the output:
  if (
    nrow( arResults$lik.anc ) != tree$Nnode ||
    TRUE %in% isBad ||
    length( fVals ) != ncol( arResults$lik.anc )
    )
  {

    cat(
      "Error in ace:\n",
      "  -- Number of tips with feature: ", nFeatures, '\n',
      "  -- Number of nodes: ", tree$Nnode, '\n',
      "  -- Numer of constructed nodes: ", nrow( arResults$lik.anc ), '\n',
      "  -- Feature levels: ", paste( fVals, collapse = ", " ) , '\n',
      "  -- Number of constructed features: ", ncol( arResults$lik.anc ), '\n',
      "  -- Invalid character: ", TRUE %in% isBad, '\n',
      sep = ''
      )
    return( NA )
  }
  arResults

  fVals = unique( features )
  fVals = fVals[ order( fVals ) ]

  nodeChrSts = sapply( 1:nrow( arResults$lik.anc ),
    function(i) fVals[ which.max( arResults$lik.anc[ i, ] ) ]
    )
  return( nodeChrSts )
}

# From the input nodes connected by an edge on the tree, determine the
# the transition status
GetEdgeTrDir = function( lhs, rhs )
{
  # Return whether there was no change on this edge, or a transition occurred
  if ( lhs == rhs ) return( NO_CHANGE )
  if ( lhs < rhs ) return( ACQUISTION )
  return( DELETION )
}

# Make all features in the matrix 0 for the out group so that there
# will be at least one acquisition
ZeroOutGroup = function( featureMat, outGroup )
{
  ogIdx = which( row.names( featureMat ) == outGroup )
  for ( j in 1:ncol( featureMat ) ) featureMat[ ogIdx, j ] = 0
  return( featureMat )
}

# ------------------------------------------------------------------------------