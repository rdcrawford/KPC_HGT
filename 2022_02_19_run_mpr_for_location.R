# ------------------------------------------------------------------------------
# Run MPR for location
# 2021/12/03
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# Using the core gene tree and the corresponding location from which each
# isolate was sourced, find the edge corresponding to the point which each
# isolate was imported into each state and facility
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( ape )

# ---- Constant declarations ---------------------------------------------------

ACQUISTION = 1 # Integer indicating an acquisition occurred on an edge
DELETION   = 2 # Integer indicating an deletion occurred on an edge
NO_CHANGE  = 0 # Integer indication that there was no change between nodes

# ---- Import the data ---------------------------------------------------------

# Read in the meta-data
isolateMetaData = read.table(
  "../data/isolateMetaData.tsv",
  sep = '\t',
  stringsAsFactors = FALSE,
  header = TRUE
  )
load( "../data/eipData.rData" )

# Read in the core gene tree
cgTree = read.tree( "../data/cgTree.tre" )

# Load the annotations and colors: "annots" and "annotCols"
load( "../data/2022_02_07_make_annotations_and_colors.rData" )

# ---- Function definitions ----------------------------------------------------

source( "R/FindTransitonEdgeFns.R" )

source( "R/TreeParsingFunctions.R" )

source( "R/GetTreeOrder.R" )

GetStrainId = function( annots )
{
  strainId = sapply( 1:nrow( annots ), function(i)
  {
    if ( annots$ST[i] != "Other" ) return( annots$ST[i] )
    return( annots$species[i] )
  })
  return( strainId )
}

GetStrainMat = function()
{
  strainId   = GetStrainId( annots )
  strainVals = unique( strainId )
  strainVals = strainVals[ order( strainVals ) ]
  strainMat  = sapply( strainVals, function(x) as.integer( strainId == x ) )
  row.names( strainMat ) = row.names( annots )
  colnames( strainMat )  = strainVals
  return( strainMat )
}

GetStateMat = function()
{
  for ( i in which( is.na( isolateMetaData$state ) ) )
    isolateMetaData$state[i] = ''
  states = c( "Connecticut", "Minnesota", "Tennessee" )
  stateMat = sapply( states,
    function(x) as.integer( isolateMetaData$state == x )
    )
  row.names( stateMat ) = isolateMetaData$genomeNames
  colnames( stateMat ) = states
  return( stateMat )
}

GetFacMat = function()
{
  facCounts = table( eipData$Fac_ID )
  facCounts = facCounts[ facCounts > 1 ]
  facMat    = sapply( names( facCounts ), function(x)
  {
    facSampleIds = eipData$Sample_ID[ eipData$Fac_ID == x ]
    return( as.integer( isolateMetaData$genomeNames %in% facSampleIds ) )
  })
  row.names( facMat ) = isolateMetaData$genomeNames
  colnames( facMat ) = names( facCounts )
  return( facMat )
}

FindOutGroup = function( tree )
{
  tree = midpoint.root( tree )

  treeOrder = GetTreeOrder( tree, tree$tip.label )

  outSideTips =
    tree$tip.label[ c( which.min( treeOrder ), which.max( treeOrder ) ) ]

  edgeIdx = sapply( outSideTips, GetTipEdge, tree = tree )

  return( outSideTips[ which.max( tree$edge.length[ edgeIdx ] ) ] )
}

GetFeatureEdgeIdxs = function( tree, edgeFeatList )
{
  nEdges   = nrow( tree$edge )
  edgeIdxs = vector( "integer", nEdges )
  for ( j in 1:length( edgeFeatList ) )
    for ( i in edgeFeatList[[j]] ) edgeIdxs[i] = j

  return( edgeIdxs )
}

GetTipTrEdges = function( cgTree, trMat )
{
  trEdge = vector( "integer", length( isolateMetaData$genomeNames ) )
  for ( j in 1:ncol( trMat ) )
  {
    acqEdges = which( trMat[ , j ] == ACQUISTION )
    delEdges = which( trMat[ , j ] == DELETION )
    for ( x in acqEdges )
    {
      brIsos  = GetOnBranchIsolates( cgTree, x, delEdges )
      isBrIso = isolateMetaData$genomeNames %in% brIsos
      for ( i in which( isBrIso ) ) trEdge[i] = x
    }
  }
  return( trEdge )
}

source( "R/PlotTree.R" )

# ---- Find the transition edges from the plasmid trees ------------------------

# Make a matrix with a binary presence absence matrix of the genome x state
# and facility
arFeatureMats = list( GetStrainMat(), GetStateMat(), GetFacMat() )
print( sapply( arFeatureMats, dim ) )

# Run MPR to determine the edge on which an isolate was imported into each
# region. For each of the features, generate and edge by feature matrix.
trMatList = lapply( arFeatureMats,
  function(x) GetTransitionEdges( cgTree, x, "mpr", "2059860.3" )
  )

# For each matrix with the predicted acquisitions and deletions, create a list
# with all of the edges each feature was present. Return the list of lists.
edgeFeatList = lapply( trMatList, function(x)
{
  edgeList = sapply( 1:ncol(x),
    function(j) GetFeatureTransitionEdges( cgTree, x[ , j ] )
    )
  return( edgeList )
})

# Collapse the individual features to generate a matrix of the feature assigned
# to each edge.
featEdgeMat = sapply( edgeFeatList,
  function(x) GetFeatureEdgeIdxs( cgTree, x )
  )

# Make a list of the colors for each feature
strainCols = c( annotCols$ST[ 1:4 ], annotCols$species )
colList = list(
  strainCols[ order( names( strainCols ) ) ],
  annotCols$sources,
  khroma::color( "smooth rainbow")( ncol( arFeatureMats[[3]] ) )
  )

# Make a matrix with the colors corresponding to the features on each edge
edgeColMat = sapply( 1:ncol( featEdgeMat ), function(j)
{
  colVals = sapply( featEdgeMat[, j], function(x)
  {
    if ( x == 0 ) return( "black" )
    return( colList[[j]][x] )
  })
  return( colVals )
})

# Collapse the edge x feature matrices to individual columns with the edge at
# which each acquired the corresponding feature.
tipTrEdgeMat = sapply( trMatList, function(x) GetTipTrEdges( cgTree, x ) )

# Set the names of the objects to export
colnames( tipTrEdgeMat )  = c( "strain", "state", "facility" )
row.names( tipTrEdgeMat ) = isolateMetaData$genomeNames
colnames( edgeColMat )  = c( "strain", "state", "facility" )
colnames( featEdgeMat ) = c( "strain", "state", "facility" )

# ---- Save the data -----------------------------------------------------------

print( pryr::mem_used() )
rData = "../data/2022_02_19_run_mpr_for_location.rData"
save(
  file = rData,
  list = c(
    "arFeatureMats", "edgeColMat", "featEdgeMat", "edgeFeatList",
    "tipTrEdgeMat"
    )
  )

# ---- Make the figures --------------------------------------------------------

figDir = "../figures/"
figPaths = paste0(
  figDir, "2022_02_28_cg_tree_mpr_results_", colnames( edgeColMat ), ".png"
  )
annots    = annots[ , 1:3 ]
annotCols = annotCols[ 1:3 ]


for ( j in 1:ncol( edgeColMat ) )
{
  tiff(
    filename = figPaths[j],
    units    = "in",
    width    = 3.75,
    height   = 5.25,
    res      = 300
    )
  PlotTree( cgTree, annots, annotCols, edgeColMat[ , j ], 1.5 )
  dev.off()
}

# ------------------------------------------------------------------------------
