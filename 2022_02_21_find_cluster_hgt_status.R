# ------------------------------------------------------------------------------
# Find cluster hgt status
# 2021/02/21
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# Use the MPR results for the clusters and geographic location to determine
# if presence of a given cluster represents HGT in that state/facility
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( ape )

# ---- Constant declarations ---------------------------------------------------

HGT_EVENT   = 1
IMPORTATION = 2
ABSENT      = 0
STATE       = 1
FACILITY    = 2

# ---- Import the data ---------------------------------------------------------

# Load the data on the EIP genomes
load( "../data/eipData.rData" )

# Load the annotations and colors: "annots" and "annotCols"
load( "../data/2022_02_07_make_annotations_and_colors.rData" )

# Load the mpr results on the edge in which an isolate was included in a region.
# includes: "arFeatureMats", "edgeColMat", "featEdgeMat", "edgeFeatList",
# and "tipTrEdgeMat"
load( "../data/2022_02_19_run_mpr_for_location.rData" )

# Load the data on the deduplicated eip gene clusters: "eipClMat",
# "eipAcqEdge", "plMbrList", and "plMat"
load( "../data/2022_02_16_analyze_hgt_cluster_distribution.rData" )

# Load the data on the KPC acquisition permutation test: "pValsTrEdge",
# "sigThresh", and "idxs"
load( "../data/2022_02_09_find_cluster_kpc_acquisition_enrichment_sig.rData" )

# Load the annotations and colors: "annots" and "annotCols"
load( "../data/2022_02_07_make_annotations_and_colors.rData" )

# Read in the core-gene tree
cgTree = read.tree( "../data/cgTree.tre" )

# ---- Function definitions ----------------------------------------------------

source( "R/TreeParsingFunctions.R" )

GetNumEdgesFromTip = function( acqMat, ancEdges )
{
  edgeIdxMat = sapply( 1:ncol( acqMat ), function(j)
  {

    edgeIdxs = sapply( 1:nrow( acqMat ), function(i)
    {
      cat( i, '\n' )
      if ( acqMat[ i, j ] == ABSENT ) return( ABSENT )
      return( which( ancEdges[[i]] == acqMat[ i, j ]  ) )
    })
    return( edgeIdxs )
  })
  dimnames( edgeIdxMat ) = dimnames( acqMat )
  return( edgeIdxMat )
}

FindEipHgtCls = function( eipAcqEdge, locEdges )
{
  clHgtStatus = sapply( 1:ncol( eipAcqEdge ), function(j)
  {
    clStatus = sapply( 1:nrow( eipAcqEdge ), function(i)
    {
      if ( locEdges[i] == ABSENT ) return( NA )
      if ( eipAcqEdge[ i, j ] == 0 ) return( ABSENT )
      if ( eipAcqEdge[ i, j ] <= locEdges[i] ) return( HGT_EVENT )
      return( IMPORTATION )
    })
    return( clStatus )
  })
  dimnames( clHgtStatus ) = dimnames( eipAcqEdge )

  return( clHgtStatus )
}

GetTypeCounts = function( inMat )
{
  countMat = sapply( 1:ncol( inMat ),
    function(j) sapply( seq( 2 ),
      function(x) sum( inMat[ , j ] == x, na.rm = TRUE )
      )
    )
  return( countMat )
}

GetHgtAndImpIsoCounts = function( eipAcqEdge, stHgtStatus )
{
  # Get the edges on which each cluster was acquired in the EIP genomes
  egdeVals = lapply( 1:ncol( eipAcqEdge ), function(j)
  {
    edgeVals = unique( eipAcqEdge[ , j ] )
    return( edgeVals[ edgeVals != 0 ] )
  })

  # Make a list of lists of logical vectors
  isClEdge = lapply( 1:ncol( eipAcqEdge ),
    function(j)
      lapply( egdeVals[[j]], function(x) which( eipAcqEdge[ , j ] == x ) )
    )

  edgeCounts = lapply( 1:ncol( eipAcqEdge ),
    function(j) sapply( isClEdge[[j]], length )
    )

  hgtEdgeCounts = lapply( 1:ncol( eipAcqEdge ),
    function(j) sapply( isClEdge[[j]],
      function(x) sum( stHgtStatus[ x, j ] == HGT_EVENT )
      )
    )
  impEdgeCounts = lapply( 1:ncol( eipAcqEdge ),
    function(j) sapply( isClEdge[[j]],
      function(x) sum( stHgtStatus[ x, j ] == IMPORTATION )
      )
    )
  return( list( hgtEdgeCounts, impEdgeCounts ) )
}

# ---- Get the clusters which are HGT in the EIP -------------------------------

# # Subset to the clusters significantly associated with KPC acquisition
# sigIdxs    = idxs[ pValsTrEdge <= sigThresh ]
clList     = strsplit( names( plMbrList ), '/' )
# isSig      = sapply( clList, function(x) TRUE %in% ( sigIdxs %in% x  )
# eipAcqEdge = eipAcqEdge[ , isSig ]
# eipClMat   = eipClMat[ , isSig ]

# Sort the data, and subset to the EIP genomes and clusters which are present
# in the EIP genomes.
geoLocAcqEdges = tipTrEdgeMat[ , c( 2, 3 ) ]
glRowIdxs      = sapply( row.names( eipClMat ),
  function(x) which( row.names( geoLocAcqEdges ) == x )
  )
geoLocAcqEdges = geoLocAcqEdges[ glRowIdxs, ]

# Check that the data are sorted properly
if ( !identical( row.names( geoLocAcqEdges ), row.names( eipClMat ) ) )
  stop( "row.names( geoLocAcqEdges )!= row.names( eipClMat )" )
if ( !identical( row.names( eipAcqEdge ), row.names( eipClMat ) ) )
  stop( "row.names( geoLocAcqEdges )!= row.names( eipClMat )" )

# Look up the index of the edges connected to each tip
tipEdges =
  sapply( row.names( eipAcqEdge ), function(x) GetTipEdge( x, cgTree ) )

# Get all of the edges down to the root
ancEdges = lapply( tipEdges, function(x) GetAncestorEdges( cgTree, x, NULL ) )

# Find the number of edges between each acquisition and the tips
geoLocEdgeIdxs = GetNumEdgesFromTip( geoLocAcqEdges, ancEdges )
clEdgeIdxs     = GetNumEdgesFromTip( eipAcqEdge, ancEdges )

# Determine if each cluster represents an importation or HGT
stHgtStatus  = FindEipHgtCls( clEdgeIdxs, geoLocEdgeIdxs[ , STATE ] )
facHgtStatus = FindEipHgtCls( clEdgeIdxs, geoLocEdgeIdxs[ , FACILITY ] )

# ----

stVals = c( "CT", "MN", "TN" )

stIdxList = lapply( stVals, function( st )
{
  isSt = row.names( eipAcqEdge ) %in% eipData$Sample_ID[ eipData$State == st ]
  return( isSt )
})

hgtImpCountList = lapply( stIdxList, function( isSt )
{
  return( GetHgtAndImpIsoCounts( eipAcqEdge[ isSt, ], stHgtStatus[ isSt, ] ) )
})

stHgtCounts = lapply( stIdxList,
  function( isSt ) GetTypeCounts( stHgtStatus[ isSt, ] )
  )

hgtCountMat = t( sapply( stHgtCounts, function(x) x[ 1, ] ) )
impCountMat = t( sapply( stHgtCounts, function(x) x[ 2, ] ) )

SubsetCounts = function( inList, idx )
{
  vals = unlist( inList[[ idx ]] )
  return( vals[ vals != 0] )
}


nStHgtEdges = sapply( stIdxList, function(isSt)
{
  nHgtEdges = sapply( 1:ncol( eipAcqEdge ), function(j)
  {
    isHgt = stHgtStatus[ , j ] == HGT_EVENT
    edVals = unique( eipAcqEdge[ isHgt & isSt, j ] )
    edVals = edVals[ edVals != 0 ]
    return( length( edVals ) )
  })
  return( nHgtEdges )
})

hgtStCounts =
  lapply( hgtImpCountList, function(x) SubsetCounts( x, HGT_EVENT ) )
impStCounts =
  lapply( hgtImpCountList, function(x) SubsetCounts( x, IMPORTATION ) )

sapply( hgtStCounts, summary )
sapply( impStCounts, summary )

# ---- Save the data -----------------------------------------------------------

print( pryr::mem_used() )
rData = "../data/2022_02_21_find_cluster_hgt_status.rData"
save( file = rData, list = c( "stHgtStatus", "facHgtStatus" ) )

# ---- Make the figures --------------------------------------------------------

figDir = "../figures/"

typeCounts = GetTypeCounts( stHgtStatus )
clCounts   = colSums( typeCounts )
typeFracs  =
  sapply( 1:ncol( typeCounts ), function(j) typeCounts[ , j] / clCounts[j] )
plOrder    = order( clCounts, decreasing = TRUE )

sum( stHgtStatus  == HGT_EVENT )
sum( stHgtStatus  == IMPORTATION )
par( mar = c( 5.1, 4.1, 4.1, 2.1 ), mfrow = c( 1, 2 ) )
# tiff(
#   filename = paste0( figDir, "2022_02_21_find_cluster_hgt_counts.tiff" ),
#   units    = "in",
#   width    = 6,
#   height   = 4,
#   res      = 300
#   )
barplot(
  typeCounts[ , plOrder ],
  col    = c( "#395983", "#952827" ),
  border = c( "#395983", "#952827" ),
  xlab   = "Gene cluster",
  ylab   = "Frequency",
  main   = "Counts of HGT vs. imported clusters"
  )
legend(
  "topright",
  legend = c( "HGT", "Importation" ),
  col    = c( "#395983", "#952827" ),
  pch    = 15
  )
barplot(
  typeFracs[ , plOrder ],
  col    = c( "#395983", "#952827" ),
  border = c( "#395983", "#952827" ),
  xlab   = "Gene cluster",
  ylab   = "Fraction",
  main   = "Fraction of HGT vs. imported clusters"
  )
# dev.off()

typeCounts = GetTypeCounts( facHgtStatus )
typeFracs  = sapply( 1:ncol( typeCounts ),
  function(j) typeCounts[ , j ] / sum( typeCounts[ , j ]  )
  )

# tiff(
#   filename = paste0( figDir, "2022_02_21_find_cluster_hgt_counts_fac.tiff" ),
#   units    = "in",
#   width    = 6,
#   height   = 4,
#   res      = 300
#   )
barplot(
  typeCounts[ , plOrder ],
  col    = c( "#395983", "#952827" ),
  border = c( "#395983", "#952827" ),
  xlab   = "Gene cluster",
  ylab   = "Frequency",
  main   = "Counts of HGT vs. imported clusters in facilities"
  )
legend(
  "topright",
  legend = c( "HGT", "Importation" ),
  col    = c( "#395983", "#952827" ),
  pch    = 15
  )
barplot(
  typeFracs[ , plOrder ],
  col    = c( "#395983", "#952827" ),
  border = c( "#395983", "#952827" ),
  xlab   = "Gene cluster",
  ylab   = "Fraction",
  main   = "Fraction of HGT vs. imported clusters in facilities"
  )
# dev.off()


xPosHgt = unlist(sapply( 1:length( hgtStCounts ),
  function(i) rep( i, length( hgtStCounts[[i]] )  )
  ))

xPosImp = unlist(sapply( 1:length( impStCounts ),
  function(i) rep( i, length( impStCounts[[i]] )  )
  ))

hgtCols = unlist(sapply( 1:length( hgtStCounts ),
  function(i) rep( annotCols$source[i], length( hgtStCounts[[i]] )  )
  ))

impCols = unlist(sapply( 1:length( impStCounts ),
  function(i) rep( annotCols$source[i], length( impStCounts[[i]] )  )
  ))

par( mar = c( 5.1, 4.1, 4.1, 2.1 ), mfrow = c( 1, 2 ) )
plot(
  x    = jitter( xPosHgt ),
  y    = unlist( hgtStCounts ),
  col  = hgtCols,
  pch  = 16,
  main = "HGT Cluster size",
  xlab = "State",
  ylab = "Number of isolates"
  )
plot(
  x    = jitter( xPosImp ),
  y    = unlist( impStCounts ),
  col  = impCols,
  pch  = 16,
  main = "Imporation Cluster size",
  xlab = "State",
  ylab = "Number of isolates"
  )

nHgt   = sapply( hgtStCounts, length )
nImp   = sapply( impStCounts, length )
stVals = c( "CT", "MN", "TN" )
names( nHgt ) = stVals
names( nImp ) = stVals
barplot(
  nHgt,
  col    = annotCols$source,
  border = annotCols$source,
  main   = "Number of HGT events",
  ylab   = "Number of events"
  )
barplot(
  nImp,
  col    = annotCols$source,
  border = annotCols$source,
  main   = "Number of imporation events",
  ylab   = "Number of events"
  )

stId     =
  sapply( row.names( eipClMat ), function(x) strsplit( x, '_' )[[1]][1] )
stCounts = sapply( stVals, function(x) sum( stId == x ) )
fracHgt  = sapply( seq( 3 ), function(i) nHgt[i] / stCounts[i] )
fracImp  = sapply( seq( 3 ), function(i) nImp[i] / stCounts[i] )

barplot(
  fracHgt,
  col    = annotCols$source,
  border = annotCols$source,
  main   = "Fraction of HGT events\nper genome",
  ylab   = "Fraction of events"
  )
barplot(
  fracImp,
  col    = annotCols$source,
  border = annotCols$source,
  main   = "Fraction of imporation events\nper genome",
  ylab   = "Fraction of events"
  )

sum( )


nStHgtEdges = t(sapply( stIdxList, function(isSt)
{
  nHgtEdges = sapply( 1:ncol( eipAcqEdge ), function(j)
  {
    isHgt = stHgtStatus[ , j ] == HGT_EVENT
    edVals = unique( eipAcqEdge[ isHgt & isSt, j ] )
    edVals = edVals[ edVals != 0 ]
    return( length( edVals ) )
  })
  return( nHgtEdges )
}))
par( mar = c( 5.1, 4.1, 4.1, 2.1 ), mfrow = c( 1, 1 ) )
barplot(
  nStHgtEdges[ , plOrder ],
  col    = annotCols$source,
  border = annotCols$source,
  main   = "HGT Events",
  ylab   = "Fraction of events"

  )

nStEvs = rowSums( nStHgtEdges )

nStIsos = sapply( stIdxList, sum )
nEvsPerIso = sapply( 1:length( nStIsos ), function(i) nStEvs[i] / nStIsos[i] )
names( nEvsPerIso ) = stVals
names( nStEvs ) = stVals
barplot(
  nStEvs,
  col    = annotCols$source,
  border = annotCols$source,
  main   = "Number of HGT Events per state",
  # ylab   = "Fraction of events"
  )
barplot(
  nEvsPerIso,
  col    = annotCols$source,
  border = annotCols$source,
  main   = "Number of HGT Events per state per isolate",
  ylab   = "Fraction of isolates"
  )

barplot(
  nStEvs,
  col    = annotCols$source,
  border = annotCols$source,
  main   = "Number of HGT Events per state",
  # ylab   = "Fraction of events"
  )

par( mar = c( 5.1, 4.1, 4.1, 2.1 ), mfrow = c( 2, 1 ) )
barplot(
  hgtCountMat[ , plOrder ],
  col    = annotCols$source,
  border = annotCols$source,
  main   = "HGT Events",
  ylab   = "Fraction of events"
  )
barplot(
  impCountMat[ , plOrder ],
  col    = annotCols$source,
  border = annotCols$source,
  main   = "Importation",
  ylab   = "Fraction of events"
  )


# ------------------------------------------------------------------------------

sw = c(
  "I: The Phantom Menace"      = 6.5,
  "II: Attack of the Clones"   = 6.5,
  "III: Revenge of the Sith"   = 7.5,
  "IV: A New Hope"             = 8.6,
  "V: The Empire Strikes Back" = 8.7,
  "VI: Return of the Jedi"     = 8.3,
  "VII: The Force Awakens"     = 7.8,
  "VIII: The Last Jedi"        = 6.9,
  "IX: The Rise of Skywalker"  = 6.5
  )

par( mar = c( 12, 4.1, 4.1, 2.1 ), mfrow = c( 1, 1 ) )
barplot( sw, col = "#395983", border = "#395983", las = 2, ylim = c( 0, 10 ) )

# ------------------------------------------------------------------------------
