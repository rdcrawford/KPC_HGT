# ------------------------------------------------------------------------------
#
# 2022/03/09
# Ryan D. Crawford
# ------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------

# ---- Constant declarations ---------------------------------------------------

HGT_EVENT   = 1
IMPORTATION = 2

# ---- Import the data ---------------------------------------------------------

# Load the data on the deduplicated eip gene clusters: "eipClMat",
# "eipAcqEdge", "plMbrList", and "plMat"
load( "../data/2022_02_16_analyze_hgt_cluster_distribution.rData" )

# Load the data on the KPC acquisition permutation test: "pValsTrEdge",
# "sigThresh", and "idxs"
load( "../data/2022_02_09_find_cluster_kpc_acquisition_enrichment_sig.rData" )

# Load the annotations and colors: "annots" and "annotCols"
load( "../data/2022_02_07_make_annotations_and_colors.rData" )

# Load the data on the deduplicated eip gene clusters: "stHgtStatus" and
# "facHgtStatus"
load( "../data/2022_02_21_find_cluster_hgt_status.rData" )

# ---- Function definitions ----------------------------------------------------

GetStrainId = function( annots )
{
  strainId = sapply( 1:nrow( annots ), function(i)
  {
    if ( annots$ST[i] != "Other" ) return( annots$ST[i] )
    return( annots$species[i] )
  })
  return( strainId )
}

GetIdxCounts = function( inMat, idxList, val )
{
  countMat = sapply( 1:ncol( inMat ),
    function(j) sapply( idxList,
      function(i) sum( inMat[ i, j ] == val, na.rm = TRUE )
      )
    )
  return( countMat )
}

NormalizeToCounts = function( inMat, counts )
{
  return( sapply( 1:ncol( inMat), function(j) inMat[  , j ] / counts[j] ) )
}

MakeHgtPlots = function( hgtStatusData, rIdxs, vals, cols, mnTxt, figHandle )
{
  # Make matricies with the counts of importation and hgt by strain
  hgtMat = GetIdxCounts( hgtStatusData, rIdxs, HGT_EVENT )
  impMat = GetIdxCounts( hgtStatusData, rIdxs, IMPORTATION )

  strHgtCounts = rowSums( hgtMat )
  strImpCounts = rowSums( impMat )
  names( strHgtCounts ) = vals
  names( strImpCounts ) = vals

  # Get the number of acquisitions for each cluster
  nAcqEvents = sapply( 1:ncol( eipAcqEdge ), function(j)
  {
    nStrEdges = sapply( rIdxs, function(i)
    {
      isHgt =  hgtStatusData[ i , j ] == HGT_EVENT
      return( length( unique( eipAcqEdge[ i[ isHgt ], j ] ) ) )
    })
    return( nStrEdges )
  })

  nVal = sapply( rIdxs, length )
  strHgtCountsNorm =
    sapply( 1:length( strHgtCounts ), function(i) strHgtCounts[i] / nVal[i] )
  strImpCountsNorm =
    sapply( 1:length( strHgtCounts ), function(i) strImpCounts[i] / nVal[i] )

  nAcqs = rowSums( nAcqEvents )
  nAcqsNorm = sapply( 1:length( nAcqs ), function(i) nAcqs[i] / nVal[i] )
  names( nAcqs )     = vals
  names( nAcqsNorm ) = vals

  figNames = c(
    "hgt_vs_importation_barplot",
    "hgt_vs_importation_iso_counts",
    "hgt_vs_importation_iso_counts_per_isolate",
    "number_of_acquisition_events_by_cluster",
    "acquisition_event_count",
    "acquisition_events_per_isolate",
    "acquisition_events_by_isolates"
    )
  figPaths = paste0( figDir, "2022_03_09_", figHandle, '_', figNames, ".png" )

  png(
    filename = figPaths[1],
    units    = "in",
    width    = 6,
    height   = 9.5,
    res      = 300
    )
  plOrder = order( colSums( hgtMat ), decreasing = TRUE )
  par( mar = c( 2, 4.1, 2, 2.1 ), mfrow = c( 2, 1 ) )
  barplot(
    hgtMat[ , plOrder ],
    col    = cols,
    border = cols,
    ylab   = "Number of isolates",
    xlab   = "Gene cluster",
    main   = paste0( "Count of isolates", mnTxt, " HGT" )
    )
  legend(
    "topright",
    legend = vals,
    col    = cols,
    pch    = 15,
    cex    = 0.75
    )
  barplot(
    impMat[ , plOrder ],
    col = cols,
    border = cols,
    ylab   = "Number of isolates",
    xlab   = "Gene cluster",
    main   = paste0( "Number of ", mnTxt, " imporation isolates" )
    )
  dev.off()

  png(
    filename = figPaths[2],
    units    = "in",
    width    = 6,
    height   = 9.5,
    res      = 300
    )
  par( mar = c( 10, 4.1, 2, 2.1 ), mfrow = c( 2, 1 ) )
  barplot(
    strHgtCounts,
    col    = cols,
    border = cols,
    las    = 2,
    ylab   = "Number of isolates",
    main   = paste0( "Number of ", mnTxt, " HGT isolates" )
    )
  barplot(
    strImpCounts,
    col    = cols,
    border = cols,
    las    = 2,
    ylab   = "Number of isolates",
    main   = paste0( "Number of ", mnTxt, " imporation isolates" )
    )
  dev.off()

  png(
    filename = figPaths[3],
    units    = "in",
    width    = 6,
    height   = 9.5,
    res      = 300
    )
  par( mar = c( 10, 4.1, 2, 2.1 ), mfrow = c( 2, 1 ) )
  barplot(
    strHgtCountsNorm,
    col    = cols,
    border = cols,
    las    = 2,
    ylab   = "Events per isolate",
    main   = paste0( "Number of ", mnTxt, " HGT isolates" )
    )
  barplot(
    strImpCountsNorm,
    col    = cols,
    border = cols,
    las    = 2,
    ylab   = "Events per isolate",
    main   = paste0( "Number of ", mnTxt, " imporation isolates" )
    )
  dev.off()

  plOrd = order( colSums( nAcqEvents ), decreasing = TRUE )

  png(
    filename = figPaths[4],
    units    = "in",
    width    = 5,
    height   = 5,
    res      = 300
    )
  par( mar = c( 4.1, 4.1, 4.1, 2.1 ), mfrow = c( 1, 1 ) )
  barplot(
    nAcqEvents[ , plOrd ],
    col    = cols,
    border = cols,
    las    = 2,
    xlab   = "Gene cluster",
    ylab   = "Number of cluster acquisitions",
    main   = paste0( "Number of ", mnTxt, " acquisition events" )
    )
  legend(
    "topright",
    legend = vals,
    col    = cols,
    pch    = 15,
    cex    = 0.75
    )
  dev.off()

  png(
    filename = figPaths[5],
    units    = "in",
    width    = 5,
    height   = 5,
    res      = 300
    )
  par( mar = c( 10, 4.1, 4.1, 2.1 ), mfrow = c( 1, 1 ) )
  barplot(
    nAcqs,
    col    = cols,
    border = cols,
    las    = 2,
    ylab   = "Number of cluster acquisitions",
    main   = paste0( "Number of ", mnTxt, " acquisition events" )
    )
  dev.off()

  png(
    filename = figPaths[6],
    units    = "in",
    width    = 5,
    height   = 5,
    res      = 300
    )
  par( mar = c( 10, 4.1, 4.1, 2.1 ), mfrow = c( 1, 1 ) )
  barplot(
    nAcqsNorm,
    col    = cols,
    border = cols,
    las    = 2,
    ylab   = "Number of cluster acquisitions per isolate",
    main   = paste0( "Number of ", mnTxt, " acquisition events" )
    )
  dev.off()

  png(
    filename = figPaths[7],
    units    = "in",
    width    = 5,
    height   = 5,
    res      = 300
    )
  par( mar = c( 5.1, 4.1, 4.1, 2.1 ), mfrow = c( 1, 1 ) )
  plot(
    nAcqsNorm,
    nVal,
    col  = cols,
    pch  = 16,
    cex  = 2,
    xlab = "Number of acquisitons observed",
    ylab = "Number of isolates "
    )
  dev.off()

  return( figPaths )
}

# ---- Find the KPC acquisition enriched clusters ------------------------------

# Make a list of the gene clusters which were merged
clList = lapply( names( plMbrList ),
  function(x) as.integer( strsplit( x, '/' )[[1]] )
  )

# Find the KPC enriched clusters
isSig   = pValsTrEdge <= sigThresh
sigIdxs = idxs[ isSig ]
isSigCl = sapply( clList, function(x) TRUE %in% ( x %in% sigIdxs ) )

# Sort the acquisition edge matrix to match the annotations
annots = annots[ row.names( annots ) %in% row.names( stHgtStatus ), ]
aOrd = sapply( row.names( annots ),
  function(x) which( row.names( facHgtStatus ) == x )
  )
annots = annots[ order( aOrd ), ]


# Subset to the clusters significantly associated with KPC
eipAcqEdge   = eipAcqEdge[ , isSigCl ]
eipClMat     = eipClMat[ , isSigCl ]
facHgtStatus = facHgtStatus[ , isSigCl ]
stHgtStatus  = stHgtStatus[ , isSigCl ]

dim( annots )
dim( stHgtStatus )
dim( facHgtStatus )
isIdent = c(
  identical( row.names( stHgtStatus ), row.names( facHgtStatus ) ),
  identical( row.names( eipAcqEdge ),  row.names( facHgtStatus ) ),
  identical( row.names( annots ),      row.names( facHgtStatus ) ),
  identical( row.names( annots ),      row.names( eipAcqEdge ) ),
  identical( row.names( annots ),      row.names( eipClMat ) ),
  identical( row.names( annots ),      row.names( stHgtStatus ) )
  )

if ( FALSE %in% isIdent )
  stop( "There data are not sorted..." )

# Get the data on the strains in the study
strainId   = GetStrainId( annots )
strainVals = unique( strainId )
strainVals = strainVals[ strainVals != "Other" ]
strainVals = strainVals[ order( strainVals ) ]
strainIdxs = lapply( strainVals, function(x) which( strainId == x ) )
colVals    = c( annotCols$ST, annotCols$species[ 1:5 ] )
strainCols = colVals[ order( names( colVals ) ) ]
strainCols = strainCols[ names( strainCols ) != "Other" ]

# Get the data on the states in the study
stVals = names( annotCols$sources )[ 1:3 ]
stCols = annotCols$sources[ 1:3 ]
stIdxs = lapply( stVals, function(x) which( annots$sources == x ) )

# ---- Save the data -----------------------------------------------------------

print( pryr::mem_used() )
rData = "../data/2022_03_09_analyze_hgt_cluster_distribution.rData"
save( file = rData, list = ls() )

# ---- Make the figures --------------------------------------------------------

figDir = "../figures/2022_03_09_analyze_hgt_cluster_distribution/"

stateByStrain = MakeHgtPlots(
  stHgtStatus, strainIdxs, strainVals, strainCols, "state", "state_by_strain"
  )
stateByState = MakeHgtPlots(
  stHgtStatus, stIdxs, stVals, stCols, "state", "state_by_state"
  )
strainByFacility = MakeHgtPlots(
  facHgtStatus, strainIdxs, strainVals, strainCols, "facility",
  "strain_by_facility"
  )
stateByFacility = MakeHgtPlots(
  facHgtStatus, stIdxs, stVals, stCols, "facility", "state_by_facility"
  )

# ---- Generate the Rmarkdown --------------------------------------------------

sink( "2022_03_09_analyze_hgt_cluster_distribution.Rmd" )

cat(
"---
title: \"Analyze HGT cluster distribution\"
author: \"Ryan D. Crawford\"
date: \"03/09/2022\"
output: pdf_document
---\n\n\n"
)

cat( "# Distribution of HGT events in state by strain\n" )
cat( paste0( "![](", stateByStrain, ")\n\n", collapse = '' ) )
cat( "# Distribution of HGT events in state by state\n" )
cat( paste0( "![](", stateByState, ")\n\n", collapse = '' ) )
cat( "# Distribution of HGT events in facility by strain\n" )
cat( paste0( "![](", strainByFacility, ")\n\n", collapse = '' ) )
cat( "# Distribution of HGT events in facility by state\n" )
cat( paste0( "![](", stateByFacility, ")\n\n", collapse = '' ) )

sink()

# ------------------------------------------------------------------------------