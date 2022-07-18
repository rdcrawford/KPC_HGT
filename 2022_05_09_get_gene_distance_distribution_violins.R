# ------------------------------------------------------------------------------
#
# 2022/03/
# Ryan D. Crawford
# ------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( ape )
library( plotly )
library( scales )
library( cluster )
library( magrittr )
library( pheatmap )

# ---- Import the data ---------------------------------------------------------

load( "../data/2022_04_14_compare_plasmid_alignment_distances.rData" )

isolateMetaData = read.table(
  "../data/isolateMetaData.tsv",
  sep = '\t',
  stringsAsFactors = FALSE,
  header = TRUE
  )

# Load the annotations and colors: "annots" and "annotCols"
load( "../data/2022_02_07_make_annotations_and_colors.rData" )

load( "../data/2022_02_21_find_cluster_hgt_status.rData" )

# ---- Constant declarations ---------------------------------------------------

# Threshold for
MIN_IN_GRP = 1

HGT_EVENT   = 1
IMPORTATION = 2

# ---- Function definitions ----------------------------------------------------

# clIdx = 1; distMat = algnDistMats[[clIdx]];
GetClData = function( distMat, clIdx )
{
  rIdxs = sapply( row.names( distMat ),
    function(x) which( row.names( metaData ) == x )
    )

  sourceId = metaData$Source[ rIdxs ]
  facId    = metaData$facility[ rIdxs ]
  toTest   = which( sourceId != "Other" )

  distData = data.frame( )
  for ( i in toTest )
  {
    distData = rbind.data.frame(
      distData, GetClStats( distMat[i, ], i, sourceId, facId, clIdx )
      )
  }

  return( distData )
}

# i = 6; distVals = distMat[i, ];
GetClStats = function( distVals, i, sourceId, facId, clIdx )
{
  qVals = seq( length( distVals ) )[ -c(i) ]
  nQ    = length( qVals )
  isSo  = sourceId[ qVals ] == sourceId[i]

  if ( is.na( facId[i] ) )
  {
    isFac = rep( NA,nQ )
  } else {
    isFac = sapply( facId[ qVals ],
      function(x) ifelse( is.na(x), FALSE, x ==  facId[i] )
      )

    if ( TRUE %in% isFac )
    {
      isFac = isFac & isSo
    } else {
      isFac = rep( NA, nQ )
    }
  }


  qDistData = data.frame(
    "gId"    = rep( names( distVals )[i], nQ ),
    "clIdx"  = rep( clIdx, nQ ),
    "source" = rep( sourceId[i], nQ ),
    "dist"   = distVals[ qVals ],
    "isSo"   = isSo,
    "isFac"  = isFac
    )
  return( qDistData )
}

GetStrainId = function( annots )
{
  strainId = sapply( 1:nrow( annots ), function(i)
  {
    if ( annots$ST[i] != "Other" ) return( annots$ST[i] )
    return( annots$Species[i] )
  })
  return( strainId )
}

GetFacId = function( annots )
{
  facIds = sapply( 1:nrow( annots ), function(i)
  {
    if ( annots$Source[i] == "Other" ) return( NA )
    return( eipData$Fac_ID[ eipData$Sample_ID == row.names( annots )[i] ] )
  })
  return( facIds )
}

# ----

# Create a list of lists where the isolates decending from each cluster
# acquisiton edge are created for each cluster
clGenomesByAcqEdge = lapply( 1:ncol( eipAcqEdge ), function(j)
{
  edgeVals = unique( eipAcqEdge[ eipAcqEdge[ , j ] != 0, j ] )
  isoList  = lapply( edgeVals,
    function(x) row.names( eipAcqEdge )[ eipAcqEdge[ , j ] == x ]
    )
  return( isoList )
})
names( clGenomesByAcqEdge ) = colnames( eipAcqEdge )

isDupIsoList = lapply( clGenomesByAcqEdge, function(x)
{
  patId = lapply( x,
    function(y) eipData$Patient_ID[ eipData$Sample_ID %in% y ]
    )
  return( lapply( patId, duplicated ) )
})

clGenomesByAcqEdge = sapply( 1:ncol( eipAcqEdge ), function(j)
{
  sapply( 1:length( clGenomesByAcqEdge[[j]] ),
    function(i) clGenomesByAcqEdge[[j]][[i]][ !isDupIsoList[[j]][[i]] ]
    )
})
names( clGenomesByAcqEdge ) = colnames( eipAcqEdge )

edgeCounts = unlist(lapply( clGenomesByAcqEdge, function(x)
{
  nIsos = sapply( x, length )
  return( unlist( lapply( nIsos, function(x) rep(x, x)) ) )
}))

clIdxs = unlist( lapply( 1:ncol( eipAcqEdge ),
  function(j) lapply( sapply( clGenomesByAcqEdge[[j]], length ),
    function(x) rep( j, x )
    )
  ))
stHgtStatus =
  stHgtStatus[ , colnames(stHgtStatus  ) %in% colnames( eipAcqEdge ) ]

isHgt = unlist(lapply( 1:ncol( stHgtStatus ), function(j)
{
  hgtIsos = row.names( stHgtStatus )[ stHgtStatus[ , j ] == HGT_EVENT ]
  return( unlist( clGenomesByAcqEdge[[j]] ) %in% hgtIsos )
}))


genomeClData = cbind.data.frame(
  "genomeIds"  = unlist( clGenomesByAcqEdge ),
  "clIdxs"     = clIdxs,
  "edgeCounts" = edgeCounts,
  "clNames"    = colnames( eipAcqEdge )[ clIdxs ],
  "isHgt"      = isHgt
  )

sum( genomeClData$edgeCounts == 1 )
sum( genomeClData$isHgt == 1 )
sum( genomeClData$isHgt == 1 & genomeClData$edgeCounts == 1 )

save(
  file = "../data/clGenomesByAcqEdge.rData",
  list = c( "clGenomesByAcqEdge", "genomeClData" )
  )

# ----


metaData = cbind.data.frame(
  annots, "strain" = GetStrainId( annots ), "facility" = GetFacId( annots )
  )

for ( i in 1:length( algnDistMats ) )
{
  isEip = row.names( algnDistMats[[i]] ) %in% eipData$Sample_ID
  isUni = row.names( algnDistMats[[i]] ) %in% unlist( clGenomesByAcqEdge[[i]] )
  isIn  = isUni | !isEip
  algnDistMats[[i]] = algnDistMats[[i]][ isIn, isIn ]
}

clDistData = data.frame( )
for ( i in 1:length( algnDistMats ) )
  clDistData = rbind.data.frame( clDistData, GetClData( algnDistMats[[i]], i ) )

notMissingFac = unique( clDistData$gId[ !is.na( clDistData$isFac ) ] )


hgtClstrData =
  genomeClData[ genomeClData$edgeCounts == 1 & genomeClData$isHgt, ]

ii =  hgtClstrData$genomeIds %in% notMissingFac
soIds = sapply( hgtClstrData$genomeIds[ii], function(x) strsplit(x, '_')[[1]][1] )
soIds = soIds[ !duplicated( hgtClstrData$genomeIds[ii] ) ]


isSngHgt = sapply( 1:nrow( clDistData ), function(i)
{
  isSmpl  = hgtClstrData$genomeIds == clDistData$gId[i]
  isClstr = hgtClstrData$clIdxs == clDistData$clIdx[i]
  return( TRUE %in% ( isSmpl & isClstr ) )
})
sum( isSngHgt )
sngHgtData = clDistData[ isSngHgt, ]


# ---- Do the statistical analysis by state

stVals = c( "CT", "MN", "TN" )

stSilRnkSumTest = lapply( stVals,
  function(x) wilcox.test( dist ~ isSo, sngHgtData[ sngHgtData$source == x, ] )
  )
sapply( stSilRnkSumTest, function(x) x$p.value )

isFac = !is.na( sngHgtData$isFac )
facWilRnkSumTest = lapply( stVals, function(x)
{
  wilcox.test( dist ~ isSo, sngHgtData[ sngHgtData$source == x & isFac, ] )
})
sapply( facWilRnkSumTest, function(x) x$p.value )

wilcox.test( dist ~ isSo, sngHgtData )

stVals = c( "CT", "MN", "TN" )

stSilRnkSumTest = lapply( stVals,
  function(x) wilcox.test( dist ~ isSo, sngHgtData[ sngHgtData$source == x, ] )
  )
sapply( stSilRnkSumTest, function(x) x$p.value )

stKsTest = lapply( stVals, function(x)
{
  ks.test(
    sngHgtData$dist[ sngHgtData$source == x & sngHgtData$isSo ],
    sngHgtData$dist[ sngHgtData$source == x & !sngHgtData$isSo ]
    )
})
sapply( stKsTest, function(x) x$p.value )

isFac = !is.na( sngHgtData$isFac )
facKsTest = lapply( stVals, function(x)
{
  ks.test(
    sngHgtData$dist[ sngHgtData$source == x & isFac & sngHgtData$isFac ],
    sngHgtData$dist[ sngHgtData$source == x & isFac & !sngHgtData$isFac ]
    )
})
sapply( facKsTest, function(x) x$p.value )

ks.test(
  sngHgtData$dist[ sngHgtData$isSo ],
  sngHgtData$dist[ !sngHgtData$isSo ]
  )

# ----
soVals = sapply( sngHgtData$gId, function(x) strsplit( x, '_' )[[1]][1] )

stClList = sapply( c( "CT", "MN", "TN" ),
  function(x) unique( sngHgtData$clIdx[ soVals == x ] )
  )

stClMat = sapply( stClList,
  function(x) sapply( stClList, function(y) sum( x %in% y ) )
  )

nCols =  max( stClMat ) - min( stClMat )
hmCols = colorRampPalette( c( "#395983", "#952827" ) )( nCols )

pheatmap(
  stClMat,
  show_rownames   = TRUE,
  show_colnames   = TRUE,
  color           = hmCols,
  number_format   = "%.0f",
  number_color    = "white",
  display_numbers = TRUE,
  cluster_rows    = TRUE,
  cluster_cols    = TRUE,
  filename        = NA#paste0( figDir, ".png" )
  )
stCls = unique( unlist( stClList ) )

stClMat = sapply( stClList,
  function(x) sapply( stCls, function(y) sum( x %in% y ) )
  )
pheatmap(
  stClMat,
  show_rownames   = TRUE,
  show_colnames   = TRUE,
  color           = c( "#e7e7e7", "#395983" ),
  # number_format   = "%.0f",
  # number_color    = "white",
  # display_numbers = TRUE,
  cluster_rows    = TRUE,
  cluster_cols    = TRUE,
  filename        = NA#paste0( figDir, ".png" )
  )


f = function( clDistData )
{
  soVals = c( "CT", "MN", "TN") # , "Other" )
  medVals = lapply( soVals, function(x)
  {
     isSoIso = clDistData$source == x
     medVals = c(
       median( clDistData$dist[ clDistData$isSo & isSoIso ] ),
       median( clDistData$dist[ !clDistData$isSo & isSoIso ] )
       )
     y = gsub( 'EIP_', '', x)
     names( medVals ) =
       c( paste0( y, " in state" ), paste0( y, " out of state" ) )
     return( medVals )
  })

  xx = unlist( medVals )

  plCols = c( "#952827", "#952827", "#2E4769", "#2E4769", "#725971", "#725971" )
  par( mar = c( 7.5, 4.1, 4.1, 2.1 ) )
  barplot( xx, las = 2, col = plCols, border = plCols )
}
f( clDistData )
f( sngHgtData )

wilcox.test(
  clDistData$dist[ clDistData$isFac & isFacData ],
  clDistData$dist[ !clDistData$isFac & isFacData ]
  )

wilcox.test(
  clDistData$dist[ clDistData$isSo ],
  clDistData$dist[ !clDistData$isSo ]
  )

# ---- Save the data -----------------------------------------------------------

print( pryr::mem_used() )
rData = "../data/2022_05_09_get_gene_distance_distribution_violins.rData"
save( file = rData, list = ls() )

# ---- Make the figures --------------------------------------------------------

figDir = "../figures/2022_05_09_get_gene_distance_distribution_violins/"

yAxisLayOut = list(
  title = "HGT Singleton alignment distance",
  range = c( 0, 0.012 )
  )
xAxisLayOut = list( title = "State" )

fig = plot_ly() %>%
  add_trace(
    data        = clDistData[ clDistData$isSo, ],
    x           = ~source,
    y           = ~dist,
    legendgroup = "Within state",
    scalegroup  = "Within state",
    name        = "Within state",
    side        = 'negative',
    type        = 'violin',
    color       = I("#395983"),
    box         = list( visible =  TRUE ),
    points      = TRUE,
    spanmode    = "hard"
    ) %>%
  add_trace(
    data        = clDistData[ !clDistData$isSo, ],
    x           = ~source,
    y           = ~dist,
    legendgroup = "Outside state",
    scalegroup  = "Outside state",
    name        = "Outside state",
    side        = 'positive',
    type        = 'violin',
    color       = I("#952827"),
    box         = list( visible = TRUE ),
    points      = TRUE,
    spanmode    = "hard"
    ) %>%
  layout(
    yaxis = yAxisLayOut,
    xaxis = xAxisLayOut
    )
fig
orca( fig, paste0( figDir, "distance_distribution_by_state.png" ) )


isFacData = !is.na( clDistData$isFac )
fig = plot_ly() %>%
  add_trace(
    data        = clDistData[ clDistData$isFac & isFacData, ],
    x           = ~source,
    y           = ~dist,
    legendgroup = "Within facility",
    scalegroup  = "Within facility",
    name        = "Within facility",
    side        = 'negative',
    type        = 'violin',
    color       = I("#395983"),
    box         = list( visible =  TRUE ),
    points      = TRUE,
    spanmode    = "hard"
    ) %>%
  add_trace(
    data        = clDistData[ !clDistData$isFac & isFacData, ],
    x           = ~source,
    y           = ~dist,
    legendgroup = "Outside facility, within state",
    scalegroup  = "Outside facility, within state",
    name        = "Outside facility, within state",
    side        = 'positive',
    type        = 'violin',
    color       = I("#952827"),
    box         = list( visible = TRUE ),
    points      = TRUE,
    spanmode    = "hard"
    ) %>%
  layout(
    yaxis = yAxisLayOut,
    xaxis = xAxisLayOut
    )
fig
orca( fig, paste0( figDir, "distance_distribution_by_facility.png" ) )

# ------------------------------------------------------------------------------

fig = plot_ly() %>%
  add_trace(
    data        = sngHgtData[ sngHgtData$isSo, ],
    x           = ~source,
    y           = ~dist,
    legendgroup = "Within state",
    scalegroup  = "Within state",
    name        = "Within state",
    side        = 'negative',
    type        = 'violin',
    color       = I("#395983"),
    box         = list( visible =  TRUE ),
    points      = TRUE,
    spanmode    = "hard"
    ) %>%
  add_trace(
    data        = sngHgtData[ !sngHgtData$isSo, ],
    x           = ~source,
    y           = ~dist,
    legendgroup = "Outside state",
    scalegroup  = "Outside state",
    name        = "Outside state",
    side        = 'positive',
    type        = 'violin',
    color       = I("#952827"),
    box         = list( visible = TRUE ),
    points      = TRUE,
    spanmode    = "hard"
    ) %>%
  layout(
    yaxis = yAxisLayOut,
    xaxis = xAxisLayOut
    )
fig
orca( fig, paste0( figDir, "distance_distribution_by_state_hgt_sing.png" ) )

isFacData = !is.na( sngHgtData$isFac )
wilcox.test(
  sngHgtData$dist[ sngHgtData$isFac & isFacData ],
  sngHgtData$dist[ !sngHgtData$isFac & isFacData ]
  )
fig = plot_ly() %>%
  add_trace(
    data        = sngHgtData[ sngHgtData$isFac & isFacData, ],
    x           = ~source,
    y           = ~dist,
    legendgroup = "Within facility",
    scalegroup  = "Within facility",
    name        = "Within facility",
    side        = 'negative',
    type        = 'violin',
    color       = I("#395983"),
    box         = list( visible =  TRUE ),
    points      = TRUE,
    spanmode    = "hard"
    ) %>%
  add_trace(
    data        = sngHgtData[ !sngHgtData$isFac & isFacData, ],
    x           = ~source,
    y           = ~dist,
    legendgroup = "Outside facility, within state",
    scalegroup  = "Outside facility, within state",
    name        = "Outside facility, within state",
    side        = 'positive',
    type        = 'violin',
    color       = I("#952827"),
    box         = list( visible = TRUE ),
    points      = TRUE,
    spanmode    = "hard"
    ) %>%
  layout(
    yaxis = yAxisLayOut,
    xaxis = xAxisLayOut
    )
fig
orca( fig, paste0( figDir, "distance_distribution_by_facility_hgt_sing.png" ) )

wilcox.test(
  sngHgtData$dist[ sngHgtData$isFac & isFacData ],
  sngHgtData$dist[ !sngHgtData$isFac & isFacData ]
  )
# ------------------------------------------------------------------------------