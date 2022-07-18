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

# ---- Constant declarations ---------------------------------------------------

# Threshold for
MIN_IN_GRP = 2

# ---- Function definitions ----------------------------------------------------

GetDistVals = function(x) return( as.vector( as.dist(x) ) )

GetClData = function( distMat, clIdx )
{
  rIdxs = sapply( row.names( distMat ),
    function(x) which( row.names( metaData ) == x )
    )

  sourceId = metaData$sources[ rIdxs ]
  facId    = metaData$facility[ rIdxs ]
  toTest   = which( sourceId != "Other" )

  distData = t( sapply( toTest, function(i)
  {
    isSo = sourceId == sourceId[i]
    soStats = GetDistStats( distMat[ i, ], isSo )
    if ( is.na( facId[i] ) )
    {
      facStats = rep( NA, 4 )
    } else {

      isFac = sapply( facId,
        function(x) ifelse( is.na(x), FALSE, x == facId[i] )
        )
      isFac = isFac & isSo
      facStats = GetDistStats( distMat[ i, -c(i) ], isFac[ -c(i) ] )
    }
    return( c( soStats, facStats ) )
  }))


  distData = cbind.data.frame(
    row.names( distMat )[ toTest ], sourceId[ toTest ], facId[ toTest ],
    rep( clIdx, length( toTest ) ), distData
    )

  colnames( distData ) = c(
    "sample_id", "sources", "facility", "cl_idx", "mean_in_st", "median_in_st",
    "mean_out_st", "median_out_st", "mean_in_fac", "median_in_fac",
    "mean_out_fac", "median_out_fac"
    )

  return( distData )
}


GetDistStats = function( distVals, isGrp )
{
  if ( sum( isGrp ) < MIN_IN_GRP ) return( rep( NA, 4 )  )

  distStats = c(
    mean( distVals[ isGrp ] ),
    median( distVals[ isGrp ] ),
    mean( distVals[ !isGrp ]),
    median( distVals[ !isGrp ] )
    )
  return( distStats )
}

GetStrainId = function( annots )
{
  strainId = sapply( 1:nrow( annots ), function(i)
  {
    if ( annots$ST[i] != "Other" ) return( annots$ST[i] )
    return( annots$species[i] )
  })
  return( strainId )
}

GetFacId = function( annots )
{
  facIds = sapply( 1:nrow( annots ), function(i)
  {
    if ( annots$sources[i] == "Other" ) return( NA )
    return( eipData$Fac_ID[ eipData$Sample_ID == row.names( annots)[i] ] )
  })
  return( facIds )
}

# ----

# metaData = cbind.data.frame(
#   annots, "strain" = GetStrainId( annots ), "facility" = GetFacId( annots )
#   )
#
# # isSigCl     = names( plMbrList ) %in% colnames( eipClMat )
# clAlgnPaths = clAlgnPaths[ isSigCl ]
#
# algnSubMats  = lapply( clAlgnPaths, CreateAlgnDistMat, method = "raw" )
# algnDistMats = lapply( clAlgnPaths, CreateAlgnDistMat, method = "shared" )



soVals = c( "EIP_CT", "EIP_MN", "EIP_TN", "Other" )

gIdBySourceList = lapply( soVals,
  function(x) row.names( annots )[ annots$sources == x ]
  )

clDistData = data.frame( )
for ( i in 1:length( algnDistMats ) )
  clDistData = rbind.data.frame( clDistData, GetClData( algnDistMats[[i]], i ) )

isSngHgt = sapply( 1:nrow( clDistData ), function(i)
{
  isSmpl  = hgtGenomeIds == clDistData$sample_id[i]
  isClstr = hgtClIds == clDistData$cl_idx[i]
  return( TRUE %in% ( isSmpl & isClstr ) )
})


# ---- Save the data -----------------------------------------------------------

print( pryr::mem_used() )
rData = "../data/.rData"
save( file = rData, list = ls() )

# ---- Make the figures --------------------------------------------------------

figDir = "../figures/2022_04_22_analyze_cluster_distance_distributions/"


isData = !is.na( clDistData$mean_in_st)

sum( isData )
fig = plot_ly() %>%
  add_trace(
    data        = clDistData[ isData, ],
    x           = ~sources,
    y           = ~median_in_st,
    legendgroup = "In",
    scalegroup  = "In",
    name        = "In",
    side        = 'negative',
    type        = 'violin',
    color       = I("#395983"),
    box         = list( visible =  TRUE ),
    points      = TRUE,
    spanmode    = "hard"
    ) %>%
  add_trace(
    data        = clDistData[ isData, ],
    x           = ~sources,
    y           = ~median_out_st,
    legendgroup = "Out",
    scalegroup  = "Out",
    name        = "Out",
    side        = 'positive',
    type        = 'violin',
    color       = I("#952827"),
    box         = list( visible = TRUE ),
    points      = TRUE,
    spanmode    = "hard"
    ) %>%
  layout(
    title = "Distance distribution by state",
    yaxis = list( title = "Median distance" ),
    xaxis = list( title = "State" )
    )
fig
orca( fig, paste0( figDir, "distance_distribution_by_state.png" ) )

isData = !is.na( clDistData$median_in_st ) & isSngHgt

fig = plot_ly() %>%
  add_trace(
    data        = clDistData[ isData, ],
    x           = ~sources,
    y           = ~median_in_st,
    legendgroup = "In",
    scalegroup  = "In",
    name        = "In",
    side        = 'negative',
    type        = 'violin',
    color       = I("#395983"),
    box         = list( visible =  TRUE ),
    points      = TRUE,
    spanmode    = "hard"
    ) %>%
  add_trace(
    data        = clDistData[ isData, ],
    x           = ~sources,
    y           = ~median_out_st,
    legendgroup = "Out",
    scalegroup  = "Out",
    name        = "Out",
    side        = 'positive',
    type        = 'violin',
    color       = I("#952827"),
    box         = list( visible = TRUE ),
    points      = TRUE,
    spanmode    = "hard"
    ) %>%
  layout(
    title = "HGT singleton distance distribution by state",
    yaxis = list( title = "Median distance" ),
    xaxis = list( title = "State" )
    )
fig
orca( fig, paste0( figDir, "hgt_singleton_by_state.png" ) )

# ------------------------------------------------------------------------------

isData = !is.na( clDistData$median_in_fac)
table( clDistData$sources[ isData ] )
fig = plot_ly() %>%
  add_trace(
    data        = clDistData[ isData, ],
    x           = ~sources,
    y           = ~median_in_fac,
    legendgroup = "In",
    scalegroup  = "In",
    name        = "In",
    side        = 'negative',
    type        = 'violin',
    color       = I("#395983"),
    box         = list( visible =  TRUE ),
    points      = TRUE,
    spanmode    = "hard"
    ) %>%
  add_trace(
    data        = clDistData[ isData, ],
    x           = ~sources,
    y           = ~median_out_fac,
    legendgroup = "Out",
    scalegroup  = "Out",
    name        = "Out",
    side        = 'positive',
    type        = 'violin',
    color       = I("#952827"),
    box         = list( visible = TRUE ),
    points      = TRUE,
    spanmode    = "hard"
    ) %>%
  layout(
    title = "Distance distribution by facility",
    yaxis = list( title = "Median distance" ),
    xaxis = list( title = "State" )
    )
fig
orca( fig, paste0( figDir, "distance_distribution_by_facility.png" ) )

isData = !is.na( clDistData$median_in_fac) & isSngHgt

table( clDistData$sources[ isData ] )

fig = plot_ly() %>%
  add_trace(
    data        = clDistData[ isData, ],
    x           = ~sources,
    y           = ~median_in_fac,
    legendgroup = "In",
    scalegroup  = "In",
    name        = "In",
    side        = 'negative',
    type        = 'violin',
    color       = I("#395983"),
    box         = list( visible =  TRUE ),
    points      = TRUE,
    spanmode    = "hard"
    ) %>%
  add_trace(
    data        = clDistData[ isData, ],
    x           = ~sources,
    y           = ~median_out_fac,
    legendgroup = "Out",
    scalegroup  = "Out",
    name        = "Out",
    side        = 'positive',
    type        = 'violin',
    color       = I("#952827"),
    box         = list( visible = TRUE ),
    points      = TRUE,
    spanmode    = "hard"
    ) %>%
  layout(
    title = "HGT singleton distance distribution by facility",
    yaxis = list( title = "Median distance" ),
    xaxis = list( title = "State" )
    )
fig
orca( fig, paste0( figDir, "distance_distribution_by_facility.png" ) )

# ------------------------------------------------------------------------------
#
# fig = clDistData %>% plot_ly( type = 'violin' )
# fig = fig %>% add_trace(
#   x = ~sources,
#   y = ~median_in_st,
#   legendgroup = "In",
#   scalegroup = "In",
#   name = "In",
#   side = 'negative',
#   # box = list( visible = T ),
#   # meanline = list( visible = T ),
#   color = I("#395983")
#   )
# fig = fig %>% add_trace(
#   x = ~sources,
#   y = ~mean_out_st,
#   legendgroup = "Out",
#   scalegroup = "Out",
#   name = "Out",
#   side = 'positive',
#   # box = list( visible = T ),
#   # meanline = list( visible = T ),
#   color = I("#952827")
#   )
# fig
#
# clDistData = clDistData[ hgtIsoIdxs, ]
# fig = clDistData %>% plot_ly( type = 'violin' )
#
#
# soSmy = data.frame()
# for ( x in  c( "EIP_CT", "EIP_MN", "EIP_TN" ) )
# {
#   isSo = clDistData$sources == x
#   soSmy = rbind.data.frame(
#     soSmy,
#     c( x, "in", quantile(  clDistData$median_in_st[ isSo ], na.rm = TRUE ) ),
#     c( x, "out", quantile(  clDistData$median_out_st[ isSo ], na.rm = TRUE ) )
#     )
# }
# colnames( soSmy ) = c( "source", "dist", "Q1", "Q2", "Q3", "Q4", "Q5" )
# soSmy$source = as.factor( soSmy$source )
# soSmy$dist   = as.factor( soSmy$dist )
# for ( j in seq( 3, ncol( soSmy ) ) )  soSmy[ , j ] = as.numeric( soSmy[ , j ] )
#
#
# cols = unlist( lapply( annotCols$sources[ 1:3 ], function(x) rep( x, 2 ) ) )
#
# medianVals = soSmy$Q3
# names( medianVals ) = paste( soSmy$source, soSmy$dist )
#
# par( mar = c( 8.1, 4.1, 4.1, 2.1 ) )
# barplot(
#   medianVals,
#   col = cols,
#   border = cols,
#   las = 2
#   )
#
#
# fig = plot_ly( type = "box" )
#
# fig = fig %>% add_trace(
#   data        = soSmy[ soSmy$dist == "in", ],
#   x           = ~source,
#   legendgroup = "In",
#   lowerfence  = ~Q1,
#   q1          = ~Q2,
#   median      = ~Q2,
#   q3          = ~Q4,
#   upperfence  = ~Q5,
#   color       = I( "#952827" )
#   )
#
# fig = fig %>% add_trace(
#   data        = soSmy[ soSmy$dist == "out", ],
#   x           = ~source,
#   legendgroup = "Out",
#   lowerfence  = ~Q1,
#   q1          = ~Q2,
#   median      = ~Q2,
#   q3          = ~Q4,
#   upperfence  = ~Q5,
#   color       = I( "#395983" )
#   )
# fig
# #
# # fig = fig %>% layout(
# #   yaxis = list( exponentformat = "SI",type="log",title = "Median distance"),
# #   xaxis = list( title = "State" ),
# #   boxmode = "group"
# #   )
#
# fig = fig %>% layout(
#    boxmode = "group",
#    yaxis   = list( title = "Median distance" ),
#    xaxis   = list( title = "State" )
#    )
# fig
#
# fig = fig %>% add_trace(
#   x           = ~sources,
#   y           = ~median_in_st,
#   legendgroup = "HGT",
#   scalegroup  = "HGT",
#   name        = "HGT",
#   side        = "negative",
#   box         = list( visible = T ),
#   meanline    = list( visible = T ),
#   color       = I( "#952827" )
#   )
# fig = fig %>% add_trace(
#   x           = ~sources,
#   y           = ~median_out_st,
#   legendgroup = "Median out of state distance",
#   scalegroup  = "Median out of state distance",
#   name        = "Median out of state distance",
#   side        = "positive",
#   box         = list( visible = T ),
#   meanline    = list( visible = T ),
#   color       = I( "#395983" )
#   )
#
# fig
#
# fig = clDistData %>% plot_ly( type = "violin" )
#
# for ( i in seq(3) )
# {
#   isSource = clDistData$source == i
#   fig = fig %>% add_trace(
#     x = ~in_group_mean[ isSource ],
#     y = ~hgt_status[ isSource ],
#     legendgroup = names( annotCols$sources[i] ),
#     scalegroup  = names( annotCols$sources[i] ),
#     name        = names( annotCols$sources[i] ),
#     box         = list( visible = T ),
#     meanline    = list( visible = T ),
#     color       = I( annotCols$sources[i] )
#     )
# }
# fig
#
# figDir = "../figures/"
# png(
#   filename = paste0( figDir, ".png" ),
#   units    = "in",
#   width    = 6.5,
#   height   = 8,
#   res      = 300
#   )
# par( mar = c( 5.1, 4.1, 4.1, 2.1 ) )
#
# dev.off()
#
#
# xx= sapply( 1:nrow( clDistData ),
#   function(i) clDistData$median_in_fac[i] <= clDistData$median_out_fac[i]
#   )
# sum( xx, na.rm = TRUE ) / sum( !is.na( xx ) )
#
# xx= sapply( 1:nrow( clDistData ),
#   function(i) clDistData$median_in_st[i] <= clDistData$median_out_st[i]
#   )
# sum( xx, na.rm = TRUE ) / sum( !is.na( xx ) )
#
#
# hmCols = colorRampPalette( c( "#395983", "#e7e7e7", "#952827" ) )( 100 )
# pheatmap(
#   ,
#   show_rownames     = FALSE,
#   show_colnames     = FALSE,
#   color             = c( "#e7e7e7", "#395983" ),#hmCols,
#   annotation_row    = annots,
#   annotation_colors = annotCols,
#   cluster_rows      = TRUE,
#   cluster_cols      = TRUE,
#   filename          = NA#paste0( figDir, ".png" )
#   )
#
# # ----
#
#
# hgtIsoIdxs = sapply( 1:length( hgtGenomeIds ),
#   function(i) which(
#     clDistData$sample_id == hgtGenomeIds[i] &
#     clDistData$cl_idx == hgtClIds[i]
#     )
#   )
#
# nHgt = length( hgtIsoIdxs )
# plotCols = sapply( hgtIsoIdxs, function(i) annotCols$sources[ clDistData$sources[i] ] )
# clIdxs = c( 5, 7, 9, 11 )
# distVals = lapply( clIdxs, function(j) clDistData[ hgtIsoIdxs, j ] )
# xPos     =
#   unlist( lapply( seq( length( clIdxs ) ), function(i) rep( i, nHgt ) ) )
# xRng     = c( 0.1, 0.9 )
# xPos     = rescale( xPos, to = c( xRng[1], xRng[2] ) )
# maxVal   = 0.01
# yInc     = 0.001
# nTicks   = maxVal / yInc
# addVal   = 1e-5
# baseVal  = 2
# yPos     = log( c( unlist( distVals ), maxVal ) + addVal )
# yPos     = rescale( yPos, c( 0, 1 ) )
# yPos     = yPos[ 1:( length( yPos ) - 1 ) ]
# axisVals = seq( 0, maxVal, yInc )
# yAxisPos = rescale(
#   log( axisVals  + addVal ),
#   to = c( 0, 1 )
#   )
# sapply( 2:length(yAxisPos), function(i) yAxisPos[i] - yAxisPos[i-1])
# grpLabs = c(
#   "Within state", "Out of state", "Within facility", "Out of facility"
#   )
#
# plot.new()
# points(
#   x    = jitter( xPos ),
#   y    = jitter( yPos ),
#   col  = cols,
#   pch  = 16,
#   # xlab = "State",
#   ylab = "Median distance"
#   )
# axis(
#   side   = 1,
#   at     = seq( xRng[1], xRng[2], length.out = 4 ),
#   labels = grpLabs,
#   las    = 2
#   )
# axis(
#   side   = 2,
#   at     = yAxisPos, #seq( 0, 1, length.out = nTicks + 1 ),
#   labels = axisVals, #seq( 0, yInc * nTicks, yInc ),
#   las    = 2
#   )
#
# # ---- Generate the Rmarkdown --------------------------------------------------
#
# sink( ".Rmd" )
#
# cat(
# "---
# title: \"\"
# author: \"Ryan D. Crawford\"
# date: \"03/09/2022\"
# output: pdf_document
# ---\n\n\n"
# )
#
# cat( paste0( "![](", , ")\n\n", collapse = '' ) )
#
# sink()

# ------------------------------------------------------------------------------