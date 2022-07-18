# ------------------------------------------------------------------------------
#
# 2022/03/
# Ryan D. Crawford
# ------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( ape )
library( Rcpp )
library( maps )
library( cognac )
library( future )
library( cluster )
library( phytools )
library( pheatmap )
library( future.apply )

# ---- Constant declarations ---------------------------------------------------


# ---- Import the data ---------------------------------------------------------

# isolateMetaData = read.table(
#   "../data/isolateMetaData.tsv",
#   sep = '\t',
#   stringsAsFactors = FALSE,
#   header = TRUE
#   )

load( "../data/2022_03_09_find_forward_hgt_cluster_transmission.rData" )

#
# # Load the annotations and colors: "annots" and "annotCols"
load( "../data/2022_02_07_make_annotations_and_colors.rData" )


# Load the paths to the alignments and the fastTrees: "clAlgnPaths", and
# "cgTreePaths"
load( "../data/2022_02_25_make_plasmid_trees.rData" )


# ---- Function definitions ----------------------------------------------------

GetDistVals = function(x) return( as.vector( as.dist(x) ) )

GetClData = function( distMat, clIdx )
{

  isSourceList =
    lapply( gIdBySourceList, function(x) row.names( distMat ) %in% x )

  clDistData = rbind.data.frame(
    MakeDistData( distMat, isSourceList[[1]], 1, clIdx ),
    MakeDistData( distMat, isSourceList[[2]], 2, clIdx ),
    MakeDistData( distMat, isSourceList[[3]], 3, clIdx )
    )

  return( clDistData )
}

FmtNames = function( rName ) return( tolower( gsub( ' ', '_', rName ) ) )

MakeDistData = function( distMat, isSource, soIdx, clIdx )
{
  if ( sum( isSource ) <= 1 ) return( NULL )
  soIdxs  = which( isSource )
  soSmry  = sapply( soIdxs,
    function(i) summary( distMat[ i, soIdxs[ soIdxs != i] ] )
    )
  outSmry = sapply( soIdxs, function(i) summary( distMat[ i, !isSource ] ) )

  row.names( soSmry ) =
    paste0( "in_group_", FmtNames( row.names( soSmry ) ) )
  row.names( outSmry ) =
    paste0( "out_group_", FmtNames( row.names( outSmry ) ) )

  rIdxs = sapply( row.names( distMat )[ isSource ],
    function(x) which( row.names( stHgtStatus ) == x )
    )
  hgtStatus = stHgtStatus[ rIdxs, clIdx ]

  strainIdxVals =
    sapply( rIdxs, function(i) which( strainVals == strainId[i] ) )

  soDistData = cbind.data.frame(
    "genome_id"  = row.names( distMat )[ isSource ],
    "cluster_id" = rep( clIdx, length( soIdxs ) ),
    "source"     = rep( soIdx, length( soIdxs ) ),
    "strain_id"  = strainIdxVals,
    "hgt_status" = hgtStatus,
    t( soSmry ),
    t( outSmry )
    )
  # row.names( soDistData ) = row.names( distMat )[ isSource ]
  return( soDistData )
}

# ----

# isSigCl     = names( plMbrList ) %in% colnames( eipClMat )
clAlgnPaths = clAlgnPaths[ isSigCl ]

algnSubMats  = lapply( clAlgnPaths, CreateAlgnDistMat, method = "raw" )
algnDistMats = lapply( clAlgnPaths, CreateAlgnDistMat, method = "shared" )



soVals = c( "EIP_CT", "EIP_MN", "EIP_TN", "Other" )

gIdBySourceList = lapply( soVals,
  function(x) row.names( annots )[ annots$sources == x ]
  )


clDistData = data.frame( )
for ( i in 1:length( algnDistMats ) )
  clDistData = rbind.data.frame( clDistData, GetClData( algnDistMats[[i]], i ) )

hgtCol = "#eaa221"
impCol = "black"
hgtStCols = ifelse( clDistData$hgt_status == 1, hgtCol, impCol )

soCols = sapply( clDistData$source, function(i) annotCols$sources[i] )
plot(
  clDistData$`out_group_median`,
  clDistData$`in group Median`,
  ylim = c( 0, 0.008 ),
  xlim = c( 0, 0.008 ),
  col  = hgtStCols,
  pch  = 16
  )
legend(
  ""
  )
abline( a = 0, b = 1 )

isInMatch = sapply( 1:nrow(clDistData) ,function(i)
  clDistData$`out_group_median`[i] >= clDistData$`in_group_median`[i]
  )
sum( isInMatch ) / length( isInMatch )



library(plotly)
fig <- clDistData %>%
  plot_ly(
    x = ~hgt_status,
    y = ~in_group_mean,
    type = 'violin',
    box = list( visible = TRUE ),
    meanline = list( visible = TRUE),
    x0 = 'Total Bill'
  )
fig


fig = clDistData %>% plot_ly( type = "violin" )

for ( i in seq(3) )
{
  isSource = clDistData$source == i
  fig = fig %>% add_trace(
    x = ~in_group_mean[ isSource ],
    y = ~hgt_status[ isSource ],
    legendgroup = names( annotCols$sources[i] ),
    scalegroup  = names( annotCols$sources[i] ),
    name        = names( annotCols$sources[i] ),
    box         = list( visible = T ),
    meanline    = list( visible = T ),
    color       = I( annotCols$sources[i] )
    )
}
fig

fig <- fig %>%
  add_trace(
    x = ~day[df$sex == 'Female'],
    y = ~total_bill[df$sex == 'Female'],
    legendgroup = 'F',
    scalegroup = 'F',
    name = 'F',
    box = list(
      visible = T
    ),
    meanline = list(
      visible = T
    ),
    color = I("pink")
  )

fig <- fig %>%
  layout(
    yaxis = list(
      zeroline = F
    ),
    violinmode = 'group'
  )

fig = clDistData %>% plot_ly(
  x = ~source,
  y = ~in_group_mean,
  split = ~source,
  type = 'violin',
  box = list( visible = T ),
  meanline = list( visible = T )
  # color    = I( c( "#952827", "#2E4769", "#725971" ) )
  )
fig

fig = clDistData %>% plot_ly(
  x = ~source,
  y = ~out_group_mean,
  split = ~source,
  type = 'violin',
  box = list( visible = T ),
  meanline = list( visible = T )
  # color    = I( c( "#952827", "#2E4769", "#725971" ) )
  )
fig

fig = clDistData %>% plot_ly( type = 'violin' )
fig = fig %>% add_trace(
  x = ~source[ clDistData$hgt_status == 1 ],
  y = ~out_group_median[ clDistData$hgt_status == 1 ],
  legendgroup = "HGT",
  scalegroup = "HGT",
  name = "HGT",
  side = 'negative',
  box = list( visible = T ),
  meanline = list( visible = T ),
  color = I("blue")
  )
fig = fig %>% add_trace(
  x = ~source[ clDistData$hgt_status == 2 ],
  y = ~out_group_median[ clDistData$hgt_status == 2 ],
  legendgroup = "Importation",
  scalegroup = "Importation",
  name = "Importation",
  side = 'positive',
  box = list( visible = T ),
  meanline = list( visible = T ),
  color = I("green")
  )
fig

fig = clDistData %>% plot_ly( type = 'violin' )
fig = fig %>% add_trace(
  x = ~sources,
  y = ~median_in_st,
  legendgroup = "In",
  scalegroup = "In",
  name = "In",
  side = 'negative',
  # box = list( visible = T ),
  # meanline = list( visible = T ),
  color = I("#395983")
  )
fig = fig %>% add_trace(
  x = ~sources,
  y = ~mean_out_st,
  legendgroup = "Out",
  scalegroup = "Out",
  name = "Out",
  side = 'positive',
  # box = list( visible = T ),
  # meanline = list( visible = T ),
  color = I("#952827")
  )
fig

fig <- fig %>%
  add_trace(
    x = ~day[df$smoker == 'No'],
    y = ~total_bill[df$smoker == 'No'],
    legendgroup = 'No',
    scalegroup = 'No',
    name = 'No',
    side = 'positive',
    box = list(
      visible = T
    ),
    meanline = list(
      visible = T
    ),
    color = I("green")
  )



inVsOutDataL = t(lapply( seq(3), function(x)
{
  isSo = clDistData$source == x
  mVals = c(
    mean( clDistData$in_group_mean[ isSo ] ),
    mean( clDistData$out_group_mean[ isSo ] )
    )
  sdVals = c(
    sd( clDistData$in_group_mean[ isSo ] ),
    sd( clDistData$out_group_mean[ isSo ] )
    )
  cbind.data.frame(
    rep( soVals[x], 2 ),
    c( "in-group", "out-group" ),
    mVals,
    sdVals
    )
}))
inVsOutData = inVsOutDataL[[1]]
for ( i in seq( 2, 3 ) )
  inVsOutData = rbind.data.frame( inVsOutData, inVsOutDataL[[i]])

colnames( inVsOutData ) = c( "source", "group", "mean", "sd" )


inVsOutDataL = t(lapply( seq(3), function(x)
{
  isSo = clDistData$source == x
  mVals = c(
    mean( clDistData$in_group_mean[ isSo ] ),
    mean( clDistData$out_group_mean[ isSo ] )
    )
  sdVals = c(
    sd( clDistData$in_group_mean[ isSo ] ),
    sd( clDistData$out_group_mean[ isSo ] )
    )
  cbind.data.frame(
    rep( soVals[x], 2 ),
    c( "in-group", "out-group" ),
    mVals,
    sdVals
    )
}))
inVsOutData = inVsOutDataL[[1]]
for ( i in seq( 2, 3 ) )
  inVsOutData = rbind.data.frame( inVsOutData, inVsOutDataL[[i]])

inVsOutData$source = as.factor( inVsOutData$source )
inVsOutData$group = as.factor( inVsOutData$group )

fig <- plot_ly(
  data = inVsOutData[ inVsOutData$group == "in-group", ],
  x = ~source,
  y = ~mean,
  type = 'bar',
  name = "in-group"
  # error_y = ~list(array = sd, color = '#000000')
  )
fig <- fig %>% add_trace(
  data = inVsOutData[ inVsOutData$group == "out-group", ],
  name = "out-group"
  )
fig


fig = clDistData %>% plot_ly( type = 'bar' )
fig = fig %>% add_trace(
  x = ~source,
  y = ~in_group_median,
  legendgroup = "in_group_median",
  scalegroup = "in_group_median",
  name = "in_group_median",
  # side = 'negative',
  # box = list( visible = T ),
  # meanline = list( visible = T ),
  color = I("blue")
  )
fig = fig %>% add_trace(
  x = ~source,
  y = ~out_group_median,
  legendgroup = "out_group_median",
  scalegroup = "out_group_median",
  name = "out_group_median",
  # side = 'positive',
  # box = list( visible = T ),
  # meanline = list( visible = T ),
  color = I("green")
  )
fig


fig <- plot_ly(data, x = ~Animals, y = ~SF_Zoo, type = 'bar', name = 'SF Zoo')
fig <- fig %>% add_trace(y = ~LA_Zoo, name = 'LA Zoo')
fig <- fig %>% layout(yaxis = list(title = 'Count'), barmode = 'group')

fig

# ---- Save the data -----------------------------------------------------------

print( pryr::mem_used() )
rData = "../data/2022_04_14_compare_plasmid_alignment_distances.rData"
save( file = rData, list = ls() )

# ---- Make the figures --------------------------------------------------------

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

# ---- Generate the Rmarkdown --------------------------------------------------

# sink( "2022_03_09_analyze_hgt_cluster_distribution.Rmd" )
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