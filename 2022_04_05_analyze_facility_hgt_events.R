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

isolateMetaData = read.table(
  "../data/isolateMetaData.tsv",
  sep = '\t',
  stringsAsFactors = FALSE,
  header = TRUE
  )
head( isolateMetaData )

# Load the annotations and colors: "annots" and "annotCols"
load( "../data/2022_02_07_make_annotations_and_colors.rData" )

load( "../data/2022_03_09_find_forward_hgt_cluster_transmission.rData" )

# ---- Function definitions ----------------------------------------------------


# ----


eipClMat = eipClMat[ , isSigCl ]

save(
  file = "../data/2022_04_12_eip_hgt_clstr_data.rData",
  list = c( "eipClMat", "hgtGenomeIds", "hgtClIds", "stHgtStatus" )
  )

isSingleHgt = which( clEdgeData$nEdge == 1 & clEdgeData$hgt == HGT_EVENT )
isClExp     = clEdgeData$nEdge > 1

singleHgtIds = sapply( isSingleHgt, function(i)
{
  isClEdge = eipAcqEdge[ , clEdgeData$clId[i] ] == clEdgeData$edgeVals[i]
  return( row.names( eipAcqEdge )[ isClEdge ] )
})

hgtGenomeIds = singleHgtIds
gIdIdxs = sapply( singleHgtIds, function(x) which( eipData$Sample_ID == x ) )



# hgtClIds = colnames( eipAcqEdge )[ clEdgeData$clId[ isSingleHgt ] ]
hgtClIds = clEdgeData$clId[ isSingleHgt ]
HGT_EVENT = 1
IMPORTATION = 2

hgtGenomeIds = unique( hgtGenomeIds )
singleHgtClList = lapply( hgtGenomeIds, function(x)
{
  rIdx = which( row.names( stHgtStatus ) == x )
  return( which( stHgtStatus[ rIdx,  ] == HGT_EVENT ) )
})


sapply( 1:length( hgtGenomeIds ),
  function(i)
    eipClMat[ which( row.names( eipClMat ) == hgtGenomeIds[i] ), hgtClIds[i] ]
  )

# sapply( c( 1, 2, 5, 10,  25 ), function(x)
# {
#   isFac   = !is.na( hgtFacIds ) &
#     hgtFacIds %in% names( facCounts )[ facCounts >= x ]
#   return( sum( isFac ) / length( isFac ) )
# })

facCounts   = table( eipData$Fac_ID )
facVals     = names( facCounts )[ facCounts >= 5 ]
hgtFacIds   = eipData$Fac_ID[ gIdIdxs ]
isFac       = !is.na( hgtFacIds ) & hgtFacIds %in% facVals
hgtFacIds   = hgtFacIds[ isFac ]
gIdIdxs     = gIdIdxs[ isFac ]
isSingleHgt = isSingleHgt[ isFac ]
clIdVals    = unique( clEdgeData$clId[ isSingleHgt ] )

isIn       = eipData$Fac_ID %in% hgtFacIds
eipData    = eipData[ isIn, ]
eipAcqEdge = eipAcqEdge[ row.names( eipAcqEdge ) %in% eipData$Sample_ID, ]

clFacList = lapply( clIdVals, function(x)
{
  isCloneCl = clEdgeData$clId == x & isClExp
  if ( !TRUE %in% isCloneCl ) return( NULL )
  isEdge = eipAcqEdge[ , x ] %in% clEdgeData$edgeVals[ isCloneCl ]
  gIds   = row.names( eipAcqEdge )[ isEdge ]
  facIds = eipData$Fac_ID[ eipData$Sample_ID %in% gIds ]
  names( facIds ) = gIds
  return( facIds )
})

clIdIdx =
  sapply( isSingleHgt, function(i) which( clIdVals == clEdgeData$clId[i] ) )

isCloneFac = sapply( 1:length( isSingleHgt ),
  function(i) sum( clFacList[[ clIdIdx[i] ]] == hgtFacIds[i], na.rm = TRUE )
  )
isCloneFac = isCloneFac > 1
sum( isCloneFac ) / length( isCloneFac )

clEdgeData$clId[ isSingleHgt ]

nClIsos =
  sapply( 1:ncol( eipAcqEdge ), function(j) sum( eipAcqEdge[ , j ] > 0) )


sapply( hgt)

fracLnkd = sapply( clIdVals, function(x)
{
  isClId = clEdgeData$clId[ isSingleHgt ] == x
  return( sum( isCloneFac[ isClId ] ) / sum( isClId ) )
})

nClIsos = nClIsos[ seq( 28 ) %in% clIdVals ]

plot(
  x   = nClIsos,
  y   = fracLnkd,
  col = "#395983",
  pch = 16
  )

xx = table( eipData$Fac_ID )
sapply( seq( 10 ), function(x) sum( xx <= x ) / length( xx ) )

permFacOvlp = sapply( seq( 10000 ), function(x)
{
  rndFacIds = sample( eipData$Fac_ID )
  clFacList = lapply( clIdVals, function(x)
  {
    isCloneCl = clEdgeData$clId == x #& isClExp
    if ( !TRUE %in% isCloneCl ) return( NULL )
    isEdge = eipAcqEdge[ , x ] %in% clEdgeData$edgeVals[ isCloneCl ]
    gIds   = row.names( eipAcqEdge )[ isEdge ]
    facIds = rndFacIds[ eipData$Sample_ID %in% gIds ]
    names( facIds ) = gIds
    return( facIds )
  })

  clIdIdx =
    sapply( isSingleHgt, function(i) which( clIdVals == clEdgeData$clId[i] ) )

  isCloneFac = sapply( 1:length( isSingleHgt ),
    function(i) sum( clFacList[[ clIdIdx[i] ]] == hgtFacIds[i], na.rm = TRUE )
    )
  isCloneFac = isCloneFac > 1
  return( sum( isCloneFac ) / length( isCloneFac ) )
})
hist( permFacOvlp )
summary( permFacOvlp )

facCounts = facCounts[ names( facCounts ) %in%  hgtFacIds ]
isFacIso  = lapply( names( facCounts ), function(x) hgtFacIds == x )
nFacIsos  = sapply( isFacIso, function(x) sum( isCloneFac[x] ) / sum(x) )

cor( cbind( nFacIsos, facCounts ))

xx = glm( nFacIsos ~ facCounts )

facCounts = as.numeric( facCounts )
plot(
  x = facCounts,
  y = nFacIsos,
  col = "#395983",
  pch = 16,
  xlab = "Number of facility isolates",
  ylab = "Fraction of HGT isolates with shared cluster"
  )

cor( facCounts, nFacIsos )

# ------------------------------------------------------------------------------


facIds = sapply( row.names( eipClMat ),
  function(x) eipData$Fac_ID[ eipData$Sample_ID == x ]
  )


hgtGenomeIds = singleHgtIds
gIdIdxs      = sapply( singleHgtIds, function(x) which( eipData$Sample_ID == x ) )
hgtClIds = clEdgeData$clId[ isSingleHgt ]
hgtFacIds

toTest = which( !is.na( hgtFacIds ) )

clHgtFracCounts = sapply( toTest,
  function(i) sum(
    hgtFacIds[i] == facIds[ eipClMat[ , hgtClIds[i] ] == 1 ],
    na.rm = TRUE
    )
  )

clHgtFracCounts



# ---- Save the data -----------------------------------------------------------

print( pryr::mem_used() )
rData = "../data/.rData"
save( file = rData, list = ls() )

# ---- Make the figures --------------------------------------------------------

figDir = "../figures/"
par( mar = c( 5.1, 4.1, 4.1, 2.1 ) )
png(
  filename = paste0( figDir, ".png" ),
  units    = "in",
  width    = 6.5,
  height   = 8,
  res      = 300
  )

dev.off()

hmCols = colorRampPalette( c( "#395983", "#e7e7e7", "#952827" ) )( 100 )
pheatmap(
  ,
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  color             = hmCols, #c( "#e7e7e7", "#395983" ),
  annotation_row    = annots,
  annotation_colors = annotCols,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  filename          = paste0( figDir, ".png" )
  )

# ---- Generate the Rmarkdown --------------------------------------------------

sink( "2022_03_09_analyze_hgt_cluster_distribution.Rmd" )

cat(
"---
title: \"\"
author: \"Ryan D. Crawford\"
date: \"03/09/2022\"
output: pdf_document
---\n\n\n"
)

cat( paste0( "![](", , ")\n\n", collapse = '' ) )

sink()
# ------------------------------------------------------------------------------