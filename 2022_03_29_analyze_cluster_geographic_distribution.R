# ------------------------------------------------------------------------------
#
# 2022/03/
# Ryan D. Crawford
# ------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( ape )
library( maps )
library( cluster )
library( phytools )
library( pheatmap )

# ---- Constant declarations ---------------------------------------------------


# ---- Import the data ---------------------------------------------------------

isolateMetaData = read.table(
  "../data/isolateMetaData.tsv",
  sep = '\t',
  stringsAsFactors = FALSE,
  header = TRUE
  )

# "geneEdgeMat", "clTrData", "clMat", "mbrList", "clAcqMat", "clAcqEdgeList",
# and "clDelEdgeList"
load( "../data/2022_02_08_get_cl_gene_counts_by_edge.rData" )

# Load the annotations and colors: "annots" and "annotCols"
load( "../data/2022_02_07_make_annotations_and_colors.rData" )

# Load the data on the deduplicated eip gene clusters: "eipClMat",
# "eipAcqEdge", "plMbrList", and "plMat"
load( "../data/2022_02_16_analyze_hgt_cluster_distribution.rData" )

# Load the data on the KPC acquisition permutation test: "pValsTrEdge",
# "sigThresh", and "idxs"
load( "../data/2022_02_09_find_cluster_kpc_acquisition_enrichment_sig.rData" )

# ---- Function definitions ----------------------------------------------------

source( "R/TreeParsingFunctions.R" )

# ---- Find the KPC acquisition enriched clusters ------------------------------

# Make a list of the gene clusters which were merged
clList = lapply( names( plMbrList ),
  function(x) as.integer( strsplit( x, '/' )[[1]] )
  )

# Find the KPC enriched clusters
isSig   = pValsTrEdge <= sigThresh
sigIdxs = idxs[ isSig ]
isSigCl = sapply( clList, function(x) TRUE %in% ( x %in% sigIdxs ) )
clList  = clList[ isSigCl ]

plAcqMat = sapply( clList, function(x)
{
  sapply( 1:nrow( clAcqMat ), function(i)
  {
    isAcq = clAcqMat[ i, x ] != 0
    if ( TRUE %in% isAcq ) return( clAcqMat[ i, x ][ isAcq ][1] )
    return(0)
  })
})
row.names( plAcqMat ) = row.names( clAcqMat )
colnames( plAcqMat )  = names( plMbrList )[ isSigCl ]

plMat = sapply( 1:ncol( plAcqMat ), function(j) plAcqMat[ , j ] != 0 )
dimnames( plMat ) = dimnames( plAcqMat )

mdOrd = sapply( row.names( plAcqMat ),
  function(x) which( isolateMetaData$genomeNames == x )
  )
plAcqMat = plAcqMat[ order( mdOrd ), ]
plMat    = plMat[ order( mdOrd ), ]
identical( row.names( plAcqMat ), isolateMetaData$genomeNames )

isEip   = grepl( "EIP", isolateMetaData$sources )
pubIdxs = which( !isEip )


pubClEdgeList = lapply( 1:ncol( plAcqMat ), function(j)
{
  eipAcqEdges = unique( plAcqMat[ isEip, j] )
  eipAcqEdges = eipAcqEdges[ eipAcqEdges != 0 ]
  sharedEdgeList = lapply( eipAcqEdges, function(x)
  {
    isEdge = plAcqMat[ pubIdxs, j ] == x
    if ( TRUE %in% isEdge ) return( pubIdxs[ isEdge ] )
    return( NULL )
  })
  return( sharedEdgeList )
})


MakeRegionMat = function( pubClEdgeList, cIdx )
{
  regionVals = unique( isolateMetaData[ , cIdx ] )
  regionVals = regionVals[ !is.na( regionVals ) ]
  regionVals = regionVals[ order( regionVals ) ]
  regionMat  = matrix( 0, length( regionVals ), length( pubClEdgeList ) )

  for ( j in 1:length( pubClEdgeList ) )
  {
    rIdxs = which(
      regionVals %in% isolateMetaData[ unlist( pubClEdgeList[[j]] ), cIdx ]
      )
    for ( i in rIdxs ) regionMat[ i, j ] = 1
  }

  row.names( regionMat ) = regionVals

  return( regionMat )
}

regionMats =
  lapply( c( 11, 12, 13 ), function(c) MakeRegionMat( pubClEdgeList, c ) )


rHcl = lapply( regionMats, function(x) agnes( x ) )
cHcl = lapply( regionMats, function(x) agnes( t( x ) ) )

sapply( regionMats, dim )
sapply( rHcl, function(x) length( x$order ) )
sapply( cHcl, function(x) length( x$order ) )

for ( i in 1:length( regionMats ) )
{
  pheatmap(
    regionMats[[i]],
    show_rownames     = TRUE,
    show_colnames     = FALSE,
    color             = c( "#e7e7e7", "#395983" ),
    cluster_rows      = TRUE,
    cluster_cols      = TRUE,
    filename          = NA #paste0( figDir, ".png" )
    )
}

eipStVals = c( "EIP_CT", "EIP_MN", "EIP_TN" )
stIdxs = lapply( eipStVals, function(x) which( isolateMetaData$sources == x ) )
stRegionMatList = lapply( 1:length( eipStVals ), function(ii)
{
  pubClEdgeList = lapply( 1:ncol( plAcqMat ), function(j)
  {
    eipAcqEdges = unique( plAcqMat[ stIdxs[[ii]], j] )
    eipAcqEdges = eipAcqEdges[ eipAcqEdges != 0 ]
    sharedEdgeList = lapply( eipAcqEdges, function(x)
    {
      isEdge = plAcqMat[ pubIdxs, j ] == x
      if ( TRUE %in% isEdge ) return( pubIdxs[ isEdge ] )
      return( NULL )
    })
    return( sharedEdgeList )
  })

  cIdxs =  c( 11, 12, 13 )
  regionMats = lapply( cIdxs, function(c) MakeRegionMat( pubClEdgeList, c ) )



  # for ( i in 1:length( regionMats ) )
  # {
  #   hmTitle = paste(
  #     gsub( '_', ' ', eipStVals[ii] ),
  #     "importation isolates by",
  #     colnames( isolateMetaData)[ cIdxs[i] ]
  #     )
  #
  #   pheatmap(
  #     regionMats[[i]][ rHcl[[i]]$order,  cHcl[[i]]$order ],
  #     show_rownames     = TRUE,
  #     show_colnames     = FALSE,
  #     color             = c( "#e7e7e7", annotCols$sources[ii] ),
  #     cluster_rows      = FALSE,
  #     cluster_cols      = FALSE,
  #     main              = hmTitle,
  #     filename          = NA #paste0( figDir, ".png" )
  #     )
  # }

  return( regionMats )
})


featVals = annots$sources
identical( row.names( annots )[ isEip ], row.names( plMat )[ isEip ] )

featVals = annots$sources

eipClMat = plAcqMat[ isEip, ]
featVals = strainId
MakeFeatMat = function( plAcqMat, featVals, stIdxs )
{

  fCats = unique( featVals )
  fCats = fCats[ order( fCats ) ]
  fMat  = sapply( 1:ncol( plAcqMat ), function(j)
  {
    edgeVals = unique( plAcqMat[ stIdxs, j ] )
    edgeVals = edgeVals[ edgeVals != 0 ]
    isCl = ( plAcqMat[ , j ] %in% edgeVals ) & !isEip
    return( as.factor( as.integer( fCats %in% featVals[ isCl ] )  ) )
    # return( as.character( fCats %in% featVals[ isCl ] ) )
  })

  fMat = data.frame( t( fMat ) )
  colnames( fMat )  = fCats
  row.names( fMat ) = colnames( plAcqMat )
  # print( rowSums( fMat ) )
  return(  fMat)
}


hmCols = lapply( strainCols, function(x) c( `1` = x, `0` = "white" ) )

colAnnots = MakeFeatMat( plAcqMat, strainId )[ , 1:9 ]

colnames( colAnnots ) = names( hmCols )

for ( i in 1:length( regionMats ) )
  colnames( regionMats[[i]] ) = colnames(eipClMat  )

mat = regionMats[[i]][ rHcl[[i]]$order,  cHcl[[i]]$order ]

for ( ii in seq(3) )
{
  colAnnots = MakeFeatMat( plAcqMat, featVals, stIdxs[[ii]] )[ , 1:9 ]
  mat =  stRegionMatList[[ii]][[3]][ rHcl[[i]]$order,  cHcl[[i]]$order ]
  colnames( mat ) = colnames(eipClMat  )
  pheatmap(
    mat,
    show_rownames     = TRUE,
    show_colnames     = FALSE,
    color             = c( "#e7e7e7", annotCols$sources[ii] ), # "#395983"
    cluster_rows      = FALSE,
    cluster_cols      = FALSE,
    annotation_col    = colAnnots,
    annotation_colors = hmCols,
    annotation_legend = FALSE,
    legend            = FALSE,
    filename          = NA #paste0( figDir, ".png" )
    )
}


identical( colnames( colAnnots ), names( hmCols ) )

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
