# ------------------------------------------------------------------------------
#
# 2022/05/
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

# ---- Import the data ---------------------------------------------------------

isolateMetaData = read.table(
  "../data/isolateMetaData.tsv",
  sep = '\t',
  stringsAsFactors = FALSE,
  header = TRUE
  )

# Load the meta-data from the EIP isolates
load( "../data/eipData.rData" )

# Load the annotations and colors: "annots" and "annotCols"
load( "../data/2022_02_07_make_annotations_and_colors.rData" )

# Load the data on the plasmid clusters: "eipClMat", "eipAcqEdge", "plMbrList",
# and "plMat"
load( "../data/2022_02_16_analyze_hgt_cluster_distribution.rData" )

# Load the data on whether cluster were acquired prior to importation to state
# or facility: "stHgtStatus", and "facHgtStatus"
load( "../data/2022_02_21_find_cluster_hgt_status.rData" )

# Load the mpr results on the edge in which an isolate was included in a
# region. Includes: "arFeatureMats", "edgeColMat", "featEdgeMat",
# "edgeFeatList", and "tipTrEdgeMat"
load( "../data/2022_02_19_run_mpr_for_location.rData" )

# "geneEdgeMat", "clTrData", "clMat", "mbrList", "clAcqMat",
# "clAcqEdgeList", and "clDelEdgeList"
load( "../data/2022_02_08_get_cl_gene_counts_by_edge.rData" )

# Read in the cognac fastTree
cgTree = midpoint.root( read.tree( "../data/cgTree.tre" ) )


# Load the annotations and colors: "annots" and "annotCols"
load( "../data/2022_02_07_make_annotations_and_colors.rData" )


# ---- Constant declarations ---------------------------------------------------

HGT_EVENT   = 1
IMPORTATION = 2

# ---- Function definitions ----------------------------------------------------



# ----

annots = annots[ row.names( annots ) %in% row.names( stHgtStatus ), ]
aOrd =  sapply( row.names( annots ),
  function(x) which( row.names( stHgtStatus ) == x)
  )
annots = annots[ order( aOrd ), ]
eipData = eipData[ eipData$Sample_ID %in% row.names( stHgtStatus ), ]
eipAcqEdge  = eipAcqEdge[ , isSigCl ]

identical( dimnames( eipAcqEdge ),  dimnames( stHgtStatus ) )
identical( row.names( eipAcqEdge ),  row.names( annots ) )
# Get the data on isolates who have a patient assigned to them in the meta-data
isMissing = is.na( eipData$Patient_ID )
patId     = eipData$Patient_ID[ !isMissing ]
names( patId  ) = eipData$Sample_ID[ !isMissing ]

strainId = GetStrainId( annots )
# Make a list of the
patVals   = unique( patId )
patIsos   = lapply( patVals, function(x) names( patId  )[ patId == x ] )
patIsoIdx = lapply( patIsos,
  function(x) sapply( x, function(y) which( row.names( stHgtStatus ) == y ) )
  )
nPatIsos  = sapply( patIsos, length )
patIsos   = patIsos[ nPatIsos > 1 ]
patIsoIdx = patIsoIdx[ nPatIsos > 1 ]

isoHgtStatus = sapply( 1:nrow( stHgtStatus ), function(i)
{
  isImp = TRUE %in% ( stHgtStatus[ i, ] == IMPORTATION )
  if ( isImp ) return( "Importation" )
  isHgt = TRUE %in% ( stHgtStatus[ i, ] == HGT_EVENT )
  if ( isHgt ) return( "HGT" )
  return( "Unknown" )
})


patHgtStatus = sapply( patIsoIdx, function(x) isoHgtStatus[x] )
patStrVals   = sapply( patIsoIdx, function(x) strainId[x] )

xx = sapply( 1:length( patHgtStatus ), function(i)
{
  isHgt = patHgtStatus[[i]] == "HGT"
  isImp = patHgtStatus[[i]] == "Importation"
  if ( TRUE %in% isHgt )
  {
    if ( sum( isHgt ) == 1 ) return( "Multiple isolates" )
    patHgtStrs = patStrVals[[i]][ patHgtStatus[[i]] == "HGT" ]
    strVals =  unique( patHgtStrs )
    nn = sapply( strVals, function(x) sum( patHgtStrs == x ) ) > 1
    if ( TRUE %in% nn ) return( "Multiple HGT" )
    return( "Multiple isolates" )
  }
  if ( TRUE %in% isImp ) return( "Importation" )
  return( "Undetermined" )
})

yy = table( xx )
yy / sum( yy )

patStatVals = sapply( patHgtStatus, function(x) unique(x) )
nStats = sapply( patStatVals, length )
isOnStat = nStats == 1
counts = c( sum( isOnStat ), sum( !isOnStat ) )
counts / sum( counts )


sapply( 1:length( isOnStat ), function(i)
{

})
unlist( isOnStat, patStatVals[  ] )
isoCounts = table( isoHgtStatus )
isoCounts / sum( isoCounts )
patStatCounts = table(  )
patStatCounts / sum( patStatCounts )

strCounts = table( strainId )
strCounts / sum( strCounts)
# ---- Save the data -----------------------------------------------------------

print( pryr::mem_used() )
save(
  file = "../data/.rData",
  list = ls()
  )

# ---- Make the figures --------------------------------------------------------

hist(
  nPatIsos,
  main = "Patient isolate counts",
  xlab = "Number of isolates per patient",
  col  = "#395983",
  border = "#395983",
  breaks = seq( max( nPatIsos ) )
  )

barplot(
  yy,
  main = "Patient isolate counts",
  # xlab = "Number of isolates per patient",
  col  = "#395983",
  border = "#395983",
  ylim = c( 0, 50)
  # breaks = seq( max( nPatIsos ) )
  )


figDir = "../figures/"
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
  color             = c( "#e7e7e7", "#395983" ),#hmCols,
  annotation_row    = annots,
  annotation_colors = annotCols,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  filename          = NA#paste0( figDir, ".png" )
  )

# ---- Generate the Rmarkdown --------------------------------------------------

sink( ".Rmd" )

cat(
"---
title: \"\"
author: \"Ryan D. Crawford\"
date: \"05//2022\"
output: pdf_document
---\n\n\n"
)

cat( paste0( "![](", figPaths, ")\n\n", collapse = '' ) )

sink()

# ------------------------------------------------------------------------------