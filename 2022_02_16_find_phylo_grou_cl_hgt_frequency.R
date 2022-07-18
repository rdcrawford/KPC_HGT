# ------------------------------------------------------------------------------
#
# 2022/02/
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

load( "../data/eipData.rData" )

# Load the annotations and colors: "annots" and "annotCols"
load( "../data/2022_02_07_make_annotations_and_colors.rData" )

# Load the data on the cluster transition edges: "geneEdgeMat", "clTrData",
# "clMat", "mbrList", "clAcqMat", "clAcqEdgeList", and "clDelEdgeList"
load( "../data/2022_02_08_get_cl_gene_counts_by_edge.rData" )

# ---- Function definitions ----------------------------------------------------

GetGrpEdgeCountMat = function( grpVals, phyloGrps, clAcqMat )
{
  isGrp = lapply( grpVals, function(x) phyloGrps == x )

  edgeCountMat = sapply( 1:ncol( clAcqMat ), function(j)
  {
    sapply( isGrp, function(x)
    {
      edgeVals = unique( clAcqMat[ x, j ] )
      return( sum( edgeVals != 0 ) )
    })
  })
  row.names( edgeCountMat ) = grpVals

  return( edgeCountMat )
}

GetGrpDistMat = function( edgeCountMat )
{
  nGrps   = nrow( edgeCountMat )
  distMat = matrix( 0, nGrps, nGrps )
  row.names( distMat ) = row.names( edgeCountMat )
  colnames( distMat )  = row.names( edgeCountMat )

  clCounts  = rowSums( edgeCountMat )
  for ( i in seq( nGrps ) )
  {
    if ( clCounts[i] > 0 )
    {
      for ( j in seq( nGrps ) )
      {
        distMat[ i, j ] =
          sum( edgeCountMat[ i, edgeCountMat[ j, ] != 0 ] ) / clCounts[i]
      }
    }
  }

  if ( 0 %in% clCounts )
  {
    isHgtGrp = clCounts != 0
    distMat  = distMat[ isHgtGrp, isHgtGrp ]
  }

  return( distMat )
}

# ---- For each phylo-group find the instance of HGT ---------------------------

# Sort the clusters to the same order as the meta-data
mdOrd = sapply( row.names( clAcqMat ),
  function(x) which( isolateMetaData$genomeNames == x )
  )
clAcqMat = clAcqMat[ order( mdOrd ), ]

# Make sure everything is sorted
if ( !identical( row.names( clAcqMat ), isolateMetaData$genomeNames ) )
  stop( "row.names( clAcqMat ) != isolateMetaData$genomeNames" )

spCounts = table( isolateMetaData$species )
spVals   = names( spCounts )[ names( spCounts ) != "Other" ]
stCounts = table( isolateMetaData$mlst )
mlstVals = names( stCounts )[ names( stCounts ) != "Other" &  stCounts > 5 ]
mlstVals = mlstVals[ grepl( "CG", mlstVals ) | grepl( "ST", mlstVals ) ]

isEip = isolateMetaData$genomeNames %in% eipData$Sample_ID

grpClCountMats = list(
  GetGrpEdgeCountMat( spVals,   isolateMetaData$species,          clAcqMat ),
  GetGrpEdgeCountMat( spVals,   isolateMetaData$species[ isEip ], clAcqMat ),
  GetGrpEdgeCountMat( mlstVals, isolateMetaData$mlst,             clAcqMat ),
  GetGrpEdgeCountMat( mlstVals, isolateMetaData$mlst[ isEip ],    clAcqMat )
  )

grpDistMats = lapply( grpClCountMats, GetGrpDistMat )

# ---- Save the data -----------------------------------------------------------

print( pryr::mem_used() )
rData = "../data/2022_02_16_find_phylo_grou_cl_hgt_frequency.rData"
save( file = rData, list = ls() )

# ---- Make the figures --------------------------------------------------------

figDir = "../figures/"
par( mar = c( 5.1, 4.1, 4.1, 2.1 ) )
tiff(
  filename = paste0( figDir, ".tiff" ),
  units    = "in",
  width    = 5,
  height   = 5,
  res      = 300
  )

dev.off()

hmCols = colorRampPalette( c( "#395983", "#e7e7e7", "#952827" ) )( 100 )

for ( i in 1:length( grpDistMats ) )
{
  hCl = agnes( grpDistMats[[i]] )
  pheatmap(
    grpDistMats[[i]][ hCl$order, hCl$order ],
    show_rownames     = TRUE,
    show_colnames     = TRUE,
    color             = hmCols,
    cluster_rows      = FALSE,
    cluster_cols      = FALSE,
    filename          = NA
    )
}


for ( i in 1:length( grpDistMats ) )
{
  mVal = max( grpClCountMats[[i]] )
  hmCols = c( "#e7e7e7", colorRampPalette( c( "#395983", "#952827" ) )( mVal) )
  pheatmap(
    grpClCountMats[[i]],
    show_rownames     = TRUE,
    show_colnames     = TRUE,
    color             = hmCols,
    cluster_rows      = TRUE,
    cluster_cols      = TRUE,
    filename          = NA
    )
}

# ------------------------------------------------------------------------------