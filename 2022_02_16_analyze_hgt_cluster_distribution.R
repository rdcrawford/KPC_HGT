# ------------------------------------------------------------------------------
# Analyze HGT cluster distribution
# 2022/02/16
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# Deduplicate the clusters, which are present in a high proportion of the same
# isolates and occur on the same transition edge, likely representing the
# same plasmids. Then investigate the distribution of clusters across the
# different species and STs
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( ape )
library( maps )
library( cluster )
library( phytools )
library( pheatmap )

# ---- Constant declarations ---------------------------------------------------

# Fraction of shared genomes for a cluster to be identical
GENOME_FRAC = 0.9

# ---- Import the data ---------------------------------------------------------

# Load the data on the EIP genomes
load( "../data/eipData.rData" )

# Load the data on the cluster transtion edges: "geneEdgeMat", "clTrData",
# "clMat", "mbrList", "clAcqMat", "clAcqEdgeList", and "clDelEdgeList"
load( "../data/2022_02_08_get_cl_gene_counts_by_edge.rData" )

# Load the annotations and colors: "annots" and "annotCols"
load( "../data/2022_02_07_make_annotations_and_colors.rData" )

# ---- Function definitions ----------------------------------------------------

GetDupliacteCls = function( eipClMat, eipAcqEdge )
{
  # Find the genomes which have each cluster
  clList =
    lapply( 1:ncol( eipClMat ), function(j) which( eipClMat[ , j ] == 1 ) )

  # Calcuate the fraction of shared genomes shared by each cluster
  clSimMat = sapply( clList,
    function(x) sapply( clList, function(y) sum( x %in% y ) / length(x) )
    )

  clIdxs = vector( "integer", length( clList ) )
  ref    = 1
  clId   = 1
  while( 0 %in% clIdxs )
  {
    clIdxs[ ref ] = clId
    qrys   = which( clIdxs == 0 )
    nGs    = length( clList[[ ref ]] )
    percCl = sapply( qrys,
      function(i) max( c( clSimMat[ i, ref ], clSimMat[ ref, i ] ) )
      )
    isOverlap =  percCl >= GENOME_FRAC
    if ( TRUE %in% isOverlap )
    {
      isAcqEdge = FindSameAcqEdgeCls( eipAcqEdge, ref, qrys, clList )
      isPl      = isOverlap & isAcqEdge
      for ( i in qrys[ isPl ] ) clIdxs[i] = clId
      ref  = qrys[ !isPl ][1]
    } else {
      ref  = qrys[1]
    }
    clId = clId + 1
  }
  return( lapply( seq( max( clIdxs ) ), function(x) which( clIdxs == x ) ) )
}

FindSameAcqEdgeCls = function( eipAcqEdge, ref, qrys, clList )
{
  idxs = unique( unlist( clList[ c( ref, qrys ) ] ) )
  isAcqEdge = sapply( qrys, function(j)
  {
    isQryIdx = sapply( idxs, function(i)
    {
      if ( eipAcqEdge[ i, j ] == 0 || eipAcqEdge[ i, ref ] == 0 ) return( TRUE )
      return( eipAcqEdge[ i, j ] == eipAcqEdge[ i, ref ] )
    })
    return( !FALSE %in% isQryIdx )
  })
  return( isAcqEdge )
}

MakeFeaturMat = function( clMat, factorVals, cIdx )
{
  plFacVals = sapply( 1:ncol( clMat ), function(j)
  {
    isG  = row.names( annots ) %in% row.names( clMat )[ clMat[ , j ] == 1 ]
    return( sapply( factorVals, function(x) sum( annots[ isG, cIdx ] == x ) ) )
  })
  return( plFacVals )
}

source( "R/GetTreeOrder.R" )

# ---- Deduplicate the the gene clusters ---------------------------------------

# Subset the data to the EIP genomes and clusters
isEip      = row.names( clMat ) %in% eipData$Sample_ID
eipClMat   = clMat[ isEip,  ]
isEipCl    = colSums( eipClMat ) > 2
eipClMat   = eipClMat[ , isEipCl ]
mbrList    = mbrList[ isEipCl ]
eipAcqEdge = clAcqMat[ isEip ,isEipCl ]
clMat      = clMat[ , isEipCl ]

# Find the clusters which are duplicated
clList   = GetDupliacteCls( eipClMat, eipAcqEdge )
clIds    = which( isEipCl )
clIdList = lapply( clList, function(x) clIds[ x ] )

# Make the deduplicated gene list
plMbrList = lapply( clList, function(x) unlist( mbrList[x] ) )
names( plMbrList ) = sapply( clList, function(x) paste0( x, collapse = '/' ) )

# Make the cluster presence absence matrix
eipClMat = sapply( clList, function(j)
{
  if ( length( j ) == 1 ) return( eipClMat[ , j ] )
  return( as.integer( rowSums( eipClMat[ , j ] ) > 0 ) )
})
row.names( eipClMat ) = row.names( clMat )[ isEip ]
colnames( eipClMat )  = names( plMbrList )

# Make the cluster presence absence matrix
eipAcqEdge = sapply( clList, function(j)
{
  if ( length( j ) == 1 ) return( eipAcqEdge[ , j ] )

  clAcqEdges = sapply( 1:nrow( eipAcqEdge ), function(i)
  {
    egVal = unique( eipAcqEdge[ i, j ] )
    isEg  = egVal != 0
    if ( !TRUE %in% isEg ) return( 0 )
    egVal = egVal[ isEg ]
    if ( length( egVal ) > 1 )
      stop( "i = ", i, "; j = ", j, ": ", egVal, '\n' )
    return( egVal )
  })

  return( clAcqEdges )
})
row.names( eipAcqEdge ) = row.names( eipClMat )
colnames( eipAcqEdge )  = names( plMbrList )

# Make a presence absence matrix of each cluster across all of the genomes
plMat = sapply( clList, function(j)
{
  if ( length( j ) == 1 ) return( clMat[ , j ] )
  return( as.integer( rowSums( clMat[ , j ] ) > 0 ) )
})
row.names( plMat ) = row.names( clMat )
colnames( plMat )  = names( plMbrList )

gPlCounts = rowSums( eipClMat )

annots = annots[ row.names( annots ) %in% row.names( plMat ), ]

identical(  row.names( annots ), row.names( eipClMat ) )
countSmry =  summary( gPlCounts )
countSmry = round( countSmry, 2 )
countSmry = sapply(1:length( countSmry ),
  function(i) paste0( "  -- ", names( countSmry[i] ), ": ", countSmry[i], '\n' )
  )

isPl  = gPlCounts > 0
isKpc = annots$KPC == "KPC+"

nNoCl = sum( !isPl & isKpc )
nNoCl = sum( !isPl & !isKpc )

xx= matrix( sample( 4, 16, replace = TRUE ), 4, 4 )
hmCols = colorRampPalette( c( "#395983", "#e7e7e7", "#952827" ) )( 4 )
pheatmap::pheatmap( xx, color = hmCols )

nClKpcNeg = sum( isPl & !isKpc )
nClKpcNeg = sum( gPlCounts != 0 & annots$KPC == "KPC-" )
cat(
  "Summary of plasmid counts in genomes:\n", countSmry,
  "  -- Genomes not assigned a clulster:", nNoCl, '\n',
  "  --", round( ( 1 -  nNoCl / nrow( eipClMat ) ) * 100 , 2 ),"% assigened\n"
  )

# ---- Save the data -----------------------------------------------------------

print( pryr::mem_used() )
rData = "../data/2022_02_16_analyze_hgt_cluster_distribution.rData"
save(
  file = rData,
  list = c( "eipClMat", "eipAcqEdge", "plMbrList", "plMat" )
  )

# ---- Make the figures --------------------------------------------------------

figDir = "../figures/2022_02_16_analyze_hgt_cluster_distribution/"

hCl = agnes( t( eipClMat ) )

# Read in the cognac fastTree
cgTree = midpoint.root( read.tree( "../data/cgTree.tre" ) )

pheatmap(
  plMat[ GetTreeOrder( cgTree, row.names( plMat ) ), hCl$order ],
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  color             = c( "#e7e7e7", "#395983" ),
  annotation_row    = annots,
  annotation_colors = annotCols,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  filename          = paste0( figDir, "2022_02_16_merged_cl_mat.png" )
  )

treeOrder = GetTreeOrder( cgTree, row.names( eipClMat ) )

hmCols = colorRampPalette( c( "#395983", "#e7e7e7", "#952827" ) )( 100 )

pheatmap(
  eipClMat[ treeOrder, hCl$order ],
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  color             = c( "#e7e7e7", "#395983" ),
  annotation_row    = annots[ , 1:3 ],
  annotation_colors = annotCols,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  filename          = paste0( figDir, "2022_02_16_eip_merged_cl_mat.png" )
  )

# Make bar-plots showing the cluster distribution in all genomes
factorMats = lapply( seq( 3 ),
  function(j) MakeFeaturMat( plMat, names( annotCols[[j]] ), j )
  )
plOrder = order( colSums( eipClMat ), decreasing = TRUE )


png(
  filename = paste0( figDir, "2022_02_16_cluster_counts_by_feature.png" ),
  units    = "in",
  width    = 8,
  height   = 5,
  res      = 300
  )
par( mfrow = c( 1, 3 ) )
for ( i in 1:length( factorMats ) )
{
  barplot(
    factorMats[[i]][ , plOrder ],
    col    = annotCols[[i]],
    border = annotCols[[i]],
    main   = colnames( annots )[i],
    xlab   = "Gene Cluster",
    ylab   = "Frequency"
    )
  legend(
    "topright",
    legend = names( annotCols[[i]] ),
    col    = annotCols[[i]],
    pch    = 15
    )
}
dev.off()

# Make bar-plots showing the cluster distribution in the EIP genomes
factorMats = lapply( seq( 3 ),
  function(j) MakeFeaturMat( eipClMat, names( annotCols[[j]] ), j )
  )

png(
  filename = paste0( figDir, "2022_02_16_eip_cluster_counts_by_feature.png" ),
  units    = "in",
  width    = 8,
  height   = 5,
  res      = 300
  )
par( mfrow = c( 1, 3 ) )
for ( i in 1:length( factorMats ) )
{
  barplot(
    factorMats[[i]][ , plOrder ],
    col    = annotCols[[i]],
    border = annotCols[[i]],
    main   = colnames( annots )[i]
    )
  legend(
    "topright",
    legend = names( annotCols[[i]] ),
    col    = annotCols[[i]],
    pch    = 15
    )
}
dev.off()


par( mfrow = c( 3, 1 ) )

eipSource = c( "CT", "MN", "TN" )
isSource = lapply( eipSource,
  function(x)
    row.names( eipClMat ) %in% row.names( annots )[ annots$Source == x ]
  )

fileNames = paste0(
    figDir, "2022_02_16_", colnames( annots ), "_cluster_counts_by_state.png"
    )
for ( j in seq( 2 ) )
{
  png(
    filename = fileNames[j],
    units    = "in",
    width    = 5,
    height   = 7.5,
    res      = 300,
    )
  par( mfrow = c( 3, 1 ) )
  for ( i in seq( 3 ) )
  {
    mat =
      MakeFeaturMat( eipClMat[ isSource[[i]], ], names( annotCols[[j]] ), j )

    mText = paste( eipSource[i], colnames( annots )[j] )
    barplot(
      mat[ , plOrder ],
      col    = annotCols[[j]],
      border = annotCols[[j]],
      main   = mText,
      xlab   = "Gene Cluster",
      ylab   = "Frequency",
      cex.lab= 1.5,
      cex.axis = 1.5,
      cex.main = 1.5
      )
    # legend(
    #   "topright",
    #   legend = names( annotCols[[j]] ),
    #   col    = annotCols[[j]],
    #   pch    = 15
    #   )
  }
  dev.off()
}


nMrged = sapply( clIdList, length )
png(
  filename = paste0( figDir, "2022_02_16_n_merged_clusters.png" ),
  units    = "in",
  width    = 5,
  height   = 5,
  res      = 300
  )
barplot( nMrged, col = "#395983", border = "#395983",  )
dev.off()

# ------------------------------------------------------------------------------
