# ------------------------------------------------------------------------------


pheatmap(
  ,
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  color             = c( "#e7e7e7", "#395983" ),#hmCols,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  filename          = NA#paste0( figDir, ".png" )
  )

xx = plGeneMat[ , mbrList[[4]] ]
dim( xx )
colSums( xx )
# ------------------------------------------------------------------------------

xx = cophenetic.phylo( plTrees[[i]] )
hmCols = colorRampPalette( c( "#395983", "#e7e7e7", "#952827" ) )( 100 )
pheatmap(
  xx,
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  color             = hmCols,
  annotation_row    = annots[ , c( 1, 2, 3 )],
  annotation_colors = annotCols[ c( 1, 2, 3 ) ],
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  filename          = NA#paste0( figDir, ".png" )
  )

# ------------------------------------------------------------------------------

pheatmap(
  plMbrMats[[1]][ kpcPlMat[ , 1 ], ],
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  color             = c( "#e7e7e7", "#395983" ),#hmCols,
  annotation_row    = annots,
  annotation_colors = annotCols,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  filename          = NA#paste0( figDir, ".png" )
  )
pheatmap(
  plMbrMats[[2]][ kpcPlMat[ , 2 ], ],
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  color             = c( "#e7e7e7", "#395983" ),#hmCols,
  annotation_row    = annots,
  annotation_colors = annotCols,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  filename          = NA#paste0( figDir, ".png" )
  )


clGenes = unlist( mbrList )
isCl = colnames( trMat ) %in% clGenes
clTrMat = trMat[ , isCl ]

nTrGs = ncol( clTrMat )

isGeneEdge = lapply( seq( nTrGs ), function(j) which( clTrMat[ , j ] != 0 ) )
shTrMat = matrix( 0, nTrGs, nTrGs )
for ( i in seq( nTrGs - 1 ) )
{
  shTrMat[ i, i ] = length( isGeneEdge[[i]] )
  for ( j in seq( i + 1, nTrGs ) )
  {
    nn = sum( isGeneEdge[[i]] %in% isGeneEdge[[j]] )
    shTrMat[ i, j ] = nn
    shTrMat[ j, i ] = nn
  }
}


clGeneIdxs = lapply( mbrList, function(x) which( colnames( clTrMat ) %in% x ) )

i = which( names( mbrList ) == "2025")
i = which( names( mbrList ) == "6")
mbrList[i]
sapply( clGeneIdxs, function(i)
{
  clstrMat = shTrMat[ clGeneIdxs[[i]], clGeneIdxs[[i]] ]
  as.dist( clstrMat )

sapply( 1:nrow( clstrMat ), function(i) sum( clstrMat[ i, ] > 1  ) )


})

# ------------------------------------------------------------------------------

ii = unique( unlist( isGeneEdge[ clGeneIdxs[[i]] ] ) )
xx = clTrMat[ ii, clGeneIdxs[[i]] ]
pheatmap(
  xx,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  color             = c( "#e7e7e7", "#395983", "#952827" ),
  # annotation_row    = annots,
  # annotation_colors = annotCols,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  filename          = NA#paste0( figDir, ".png" )
  )

# ------------------------------------------------------------------------------

swathPcrData$Sample.Number
"801a"

ii = row.names( swathPcrData ) %in% isolateMetaData[ ptIdxList[[3]], 1 ]
swathPcrData[ ii, ]
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------