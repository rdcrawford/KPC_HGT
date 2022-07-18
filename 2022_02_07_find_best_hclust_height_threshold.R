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
library( scales )
library( cognac )
library( khroma )
library( future )
library( cluster )
library( pheatmap )
library( phytools )
library( future.apply )

# ---- Import the data ---------------------------------------------------------

# Load the mpr results
load( "../data/2022_02_07_run_mpr_kpc_plasmid_genes.rData" )

# Load the data on the KPC plasmids
load( "../data/2022_02_06_find_kpc_contig_genes.rData" )

# ---- Function definitions ----------------------------------------------------

GetMbrLists = function( clstrData )
{
  mbrListList = lapply( 1:ncol( clstrData ) , function(j)
  {
    clIds = unique( clstrData[ , j ] )
    clIds = clIds[ order( clIds ) ]
    mbrList = lapply( clIds,
      function(x) row.names( clstrData )[ clstrData[ , j ] == x ]
      )
    nGs = sapply( mbrList, length )
    print( j )
    print( summary( nGs ) )
    print( nGs[ order( nGs, decreasing = TRUE ) ][1:5] )
    return( mbrList[ 10 <= nGs & nGs < 150 ] )
  })
  return( mbrListList )
}

CalConfMats = function( mbrListList, contigGeneList )
{
  confMatList = lapply( mbrListList, function( mbrList )
  {
    confMats = lapply( contigGeneList, function(y)
    {

      gCounts = sapply( mbrList, function(x) sum( x %in% y ) )
      idx     = which.max( gCounts )
      return( CalcConfMat( mbrList[[ idx ]], y, nTotal ) )
    })
  })
  return( confMatList )
}

# Given a set of predicted positives and a set of true positives
# generate a confusion matrix
CalcConfMat = function( predPos, truePos, nTotal )
{
  tp = sum( predPos %in% truePos )
  fp = length( predPos ) - tp
  fn = length( truePos ) - tp
  tn = nTotal - ( tp + fp + fn )
  return( matrix( c( tp, fp, fn, tn ), nrow = 2, byrow = TRUE ) )
}

# Calculate the true positive rate from the input confusion matrix (sensitivity)
CalcTpr = function( confMat )
{
  nP = sum( confMat[ , 1 ] )
  if ( nP == 0 ) return( 0 )
  return( confMat[ 1, 1 ] / nP )
}

# Calculate the positive predicted value from the input confusion
# matrix (precision)
CalcPpv = function( confMat )
{
  return( confMat[ 1 , 1 ] / sum( confMat[ 1 , ] ) )
}

FindNumClsPerPl = function( mbrListList )
{
  nClsPerPlList = lapply( mbrListList, function( mbrList )
  {
    clMat = sapply( mbrList, function(x)
    {
      nGs = length(x)
      isCl = colnames( kpcContigMat ) %in% x
      nClGenes = sapply( 1:nrow( kpcContigMat ),
        function(i) sum( kpcContigMat[ i, isCl ] )
        )
      return( as.integer( nClGenes >= ceiling( nGs * 0.51 ) ) )
    })
    return( rowSums( clMat ) )
  })
  return( nClsPerPlList )
}

# ---- Calculate the stats on each height threshold value ----------------------

# Make a vector with the gene Ids of the kpc plasmid genes
genes  = as.integer( row.names( clstrData ) )
nTotal = nrow( clstrData )

# Subset the gene list to the genes for which AR was preformed
kpcPlGeneList = lapply( kpcPlGeneList, function(x) x[ x %in% genes ] )

# Make a list of the gene clusters at each height threshold
mbrListList = GetMbrLists( clstrData )

# Make a list of the confusion matrices for each plasmid
confMatList = CalConfMats( mbrListList, kpcPlGeneList )

# Calculate the statistics for each plasmid
tprList = lapply( confMatList, function(x) sapply( x, CalcTpr ) )
ppvList = lapply( confMatList, function(x) sapply( x, CalcPpv ) )
meanTpr = sapply( tprList, mean )
meanPpv = sapply( ppvList, mean )
nMbrs   = sapply( mbrListList, function(x) sapply( x, length ) )
nCls    = sapply( nMbrs, function(x) sum( x >= 10 ) )


plClCountList = FindNumClsPerPl( mbrListList )
nPlsWithCl    = sapply( plClCountList, function(x) sum( x != 0 ) )

nGeneInCl = sapply( mbrListList, function(x) length( unlist(x) ) )

clStats = cbind( meanPpv, meanTpr, nCls, nPlsWithCl )
print( clStats )

rankVals = sapply( 1:ncol( clStats ), function(j) rescale( clStats[ , j ] ) )

nn = sapply( 1:nrow( rankVals ), function(i) sum( rankVals[i, ] ) )
bestIdx = order( nn, decreasing = TRUE )[1]
print( bestIdx )
mbrList = mbrListList[[ bestIdx ]]

# Set names to keep track of the clusters that were removed
names( mbrList ) = seq( length( mbrList ) )

# Subset to the relevant clusters
nClMbrs = sapply( mbrList, length )
isGood  = which( nClMbrs >= 10 & nClMbrs != max( nClMbrs ) )
mbrList = mbrList[ isGood ]

# ---- Save the data -----------------------------------------------------------

print( pryr::mem_used() )
rData = "../data/2022_02_07_find_best_hclust_height_threshold.rData"
save( file = rData, list = "mbrList" )

# ---- Make the figures --------------------------------------------------------

figDir = "../figures/2022_02_07_find_best_hclust_height_threshold/"



plot(
  x   = meanTpr,
  y   = meanPpv,
  col = "#2E4769",
  pch = 16
  )
plot(
  x   = nCls,
  y   = meanPpv,
  col = "#2E4769",
  pch = 16
  )

fracInCl = sapply( kpcPlGeneList, function(x)
{
  nOnContig = sapply( trGeneList, function(y) sum( y %in% x ) )
  return( max( nOnContig ) / length(x) )
})
names( fracInCl ) = NULL

barplot(
  fracInCl[ order( fracInCl ) ],
  col    = "#2E4769",
  border = "#2E4769",
  xlab   = "Plasmid",
  ylab   = "Fraction of plasmid genes on transition edge"
  )

# distMat
# corrMat

clMbrId = vector( "integer", ncol( distMat ) )
for ( i in 1:length( mbrList ) )
  for ( j in which( colnames( distMat ) %in% mbrList[[i]] ) ) clMbrId[j] = i
clstrData = data.frame( clMbrId )
row.names( clstrData ) = colnames( distMat )
colnames( clstrData )  = "cluster"

# Get the annotations for each color
vals = unique( clstrData[ , 1 ] )
vals = vals[ order( vals ) ]
cols = c("white", color( "smooth rainbow" )( length( vals ) - 1 ) )
names( cols ) = vals
clCols = list( "cluster" = cols )


print( clMbrId[ 1:100 ] )
print( head( clstrData ) )
print( clCols )

figPath = paste0( figDir, "2022_02_07_dist_mat.png" )

# Make a vector of heatmap colors
hmCols = MakeHmCols( distMat )

# Write the heatmap to the supplied path
pheatmap(
  distMat[ hCl$order, hCl$order ],
  color             = hmCols,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  annotation_col    = clstrData,
  annotation_row    = clstrData,
  annotation_colors = clCols,
  annotation_legend = TRUE,
  file              = figPath
  )

# Make a vector of heatmap colors
hmCols = MakeHmCols( corrMat )

figPath = paste0( figDir, "2022_02_07_corr_mat.png" )

# Write the heatmap to the supplied path
pheatmap(
  corrMat[ hCl$order, hCl$order ],
  color             = hmCols,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  annotation_col    = clstrData,
  annotation_row    = clstrData,
  annotation_colors = clCols,
  annotation_legend = FALSE,
  file              = figPath
  )

# ------------------------------------------------------------------------------