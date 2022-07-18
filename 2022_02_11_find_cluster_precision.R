# ------------------------------------------------------------------------------
#
# 2022/02/11
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

# Column indices of the gff files
G_ID       = 1
G_ANNOT    = 2
CONTIG     = 3
G_START    = 4
G_END      = 5
G_CLUST    = 7
G_NT_CLUST = 8

# Column indices of the cluster stat data
G_ID       = 1
C_ID       = 2
CL_CONTIG  = 3
N_CL_GENES = 4
N_TP       = 5
N_FP       = 6
PPV        = 7

# ---- Import the data ---------------------------------------------------------

# Load the parsed gff files with the nt clusters
load( "../data/2022_02_04_parse_plasmid_gene_algns.rData" )

# Load the data on kpc contigs
load( "../data/2022_02_06_find_kpc_contig_genes.rData" )

# Load the data on the gene clusters
load( "../data/2022_02_08_get_cl_gene_counts_by_edge.rData" )

#
load( "../data/2022_02_16_analyze_hgt_cluster_distribution.rData" )

load( "../data/eipData.rData" )

isolateMetaData = read.table(
  "../data/isolateMetaData.tsv",
  sep              = '\t',
  stringsAsFactors = FALSE,
  header           = TRUE
  )

# ---- Function definitions ----------------------------------------------------

GetChrContigs = function( geneEnv, idxs )
{
  chrContigId = sapply( geneEnv$gfList[ idxs ], function(gff)
  {
    cVals   = unique( gff[ , CONTIG ] )
    nContig = sapply( cVals, function(x) sum( gff[ , CONTIG ] == x ) )
    return( cVals[ which.max( nContig ) ] )
  })
  return( chrContigId )
}

# Find the contig with the most cluster genes
GetClContig = function( gff, isClGene )
{
  contigVals = unique( gff[ isClGene, CONTIG ] )
  if ( length( contigVals ) == 1 ) return( contigVals )
  nContig = sapply( contigVals,
    function(x) sum( gff[ isClGene, CONTIG ] == x )
    )
  return( contigVals[ which.max( nContig ) ] )
}

# Get the number of genes on the best match contig
GetNumTp = function( gff, isClGene, clContig )
{
  isClContig = gff[ isClGene, CONTIG ] == clContig
  return( sum( isClContig ) )
}

# Get the statistics for the input cluster genes
GetClMbrStats = function( gff, clMbrs )
{
  # Find the entries corresponding to this cluster
  isClGene = gff[ , G_NT_CLUST ] %in% clMbrs
  if ( sum( isClGene ) == 0 )
  {
    # save( file = "../data/gff.rData", list = c( "gff", "clMbrs" ) )
    warning( "No cluster gene identified..." )
  }

  # Find the contig that this cluster is on
  clContig = GetClContig( gff, isClGene )

  # Get the number of genes on the best contig and calculate precision
  nTp = GetNumTp( gff, isClGene, clContig )
  nCl = sum( isClGene )
  nFp = nCl - nTp
  ppv = nTp / nCl

  return( c( clContig, nCl, nTp, nFp, ppv ) )
}

# clIdxs = which( clMat[ i, ] == 1 )
# gClstrList = lapply( clIdxs,
#   function(j) c( i, j, GetClMbrStats( gff, mbrList[[j]] ) )
#   )
# The the data on the number of genes which are correctly places
MakeClStatDf = function( geneEnv, clMat, isCompKpc )
{
  # For each genome, find the clusters that this genome is predicted to have
  # each cluster
  clStatList = lapply( isCompKpc, function(i)
  {
    clIdxs = which( clMat[ i, ] == 1 )
    gClstrList = lapply( clIdxs,
      function(j) c( i, j, GetClMbrStats( geneEnv$gfList[[i]], mbrList[[j]] ) )
      )
    return( gClstrList )
  })
  statMat =
    sapply( unlist( clStatList, recursive = FALSE ), function(x) return(x) )


  clStats = data.frame(
    "genomeNames" = geneEnv$genomeNames[ statMat[ G_ID, ] ],
    "clId"        = as.integer( statMat[ C_ID, ] ),
    "clContig"    = as.integer( statMat[ CL_CONTIG, ] ),
    "nClGenes"    = as.integer( statMat[ N_CL_GENES, ] ),
    "nTp"         = as.integer( statMat[ N_TP, ] ),
    "nFp"         = as.integer( statMat[ N_FP, ] ),
    "ppv"         = as.numeric( statMat[ PPV, ] )
    )
  return( clStats )
}

GetGeneContigs = function( gff, geneId )
{
  kpcIdxs = which( gff[ , G_NT_CLUST ] == geneId )
  return( unique( gff[ kpcIdxs, CONTIG ] ) )
}

# ---- Find the precision using complete genomes -------------------------------

clMat = plMat
mbrList = plMbrList

# Sort the clusters to the same order as the meta-data
mdOrd = sapply( row.names( clMat ),
  function(x) which( isolateMetaData$genomeNames == x )
  )
clMat = clMat[ order( mdOrd ), ]

# Make sure everything is sorted
if ( !identical( row.names( clMat ), isolateMetaData$genomeNames ) )
  stop( "row.names( clMat ) != isolateMetaData$genomeNames" )
if ( !identical( geneEnv$genomeNames, isolateMetaData$genomeNames ) )
  stop( "geneEnv$genomeNames != isolateMetaData$genomeNames" )

# Find the index of each complete KPC+ genome
isCompKpc =
  which( isolateMetaData$kpc == 1 & isolateMetaData$genomeStatus == "complete" )

# Find the contig corresponding to the chromosome
chrContigs = GetChrContigs( geneEnv, isCompKpc )

# Find the KPC contigs
kpcContigs = lapply( isCompKpc,
  function(i) GetGeneContigs( geneEnv$gfList[[i]], kpc )
  )

# Make a data-frame to store the cluster statistics
clStats = MakeClStatDf( geneEnv, clMat, isCompKpc )

# ----

kpcPlContigs = sapply( 1:length( chrContigs ),
  function(i) kpcContigs[[i]][ kpcContigs[[i]] != chrContigs[i] ]
  )

gIdList = lapply( 1:length( isCompKpc ), function(i)
{
  gId = isolateMetaData$genomeNames[ isCompKpc[i] ]
  return( rep( gId, length( kpcPlContigs[[i]] ) ) )
})

plGIds    = unlist( gIdList )
plContIds = unlist( kpcPlContigs )
clVals    = seq( length( mbrList ) )
nPls      = length( plGIds )

kpcPlClstrMat = t( sapply( seq( nPls ), function(i)
{
  isKpcG =
      clStats[ , G_ID ] == plGIds[i] & clStats[ , CL_CONTIG ] == plContIds[i]
  return( as.integer( clVals %in%  clStats[ isKpcG, C_ID ] ) )
}))
clCounts = rowSums( kpcPlClstrMat )
table( clCounts )

noClGIds = plGIds[ clCounts == 0 ]

nGenesInCl = sapply( seq( nPls ), function(i)
{
  clIdxs = which( kpcPlClstrMat[ i, ] == 1 )
  return( sum( kpcPlGeneList[[i]] %in% unlist( mbrList[ clIdxs ] ) ) )
})

nPlGenes = sapply( kpcPlGeneList, length )
plFrac   = sapply( seq( nPls ), function(i) nGenesInCl[i] / nPlGenes[i] )

fracClMbrMat = sapply( mbrList,
  function(x)
    rowSums( kpcContigMat[ , colnames( kpcContigMat ) %in% x ] ) / length(x)
  )

fracGenes = lapply( which( clCounts == 0 ),
  function(i) fracClMbrMat[ i, fracClMbrMat[ i, ] != 0 ]
  )
fnCls = lapply( which( clCounts == 0 ),
  function(i) which( fracClMbrMat[ i, ] >= 0.51 )
  )

# ---- Save the data -----------------------------------------------------------

print( pryr::mem_used() )
rData = "../data/2022_02_11_find_cluster_precision.rData"
save(
  file = rData,
  list = c( "clStats", "chrContigs", "kpcContigs", "kpcPlClstrMat" )
  )

# ---- Make the figures --------------------------------------------------------

figDir = "../figures/2022_02_11_find_cluster_precision/"

posCounts = c( sum( clStats[ , N_TP ] ), sum( clStats[ , N_FP ] ) )
names( posCounts ) = c( "True positives", "False positives" )

cat( "Total precision:", posCounts[1] / sum( posCounts ), '\n' )

png(
  filename = paste0( figDir, "2022_02_11_tp_vs_fp_barplot.png" ),
  units    = "in",
  width    = 4,
  height   = 7,
  res      = 300
  )
par( mar = c( 10.1, 4.1, 4.1, 2.1 ) )
barplot(
  posCounts,
  las    = 2,
  col    = "#395983",
  border = "#395983",
  ylim   = c( 0, 40000 )
  )
dev.off()


plFrac = plFrac[ order( plFrac ) ]
names( plFrac ) = NULL
png(
  filename = paste0( figDir, "2022_02_11_tp_vs_fp_barplot.png" ),
  units    = "in",
  width    = 4,
  height   = 7,
  res      = 300
  )
par( mar = c( 3, 4.1, 4.1, 2.1 ) )
barplot( plFrac, las = 2, col = "#395983", border = "#395983", ylim = c( 0, 1) )
dev.off()


pheatmap(
  kpcPlClstrMat,
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  color             = c( "#e7e7e7", "#395983" ),
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  filename          = paste0( figDir, "2022_02_11_kpc_pl_clstr_mat.png" )
  )


clIds = unique( clStats$clId )

isCl = lapply( clIds, function(x) clStats$clId == x )



ppvByCl = sapply( isCl, function(x)
{
  posCounts = c( sum( clStats[ x, N_TP ] ), sum( clStats[ x, N_FP ] ) )
  return( posCounts[1] / sum( posCounts ) )
})

nClGenomes = sapply( isCl, sum )

isEip = row.names( clMat ) %in% eipData$Sample_ID

nFp = sapply( isCl, function(x) mean( clStats[ x, N_FP ] ) )


nEipGenomes = colSums( clMat[ isEip, clIds ] )
names( nEipGenomes ) = NULL

isGood = ppvByCl >= 0.95
sum( isGood / length( isGood ) )
sum( nEipGenomes[ !isGood ] )


png(
  filename = paste0( figDir, "2022_02_11_precision_by_cluster.png" ),
  units    = "in",
  width    = 8,
  height   = 8,
  res      = 300
  )
par( mar = c( 2, 4.1, 4.1, 2.1 ), mfrow = c( 2, 2 ) )
plOrd = order( ppvByCl )
cols = khroma::color( "smooth rainbow" )( 101 )
plCols = sapply( round( scales::rescale( ppvByCl[ plOrd ] ) * 100 + 1, 0 ),
  function(x) cols[x]
  )
barplot(
  ppvByCl[ plOrd ],
  col    = plCols,
  border = plCols,
  main   = "Precision"
  )
barplot(
  nFp[ plOrd ],
  col    = plCols,
  border = plCols,
  main   = "Average number of false positives"
  )
legend(
  "topright",
  legend = round( rev( seq( 0.65, 1, length.out = 11 ) ), 2 ),
  col    = rev( cols[ seq( 1, 101, 10 ) ] ),
  pch    = 15
  )
barplot(
  nClGenomes[ plOrd ],
  col    = plCols,
  border = plCols,
  main   = "Number of complete genomes"
  )
barplot(
  nEipGenomes[ plOrd ],
  col    = plCols,
  border = plCols,
  main   = "Number of EIP genomes"
  )
dev.off()



gnPpv = sapply( 1:nrow( clStats ),
  function(i) clStats[ i, N_TP ] / sum( clStats[ i, c( N_TP, N_FP ) ] )
  )
summary( gnPpv )

clGnPpv = lapply( isCl, function(x) gnPpv[x] )

isPerfect = sapply( clGnPpv, function(x) sum( x != 1 ) == 0 )
sum( isPerfect ) / length( isPerfect )

clGnPpvList = lapply( isCl, function(x) gnPpv[x] )
ppvSummary = sapply( clGnPpvList, summary )

summary( ppvSummary[ 3, ] )

sum( gnPpv == 1 ) / length( gnPpv )
sum( gnPpv >= 0.95 ) / length( gnPpv )

# How many clusters are observed on KPC+ plasmids, KPC- plasmids, and the
# chromosome
kpcCongtigCounts = sapply( 1:length( isCompKpc ), function(i)
{
  isG       = clStats[ , G_ID ] == isolateMetaData$genomeNames[ isCompKpc[i] ]
  isKpcCont = clStats[ isG, CL_CONTIG ] %in% kpcContigs[[i]]
  isChr     = clStats[ isG, CL_CONTIG ] == chrContigs[i]
  return( c( sum( isKpcCont ), sum( !isKpcCont & !isChr ), sum( isChr ) ) )
})
nOnKpcContig = c(
  "Cluster on KPC+ Plasmid" = sum( kpcCongtigCounts[ 1, ] ),
  "Cluster on KPC- Plasmid" = sum( kpcCongtigCounts[ 2, ] ),
  "Cluster on chromosome"   = sum( kpcCongtigCounts[ 3, ] )
  )
par( mar = c( 12, 4.1, 4.1, 2.1 ) )
png(
  filename = paste0( figDir, "2022_02_11_number_of_clusters_on_KPC_pl.png" ),
  units    = "in",
  width    = 4,
  height   = 7,
  res      = 300
  )
barplot( nOnKpcContig, las = 2, col = "#395983", border = "#395983" )
dev.off()

# How many clusters are observed on KPC+ plasmids, KPC- plasmids, and the
# chromosome
kpcCongtigCounts = sapply( 1:length( isCompKpc ), function(i)
{
  isG       = clStats[ , G_ID ] == isolateMetaData$genomeNames[ isCompKpc[i] ]
  isKpcCont = clStats[ isG, CL_CONTIG ] %in% kpcContigs[[i]]
  isChr     = clStats[ isG, CL_CONTIG ] == chrContigs[i]
  return( c( sum( isKpcCont ), sum( !isKpcCont & !isChr ), sum( isChr ) ) )
})
nOnKpcContig = c(
  "Cluster on KPC+ Plasmid" = sum( kpcCongtigCounts[ 1, ] ),
  "Cluster on KPC- Plasmid" = sum( kpcCongtigCounts[ 2, ] ),
  "Cluster on chromosome"   = sum( kpcCongtigCounts[ 3, ] )
  )
png(
  filename = paste0( figDir, "2022_02_11_sig_clusters_on_KPC_pl.png" ),
  units    = "in",
  width    = 4,
  height   = 7,
  res      = 300
  )
par( mar = c( 12, 4.1, 4.1, 2.1 ) )
barplot( nOnKpcContig, las = 2, col = "#395983", border = "#395983" )
dev.off()

# ------------------------------------------------------------------------------