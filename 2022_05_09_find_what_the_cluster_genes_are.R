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

# Column indices of the gff files
G_ID       = 1
G_ANNOT    = 2
CONTIG     = 3
G_START    = 4
G_END      = 5
G_CLUST    = 7
G_NT_CLUST = 8

# ---- Import the data ---------------------------------------------------------

isolateMetaData = read.table(
  "../data/isolateMetaData.tsv",
  sep = '\t',
  stringsAsFactors = FALSE,
  header = TRUE
  )

# Load the annotations and colors: "annots" and "annotCols"
load( "../data/2022_02_07_make_annotations_and_colors.rData" )

# Includes the data on the high identity genes: "clustList", "plGeneMat",
# "geneEnv", "algnPaths", and "ntToAaClId"
load( "../data/2022_02_04_parse_plasmid_gene_algns.rData" )

# Load the data on the gene clusters
# load( "../data/2022_02_08_get_cl_gene_counts_by_edge.rData" )
load( "../data/2022_04_14_compare_plasmid_alignment_distances.rData" )

# ---- Function definitions ----------------------------------------------------

GetClGeneAnnots = function( geneId )
{
  isFound = FALSE
  curIdx  = 1
  for ( i in 1:length( geneEnv$gfList ) )
  {
    idx = which( geneEnv$gfList[[i]][ , G_NT_CLUST ] == geneId )
    if ( length( idx ) != 0 ) return( geneEnv$gfList[[i]][ idx[1], G_ANNOT ] )
  }
  stop( "geneId: ", geneId, " not found..." )
}

# ----

# Subset to the
plMbrList = plMbrList[ isSigCl ]


#
geneAnnotList = lapply( 1:length( plMbrList ),
  function(j) sapply( plMbrList[[j]], GetClGeneAnnots )
  )

# ---- Save the data -----------------------------------------------------------

print( pryr::mem_used() )
rData = "../data/2022_05_09_find_what_the_cluster_genes_are.rData"
save( file = rData, list = "geneAnnotList" )

# ---- Make the figures --------------------------------------------------------

# load( "../data/2022_05_09_find_what_the_cluster_genes_are.rData" )
annots      = unlist( geneAnnotList )
annotCounst = table( annots )
annotCounst = annotCounst[ order( annotCounst, decreasing = TRUE ) ]
geneWdList  = lapply( geneAnnotList, function(x) strsplit( x, ' ') )
clWords     = unlist( geneWdList )
clWords     = clWords[ !is.na( clWords ) ]
wordCounts  = table( clWords )
wordCounts  = wordCounts[ order( wordCounts, decreasing = TRUE ) ]

par( mar = c(5.1, 10, 4.1, 2.1 ), las = 2)
barplot(
  wordCounts[1:50],
  las = 2,
  col = "#395983",
  border = "#395983",
  horiz  = TRUE
  )
par( mar = c( 5.1, 30, 4.1, 2.1 ), las = 2 )

barplot(
  annotCounst[2:51],
  las    = 2,
  col    = "#395983",
  border = "#395983",
  horiz  = TRUE
  )

toCount = list(
  "Hypothetical" = c( "hypothetical protein" ),
  "Heavy metal\nResistance" = c( "resistance" ),
  "DNA processing" = c(
    "partitioning", "polymerase", "replication", "mobile element",
    "transposase", "dna", "methyltransferase", "helicase", "replication",
    "resolvase", "endonuclease"
    ),
  "Transcriptional\nregulator" = c( "transcriptional regulator" ),
  "Conjugation"    = c( " tra","type iv", "pilin", "conjug", "mobilization"),
  "Toxin"          = c( "colicin", "toxin", "type vi" ),
  "Phage"         = c( "phage")
 )

annots = tolower( annots )

keyWordIdxs = lapply( toCount, function(x)
{
  keyWordIdxs = lapply( x, function(y) grep( y, annots, ignore.case = TRUE ) )
  return( unique( unlist( keyWordIdxs ) ) )
})

annotIdxs     = unlist( keyWordIdxs )
isNotIn       = !( seq( length( annots ) ) %in% annotIdxs )
keyWordCounts = sapply( keyWordIdxs, length )
keyWordCounts = keyWordCounts[ order( keyWordCounts ) ]
keyWordCounts = c( "Other" = sum( isNotIn ), keyWordCounts )
keyWordCounts / length( annots )

png(
  filename = "../figures/2022_05_09_find_what_the_cluster_genes_are.png",
  units    = "in",
  width    = 8,
  height   = 6.5,
  res      = 300
  )
par( mar = c( 5.1, 7.5, 4.1, 2.1 ) )
barplot(
  keyWordCounts,
  las    = 1,
  col    = "#395983",
  border = "#395983",
  horiz  = TRUE
  )
dev.off()

terIdxs = sapply( geneAnnotList, function(x) grep( "Ter", x) )
which( sapply( terIdxs, length ) != 0 )

# ------------------------------------------------------------------------------