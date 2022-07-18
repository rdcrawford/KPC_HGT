
# ------------------------------------------------------------------------------
#
#
# Ryan D. Crawford
# ------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( cognac )
library( ape )
library( maps )
library( phytools )

# ---- Function definitions ----------------------------------------------------

source( "R/PlotTree.R" )

# ---- Variable initializations ------------------------------------------------

isolateMetaData = read.table(
  "../data/isolateMetaData.tsv",
  sep = '\t',
  stringsAsFactors = FALSE,
  header = TRUE
  )
# load( "../data/eipData.rData")

for ( i in which( is.na( isolateMetaData$kpc ) ) ) isolateMetaData$kpc[i] = 0

# ------------------------------------------------------------------------------

sts = c(
  "K. pneumoniae ST307",
  "K. pneumoniae CG258",
  "E. cloacae ST114",
  "E. cloacae ST171"
  )

mlst = sapply( isolateMetaData$mlst, function(x)
{
  if ( x %in% sts ) return( x)
  return( "Other" )
})

# ------------------------------------------------------------------------------

enterobacterSp = c(
  "Enterobacter asburiae", "Enterobacter kobei", "Enterobacter roggenkampii",
  "Enterobacter xiangfangensis"
  )
klebsiellaSp = c(
  "Klebsiella aerogenes", "Klebsiella michiganensis", "Klebsiella oxytoca",
  "Klebsiella variicola"
  )

species = sapply( isolateMetaData$species, function(x)
{
  if ( x == "Leclercia sp." || x == "Raoultella sp." ) return( "Other" )
  if ( x %in% enterobacterSp ) return( "Enterobacter sp." )
  if ( x %in% klebsiellaSp ) return( "Klebsiella sp." )
  return( x )
})
table( species )

# ------------------------------------------------------------------------------

sources = sapply( isolateMetaData$sources , function(x)
{
  x = gsub("EIP_", '', x )
  if ( x %in% c( "CT", "MN", "TN" ) ) return( x )
  return( "Other" )
})

# ------------------------------------------------------------------------------

annots = cbind.data.frame(
  "Species" = species,
  "ST"      = mlst,
  "Source"  = sources,
  "Status"  = isolateMetaData$genomeStatus,
  "KPC"     = ifelse( isolateMetaData$kpc == 1, "KPC+", "KPC-")
  )
row.names( annots ) = isolateMetaData$genomeNames

# ------------------------------------------------------------------------------

cols = c(
  "#ba3231", "#952827", "#77201f", # Red  "#A11B1B", "#d76a69" "#d04e4d",
  "#f37a48",                       # Orange
  "#f0bc5e", "#eaa221",            # Yellow
  "#A6DEAE", "#80B28D", "#5b876d", # Muted Green
  "#00eb5a", "#00bc48", "#00963a", # Bright Green
  "#717d6d", "#5a6457", "#485046", # Dark green brown
  "#7c81bd", "#545aa7", "#434886", # Purple Blue
  "#395983", "#1697B7", "#2E4769",            # Blue
  "#CFC1CE", "#725971", "#702670", # Purple
  "#e7e7e7", "#a1a1a1", "#767676", # Gray
  "#ffffff"                        # White
  )
barplot( rep( 10, length( cols ) ), col = cols, border = cols )


stCols = c(
  `K. pneumoniae ST307` = "#b0b3d7",
  `K. pneumoniae CG258` = "#545aa7",
  `E. cloacae ST114`    = "#A6DEAE",
  `E. cloacae ST171`    = "#80B28D",
  `Other`               = "#a1a1a1"
  )

statusCols = c(
  "wgs"       = "#a1a1a1",
  "complete"  = "#395983",
  "long-read" = "#ba3231"
  )

sourceCols = c(
  "CT"    = "#952827",
  "MN"    = "#2E4769",
  "TN"    = "#725971",
  "Other" = "#a1a1a1"
  )

spCols = c(
  "K. pneumoniae"    = "#2b2e56",
  "Klebsiella sp."   = "#7f4485",
  "E. cloacae"       = "#5b876d",
  #"Enterobacter sp." = "#485046",
  "Citrobacter sp."  = "#53e9ef", # "#CFC1CE",
  "E. coli"          = "#1697B7",
  "Other"            = "#a1a1a1"
  )

kpcCols = c(
  "KPC+" = "#00963a",
  "KPC-" = "#a1a1a1"
  )


annotCols = list( spCols, stCols, sourceCols, statusCols, kpcCols )
names( annotCols ) = colnames( annots )

ac = unlist( annotCols)
barplot( rep( 10, length( ac ) ), col = ac, border = ac )

# ---- Save the data -----------------------------------------------------------

rData = "../data/2022_02_07_make_annotations_and_colors.rData"
save( file = rData, list = c( "annotCols", "annots" ) )

# ---- Make the figures --------------------------------------------------------

figDir = "../figures/2022_02_07_make_annotations_and_colors/"

for ( j in 1:ncol( annots ) )
{
  cc = sapply( names( annotCols[[j]] ), function(x) sum( annots[ , j ] == x ) )
  png(
    filename = paste0( figDir, "2022_02_07_", colnames( annots )[j], ".png" ),
    units    = "in",
    width    = 5,
    height   = 5,
    res      = 300
    )
  par( mar = c( 10, 4.1, 4.1, 2.1 ) )
  barplot(
    cc,
    col    = annotCols[[j]],
    border = annotCols[[j]],
    las    = 2,
    ylab   = "Frequency"
    )
  dev.off()
}

# Read in the coganc fastTree
cognacDir  = "../analysis/2022_02_04_eip_cognac_analysis/"
cgTreePath = paste0( cognacDir, "concatenated_gene_aa_alignment_fasttree.tre" )
cgTree     = read.tree( cgTreePath )
cgTree     = midpoint.root( cgTree )

png(
  filename = paste0( figDir, "2022_02_07_core_gene_tree.png" ),
  units    = "in",
  width    = 5,
  height   = 7.5,
  res      = 300
  )
PlotTree( cgTree, annots, annotCols  )
dev.off()

# ------------------------------------------------------------------------------