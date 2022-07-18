# ------------------------------------------------------------------------------
# Select analysis isolates
# 2022/02/04
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# Pick the genomes to be included in the analysis to identify HGT of
# ------------------------------------------------------------------------------

# ---- Function definitions ----------------------------------------------------

GetNbrs = function( distMat )
{

  minDistIsos = lapply( 1:nrow( distMat ), function(i)
  {
    minVal = min( distMat[ i, ] )
    return( colnames( distMat )[ distMat[ i, ] == minVal ] )
  })

  return( minDistIsos )
}

GetInDistIsos = function( distMat, thrVal )
{

  withinDistIsos = lapply( 1:nrow( distMat ),
    function(i) colnames( distMat )[ distMat[ i, ] <= thrVal ]
    )

  return( withinDistIsos )
}

# ---- Variable initializations ------------------------------------------------

# Directory with output files from the cognac analysis
outDir = "../../genomes/analysis/2021_04_26_cognac_analysis/"

# Big matrix
load( paste0( outDir, "cognac_analysis_4_1.rData" ) )
load( "../../genomes/data/2021_03_26_download_patric_genomes.rData" )
load( "../data/eipData.rData" )


# Load metadata
isolateMetaData = read.table(
  "../../genomes/data/isolateMetaData.tsv",
  sep = '\t',
  row.names = 1
  )

cdcIsos = read.table(
  "../data/2022_02_03_PATRIC_CDC_isolate_data.tsv",
  sep = '\t',
  row.names = 1
  )

isEip   = row.names( distMat ) %in% eipData$Sample_ID
distMat = distMat[ isEip, !isEip ]

isDup = colnames( distMat ) %in% cdcIsos$Patric.Genome.ID[ cdcIsos$isDuplicate ]
distMat = distMat[ , !isDup ]

for ( i in which( is.na( isolateMetaData$kpc ) ) )
  isolateMetaData$kpc[i] = FALSE
isKpc       = isolateMetaData$kpc == 1
kpcIsos     = isolateMetaData$genomeNames[ isKpc ]
isComp      = isolateMetaData$genomeStatus == "complete"
compGenomes = isolateMetaData$genomeNames[ isKpc & isComp ]


compDists = t( distMat[ , colnames( distMat ) %in% compGenomes ] )
minDists  = sapply( 1:nrow( compDists ), function(i) min( compDists[ i, ] ) )

# ------------------------------------------------------------------------------

dist5Isos = unique( unlist( GetInDistIsos( distMat, 5  ) ) )
dist10Isos = unique( unlist( GetInDistIsos( distMat, 10 ) ) )
dist25Isos = unique( unlist( GetInDistIsos( distMat, 25 ) ) )
length( dist5Isos )
length( dist10Isos )
length( dist25Isos )

# ------------------------------------------------------------------------------


isMn   = patricData$BioProject.Accession == "PRJNA272863"
mnIsos = patricData$Genome.ID[ isMn ]

cdcGenomes     = c( cdcIsos$Patric.Genome.ID, mnIsos )
isCdc          = colnames( distMat ) %in% cdcGenomes
patDistMat     = distMat[ , !isCdc ]
isKpc          = colnames( patDistMat ) %in% kpcIsos
kpcNbrIsos     = unique( unlist( GetNbrs( patDistMat[ , isKpc ] ) ) )
kpcNegNbrIsos  = unique( unlist( GetNbrs( patDistMat[ , !isKpc ] ) ) )
nearestNbrIsos = unique( unlist( GetNbrs( patDistMat ) ) )

# ------------------------------------------------------------------------------

analysisIsos =
  unique( c( eipData$Sample_ID, kpcNbrIsos, kpcNegNbrIsos , compGenomes  ) )
"2059860.3" %in% analysisIsos

sum( analysisIsos %in% kpcIsos ) / length( analysisIsos )

colIdxs     = c( 1, 2, 3, 4, 5, 6, 8, 15, 16, 17, 31, 32, 33 )
isAnalysis  = isolateMetaData$genomeNames %in% analysisIsos
eipMetaData = isolateMetaData[ isAnalysis, colIdxs ]

write.table(
  eipMetaData,
  "../data/isolateMetaData.tsv",
  sep = '\t'
  )

# ------------------------------------------------------------------------------