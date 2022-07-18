# ------------------------------------------------------------------------------
# Look up MLST
# 2022/02/07
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# Use the tree to find the MLST and species of each genome in the analysis.
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( ape )
library( maps )
library( phytools )

# ---- Import the data ---------------------------------------------------------

# Load the metadata from the patric and EIP genomes
load( "../../genomes/data/2021_03_26_download_patric_genomes.rData" )
load( "../data/eipData.rData" )

# Read in the relevant meta-data to the eip genomes
isolateMetaData = read.table(
  "../data/2022_02_03_isolateMetaData.tsv",
  sep              = '\t',
  header           = TRUE,
  stringsAsFactors = FALSE
  )

# Read in the coganc fastTree
cognacDir  = "../analysis/2022_02_04_eip_cognac_analysis/"
cgTreePath = paste0( cognacDir, "concatenated_gene_aa_alignment_fasttree.tre" )
cgTree     = midpoint.root( read.tree( cgTreePath ) )

# ---- Source Functions --------------------------------------------------------

GetMlst = function( patMlst, isoSp )
{
  if ( !grepl( "MLST", patMlst ) ) return( "Other" )
  if ( grepl( ',', patMlst, fixed = TRUE ) )
  {
    strVals = strsplit( patMlst, ',' )[[1]]
    patMlst    = strVals[ grep( "_1", strVals ) ]
  }
  strVals = strsplit( patMlst, '[.]' )[[1]]
  st = strVals[ length( strVals ) ]
  return( paste0( isoSp, " ST", st ) )
}

GetSp = function(x)
{
  strVals = strsplit( x, ' ' )[[1]][ c( 1, 2 ) ]
  return( paste( strVals, collapse = ' ' ) )
}

MakeOthrMlst = function(x)
{
  strVals = c( strsplit( x, ' ' )[[1]][ c( 1, 2 ) ], "Other" )
  return( paste( strVals, collapse = ' ' ) )
}

GetCounts = function( inVec )
{
  counts = table( inVec )
  counts = counts[ order( counts, decreasing = TRUE ) ]
  return( counts )
}

source( "R/TreeParsingFunctions.R" )

source( "R/FindTransitonEdgeFns.R" )

GetF1Score = function( tp, fp, fn )
{
  return( tp / ( tp + 0.5 * ( fp + fn ) ) )
}

# arVals =  names( mlstCounts ); features = arMlst; defaultVal = "Other"
# arVals = genusVals; features = genusId;  defaultVal = "Other"
GetLabelEdges = function( arVals, features, defaultVal )
{
  # For each edge, find the isolate which occur downstream, to the tips
  brIsoList = lapply( 1:nrow( cgTree$edge ), function(x)
  {
    brIsos = GetOnBranchIsolates( cgTree, x, NULL )
    return( isolateMetaData$genomeNames %in% brIsos )
  })

  # For E
  trEdge = sapply( arVals, function(x)
  {
    trF1Score = sapply( brIsoList, function( isBr )
    {
      isVal = features[ isBr ] == x
      isOth = features[ isBr ] == "Other"
      tp    = sum( isVal )
      nn    = sum( isBr ) - sum( isOth )
      fp    = nn - tp
      fn    = sum( features[ !isBr ] == x )
      return( GetF1Score( tp, fp, fn ) )
    })
    return( which.max( trF1Score ) )
  })

  brVals = rep( defaultVal, length( features ) )
  for ( j in 1:length( trEdge ) )
  {
    cat( j, '\n' )
    brIsos =
      GetOnBranchIsolates( cgTree, trEdge[j], trEdge[ trEdge != trEdge[j] ] )
    for ( i in which( isolateMetaData$genomeNames %in% brIsos ) )
    {
      if ( brVals[i] != defaultVal ) stop( brVals[i], ' != ', arVals[j], '\n' )
      brVals[i] = arVals[j]
    }
  }
  return( brVals )
}

# ---- Use the tree to find the best genus -------------------------------------

isolateMetaData$species =
  sapply( isolateMetaData$species, function(x) ifelse( is.na(x), "Other", x ) )

genusId =
  sapply( isolateMetaData$species, function(x) strsplit( x, ' ' )[[1]][1] )
gCounts = GetCounts( genusId )

genusVals = names( gCounts )[ gCounts > 1 & names( gCounts ) != "Other" ]
for ( i in which( !genusId %in% genusVals ) ) genusId[i] =  "Other"
bestGenus = GetLabelEdges( genusVals, genusId, "Other" )
GetCounts( bestGenus )

isolateMetaData$species = sapply( 1:nrow( isolateMetaData ),
  function(i)
    ifelse( genusId[i] != bestGenus[i], "Other", isolateMetaData$species[i] )
  )

# ---- Use the tree to find the best species -----------------------------------

sps = c( "Klebsiella pneumoniae", "Escherichia coli", "Enterobacter cloacae" )
klebSp = c(
  "Klebsiella variicola", "Klebsiella oxytoca", "Klebsiella aerogenes"
  )
speciesId = sapply( isolateMetaData$species, function(x)
{
  if ( is.na( x ) ) return( "Other" )
  if ( x %in% sps ) return( x )
  if ( x %in% klebSp ) return(x) #return( "Klebsiella sp." )
  if ( x == "Enterobacter aerogenes" ) return( "Klebsiella aerogenes" )
  if ( x == "Enterobacter hormaechei" ) return( "Enterobacter cloacae" )
  if ( grepl( "Citrobacter", x ) ) return( "Citrobacter sp." )
  if ( grepl( "Raoultella", x ) ) return( "Raoultella sp.")
  if ( grepl( "Enterobacter", x ) ) return(x) #return( " sp.")
  return(x)
  #return( "Other" )
})
spCounts = GetCounts( speciesId )
spVals = names( spCounts )[ spCounts > 5 & names( spCounts ) != "Other" ]

for ( i in which( !speciesId %in% spVals ) ) speciesId[i] =  "Other"

bestSp = GetLabelEdges( spVals, speciesId, "Other" )


bestSp = sapply( 1:length( bestSp ), function(i)
{
  if ( strsplit( bestSp[i], ' ' )[[1]][1] == bestGenus[i] ) return( bestSp[i] )
  # if ( bestSp[i] != "Other" ) return( bestSp[i] )
  if ( bestGenus[i] == "Other" ) return( "Other" )
  return( paste( bestGenus[i], "sp." ) )
})

bestSp = sapply( bestSp, function(x)
{
  if ( x == "Klebsiella pneumoniae" ) return( "K. pneumoniae")
  if ( x == "Escherichia coli" )  return( "E. coli")
  if ( x == "Enterobacter cloacae" ) return( "E. cloacae" )
  return( x )
})

bestSpCounts = GetCounts( bestSp )
bestSpCounts

# ---- Look up the assigned MLST -----------------------------------------------

isIn = patricData$Genome.ID %in% isolateMetaData$genomeNames
patricData = patricData[ isIn, ]

# Define the relevant clonal groups
cg258 = c(
  "K. pneumoniae ST258",
  "K. pneumoniae ST512",
  "K. pneumoniae ST11",
  "K. pneumoniae ST418",
  "K. pneumoniae ST340",
  "K. pneumoniae ST554",
  "K. pneumoniae ST437",
  "K. pneumoniae ST3359"
  )
cg15  = c( "K. pneumoniae ST15", "K. pneumoniae ST14" )
st131 = c( "E. coli ST43" )
st307 = c( "K. pneumoniae ST2975" )

mlst = vector( "character", nrow( isolateMetaData ) )

for ( i in seq( nrow( isolateMetaData ) ) )
{
  isEip = eipData$Sample_ID == isolateMetaData$genomeNames[i]
  if ( TRUE %in% isEip )
  {
    eIdx =isEip
    mlst[i] = eipData$ST[ eIdx ]
  } else if ( isolateMetaData$sources[i] == "PATRIC" ) {
    pIdx = which( patricData$Genome.ID == isolateMetaData$genomeNames[i] )
    # if ( patricData$MLST[ pIdx ] != '' ) stop( i)
    mlst[i] = GetMlst( patricData$MLST[ pIdx ], bestSp[i] )
  } else {
    mlst[i] = "Other"
  }

  if ( mlst[i] %in% cg258 ) mlst[i] = "K. pneumoniae CG258"
  if ( mlst[i] %in% cg15 )  mlst[i] = "K. pneumoniae CG15"
  if ( mlst[i] %in% st131 ) mlst[i] = "E. coli ST131"
  if ( mlst[i] %in% st307 ) mlst[i] = "K. pneumoniae ST307"
}

# ---- Use the tree to find the best MLST --------------------------------------

arMlst = sapply( 1:length( bestSp ), function(i)
{
  if ( bestSp[i] %in% c( "K. pneumoniae", "E. coli", "E. cloacae" ) )
  {
    if ( mlst[i] == "Other" ) return( "Other" )
    isSp =
      strsplit( bestSp[i], ' ' )[[1]][2] == strsplit( mlst[i], ' ' )[[1]][2]
    if ( isSp ) return( mlst[i] )
    return( "Other" )
  }
  return( bestSp[i] )
})
mlstCounts = GetCounts( arMlst )
isMlst = mlstCounts > 1 & names( mlstCounts ) != "Other"
mlstVals = names( mlstCounts )[ isMlst ]

bestMlst    = GetLabelEdges( mlstVals, arMlst, "Other" )
bMlstCounts = GetCounts( bestMlst )
bMlstCounts

# ---- Update the meta-data ----------------------------------------------------

isolateMetaData$species = bestSp
isolateMetaData = cbind.data.frame( isolateMetaData, "mlst" = bestMlst )

write.table(
  isolateMetaData,
  "../data/isolateMetaData.tsv",
  sep = '\t'
  )

# ------------------------------------------------------------------------------
#
# # Get the distance from the tree
# distMat = cophenetic.phylo( cgTree )
#
#
# rIdxs = sapply( isolateMetaData$genomeNames,
#   function(x) which( row.names( distMat ) == x )
#   )
#
# mlstIdxs = lapply( mlstVals, function(x) rIdxs[ mlst == x ] )
#
# minStIdxs = lapply( rIdxs, function(i)
# {
#   mDists = sapply( mlstIdxs, function(x) median( distMat[ i, x ] ) )
#   minIdx = which.min( mDists )
#   return( list( minIdx, mDists[ minIdx ] ) )
# })
#
# OUTLIER = 0.001
# bestMlst = sapply( minStIdxs, function(x)
# {
#   stIdx = x[[1]]
#   if ( x[[2]] < OUTLIER ) return( mlstVals[ stIdx ] )
#   return( MakeOthrMlst( mlstVals[ stIdx ] ) )
# })
#
#
# bestMlstIdxs = sapply( unique( bestMlst ),
#   function(x)
#     which( row.names( distMat ) %in% names( bestMlst )[ bestMlst == x ] )
#   )
#
#
# xx = sapply( bestMlstIdxs,
#   function(x) sapply( bestMlstIdxs, function(y) mean( distMat[ x, y ] ) )
#   )
#
# row.names( xx ) = unique( bestMlst )
# colnames( xx )  = unique( bestMlst )
# pheatmap( xx )
#
#
# bestSp = sapply( bestMlst, GetSp )
#
# table( bestSp )
