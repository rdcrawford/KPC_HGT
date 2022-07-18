# ------------------------------------------------------------------------------
# FilterRecConcatAlgnGenes
# 2021/08/19
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# Parse the gubbins data. Identify any genes that were identifed as 
# recombinant, and mask the gene in the corresponding genes in the alignment. 
# ------------------------------------------------------------------------------

ReadGubbins = function( gffPath )
{
  gff = read.table(
    gffPath,
    sep ='\t',
    stringsAsFactors = FALSE,
    comment.char = '#',
    header = FALSE
    )
  gffListElm = list( gff )
 
  return( gffListElm )
}

MakeGeneRecMat = function( gIds, recLociList, algnPartitions )
{
  # make the matrix
  nGenomes = length( gIds )
  nGenes   = nrow( algnPartitions )
  recMat   = matrix( 0, nrow = nGenomes, ncol = nGenes )
  row.names( recMat ) = gIds
  
  # Get the genes involved in each recombination event
  geneIdxs = GetLocusGenes( recLociList, algnPartitions )
  
  # 
  for ( k in 1:length( recLociList ) )
  {
    isGenome = which( gIds %in% recLociList[[k]][[3]] )
    
    for ( j in geneIdxs[[k]] )
    {
      for ( i in isGenome ) recMat[ i, j ] = recMat[ i, j ] + 1
    }
  }
  return( recMat )
}

GetLocusGenes = function( recList, algnPartitions )
{
  nLoci    = length( recList )
  nGenes   = nrow( algnPartitions )
  gIdx     = 1
  recGenes = vector( "list", nLoci )
  
  for ( rIdx in seq( nLoci ) )
  {
    # 
    while( recList[[rIdx]][[1]] > algnPartitions[ gIdx, 1 ] ) 
      gIdx = gIdx + 1
    
    # 
    endGIdx = gIdx
    while( recList[[rIdx]][[2]] > algnPartitions[ endGIdx, 2 ] ) 
      endGIdx = endGIdx + 1
    
    # 
    recGenes[[ rIdx ]] = seq( gIdx, endGIdx )
  }
  return( recGenes )
}

MakeGenomeRecList = function( gubbinsGffs )
{
  nRecLoci = sum( sapply( gubbinsGffs, nrow ) )
  
  recList = vector( "list", nRecLoci )
  idx = 1
  for ( j in 1:length( gubbinsGffs ) )
  {
    gff = gubbinsGffs[[j]]
    for ( i in 1:nrow( gff ) )
    {
      gIds = GetGenomeIdsFromDescr( gff[ i, 9 ] )
      recList[[ idx ]] = list( gff[ i, 4], gff[ i, 5], gIds )
      idx = idx + 1
    }
  }
  return( recList )
}

GetGenomeIdsFromDescr = function( descr )
{
  elms   = strsplit( descr, ';' )[[1]]
  taxIdx = grep( "taxa=", elms )
  taxa   = gsub( "taxa=", '',  elms[ taxIdx ] )
  return( strsplit( taxa, ' ' )[[1]] )
}

SortLocusList = function( recList )
{
  return( recList[ order( sapply( recList, function(x) x[[1]][[1]] ) ) ] )
}

GetAlgnPartitions = function( cognacData )
{
  algnData = read.table(
    cognacData,
    sep = '\t',
    stringsAsFactors = FALSE,
    header = TRUE
    )

  algnPartitions = t(sapply( algnData$ntGenePartitions, 
    function(x) as.integer( strsplit( x, '-' )[[1]] ), USE.NAMES = FALSE
    ))
  
  return( algnPartitions )
}

sourceCpp( "R/MaskRecPositions.cpp" )

MaskRecAlgnGenes = function( algn, algnPartitions, recMat )
{
  seqIds = names( algn )
  for ( i in 1:length( algn ) )
  {
    rIdx = which( row.names( recMat ) == names( algn )[i] )
    
    isRec = recMat[ rIdx, ] != 0
    
    if ( TRUE %in% isRec )
    {
      algn[i] = MaskRecPositions(
        algn[i],
        algnPartitions[ isRec, 1 ],
        algnPartitions[ isRec, 2 ]
        )
    }
  }
  names( algn ) = seqIds 
  return( algn )
}

WriteAlgn = function( algn, path )
{
  sink( path )
  for ( i in 1:length( algn ) ) 
    cat( '>', names( algn )[i], '\n', algn[i], '\n', sep = '' ) 
  sink()
}

FilterRecConcatAlgnGenes = function( cognacData, algnPath, gffPaths )
{
  # Read in the alignment
  algn = ParseFasta( algnPath )
  
  # Parse the cognac alignment meta-data to get the partitions
  algnPartitions = GetAlgnPartitions( cognacData )
  
  # Read in the gubbins recombination predictions 
  gubbinsGffs = sapply( gffPaths, ReadGubbins, USE.NAMES = FALSE )
  
  # Make a list of the loci identified as recombinant
  recLociList = MakeGenomeRecList( gubbinsGffs )
  
  # Make genome x gene matrix indicating if a gene is recombinant
  recMat = MakeGeneRecMat( names( algn ), recLociList, algnPartitions )
  
  # Mask the alignment positions in the alignment
  algn = MaskRecAlgnGenes( algn, algnPartitions, recMat )
  
  # Write the mask alignment to the new path
  maskAlgnPath = gsub( ".fasta", "_rec_filtered.fasta", algnPath )
  WriteAlgn( algn, maskAlgnPath )
  
  return( maskAlgnPath )
}

# ------------------------------------------------------------------------------