# ------------------------------------------------------------------------------
# Run Mafft functions
# 2022/02/04
# Ryan D. Crawford
# ------------------------------------------------------------------------------

# Write the fasta file with the identifier for the sequence in the output
# directory with the extension fasta
WriteFa = function( outDir, seqId, seqs )
{
  # Create the path to write to
  faPath = paste0( outDir, seqId, ".fasta" )

  # Write the input fasta file
  sink( faPath )
  for ( i in 1:length( seqs ) )
    cat( '>', names( seqs )[i], '\n', seqs[i], '\n', sep = '' )
  sink()

  return( faPath )
}

# Run mafft for the fasta file passed in. Return the path to the alignment
RunMafft = function( outDir, seqId, faPath, mafftOpts )
{
  # Create the path to write to
  algnPath = paste0( outDir, seqId, "_mafft_alignment.fasta" )

   # Generate the mafft command
  if ( missing( mafftOpts ) )
    mafftOpts = "--retree 2 --maxiterate 2 --quiet"

  mafftCmd = paste( "mafft", mafftOpts, faPath, '>', algnPath )

  # Run Mafft
  system( mafftCmd )

  return( algnPath )
}

# For each set of sequences in the input list, generate an alignment
MakeGeneAlgns = function( seqList, outDir )
{
  algnPaths = sapply( 1:length( seqList ), function(i)
  {
    cat( i, ": Aligning ", length( seqList[[i]] ), " sequences\n", sep = '' )
    return( AlgnGeneSeqs( outDir, names( seqList )[i], seqList[[i]] ) )
  })
  names( algnPaths ) = names( seqList )
  return( algnPaths )
}

AlgnGeneSeqs = function( outDir, seqId, seqs, mafftOpts )
{
  # For each unique sequence, find the identical sequences and
  # store them in a list.
  identList = FindIdenticalGenes( seqs, names( seqs ) )

  # If there is no variation in the gene, it doesn't matter. Just return an
  # empty character vector
  if ( length( identList ) == 1 ) return( WriteFa( outDir, seqId, seqs ) )

  # Set the path to the output files for the genes and the gene alignment
  isAlgn   = names( seqs ) %in% names( identList )
  faPath   = WriteFa( outDir, seqId, seqs[ isAlgn ] )
  algnPath = RunMafft( outDir, seqId, faPath, mafftOpts )

  # Parse the input gene alignment
  algn = ParseGeneAlgn( algnPath, identList )

  # Return the path to the full alignment
  algnId = paste0( seqId, "_all_seqs_mafft_alignnt" )
  return( WriteFa( outDir, algnId, algn ) )
}

ParseGeneAlgn = function( algnPath, identList )
{
  # Read in the alignment
  algn = ParseFasta( algnPath )

  # If there were any duplicated genes
  for ( i in 1:length( identList ) )
  {
    nGenes = length( identList[[i]] )
    if ( nGenes )
    {
      algnVecLen  = length( algn ) + 1
      newSeqRange = algnVecLen:( algnVecLen + nGenes - 1 )
      algn[ newSeqRange ] = rep( algn[i], nGenes )
      names( algn )[ newSeqRange ] = identList[[i]]
    }
  }
  return( algn[ order( names( algn ) ) ] )
}

# ------------------------------------------------------------------------------
