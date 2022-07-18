# ------------------------------------------------------------------------------
# Algn Unique Seqs
# 2019/12/02
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# This function identifies unique sequences and aligns them using mafft
# ------------------------------------------------------------------------------

AlgnUniqueSeqs = function(
  mafftInPath, # Fasta file with the full data-set of sequences
  algnPath,    # Path  to the output alignment
  algnDir,     # Directory to write the alignment and input file
  maffftOps    # Optional. Argumnets for mafft. mafft [maffftOps] in > out
  )
{
  # Read in the entire number of sequences 
  algnSeqs = ParseFasta( mafftInPath )
  
  # For each unique sequence, find the identical sequences and
  # store them in a list.
  identList = FindIdenticalGenes( algnSeqs, names( algnSeqs ) )

  cat( 
    "  -- ", length( algnSeqs ), " sequences were input\n",
    "  -- ", length( identList ), " sequences are unique\n",
    sep = ''
    )
  
  # If there is no variation in the gene, it doesnt matter. Just return an
  # empty character vector
  if ( length( identList ) == 1 ) stop( "All sequences are identical..." )

  # If there are no duplication, align all of the genes
  isGeneRep  = names( algnSeqs ) %in% names( identList )
  toAlgnIdxs = which( isGeneRep )

  # Create the path to the unique sequence
  uniqueSeqsPath = paste0( algnDir, "temp_unique_seqs_mafft_input.fasta" )
  
  # Write the input fasta file
  sink( uniqueSeqsPath )
  for ( i in toAlgnIdxs )
    cat( '>', names( algnSeqs )[i], '\n', algnSeqs[i], '\n', sep = '' )
  sink()

  # Generate the mafft command
  if ( missing( maffftOps ) )
    maffftOps = "--retree 2 --maxiterate 2 --quiet --leavegappyregion"
  mafftCmd = paste( "mafft", maffftOps, uniqueSeqsPath, '>', algnPath )

  # Run Mafft
  system( mafftCmd )
  
  # Read in the alignment
  algn = ParseFasta( algnPath )

  # If there were any duplicated genes
  for ( i in 1:length( identList ) )
  {
    nGenes = length( identList[[i]] )
    if ( nGenes )
    {
      algnVecLen = length( algn ) + 1
      newSeqRange = algnVecLen:( algnVecLen + nGenes - 1 )
      algn[ newSeqRange ] = rep( algn[i], nGenes )
      names( algn )[ newSeqRange ] = identList[[i]]
    }
  }

  # Write the full alignment with the redundant sequnces filled in
  sink( algnPath )
  for ( i in 1:length( algn ) )
    cat( '>', names( algn )[i], '\n', algn[i], '\n', sep = '' )
  sink()
}

# ------------------------------------------------------------------------------

