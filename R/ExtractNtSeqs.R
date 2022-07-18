# ------------------------------------------------------------------------------
# ReverseTranslateAlgn
# 2020//13
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# This function inputs a concatenated gene alignment
# ------------------------------------------------------------------------------

ExtractNtSeqs = function( geneEnv )
{
  # Column index of the gene Ids in the gff file
  ANNOT_ID = 1

  # For each observation in the concatenated alignment, look up the sequences
  # in the gff file
  nGenomes = length( geneEnv$genomeNames )
  ntSeqList = lapply( seq( nGenomes ), function(i)
  {
    PrintStatus(i)
    genomeSeq     = ParseFasta( geneEnv$fastaFiles[i] )
    seqs          = GetNtSeqs( geneEnv$gfList[[i]], genomeSeq )
    names( seqs ) = geneEnv$gfList[[i]][ , ANNOT_ID ]
    return( seqs )
  })

  geneEnv$ntSeqs = unlist( ntSeqList )
  geneEnv$geneIds = names( geneEnv$ntSeqs )
}

PrintStatus = function(i)
{
  if ( !i %% 100 )
  {
    perc = round( i / length( geneEnv$genomeNames ) * 100, 2 )
    cat( i, ": ", perc, "%\n", sep = '' )
  }
}


# ------------------------------------------------------------------------------
