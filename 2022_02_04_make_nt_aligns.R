# ------------------------------------------------------------------------------
# Make NT aligns
# 2022/02/04
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# Extract the neucleotide sequences for each gene in the analysis. Then for
# each cd-hit cluster, which represent orthologous genes, generate an alignment
# with mafft.
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( cognac )
library( Rcpp )

# ---- Load the data -----------------------------------------------------------

# Load the "geneEnv" with the parsed gff files and cd-hit results
load( "../data/2022_02_04_2022_02_04_parse_cd_hit_data.rData" )

# ---- Function definitions ----------------------------------------------------

source( "R/ExtractNtSeqs.R" )

source( "R/RunMafftFns.R" )

sourceCpp( "R/GetNtSeqs.cpp" )

MakeSeqList = function( geneEnv )
{
  if ( length( geneEnv$geneIds ) != length( geneEnv$ntSeqs ) )
    stop(
      "n gene ids: ", length( geneEnv$geneIds ),
      " n nt seqs: ", length( geneEnv$ntSeqs ), '\n'
      )

  nGenes = length( geneEnv$clustList )
  seqList = lapply( seq( nGenes ), function(i)
  {
    isCl = geneEnv$geneIds %in% geneEnv$clustList[[i]]
    seqs = geneEnv$ntSeqs[ isCl ]
    names( seqs ) = geneEnv$geneIds[ isCl ]
    return( seqs )
  })

  names( seqList ) = paste0( "gene_cluster_", seq( nGenes ) )
  return( seqList )
}

# ---- Make the gene alignments ------------------------------------------------

# Parse the fasta files to get the neucleotide gene sequences
ExtractNtSeqs( geneEnv )

# Make a list of the sequences to align
seqList = MakeSeqList( geneEnv )

# Run Mafft on the nt sequences
outDir    = "../analysis/2022_02_04_eip_cognac_analysis/nt_gene_algns/"
algnPaths = MakeGeneAlgns( seqList, outDir )

# ---- Save the data -----------------------------------------------------------

print( pryr::mem_used() )
rData = "../data/2022_02_04_make_nt_aligns.rData"
save( file = rData, list = ls() )

# ------------------------------------------------------------------------------