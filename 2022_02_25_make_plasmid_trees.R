# ------------------------------------------------------------------------------
# Make Plasmid gene trees
# 2022/02/25
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# For each plasmid gene cluster, generate a concatenated gene alignment with
# the respective nucleotide sequences for the cluster genes. Then generate
# an ~ML tree with fast tree.
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( Rcpp )
library( cognac )

# ---- Function definitions ----------------------------------------------------

GetRefIdx = function( distMat, refGenes )
{
  isGene = which( row.names( distMat ) %in% refGenes )
  if ( length( isGene ) == 1 ) return( isGene )
  return( isGene[ which.min( rowSums( distMat[ isGene, ] ) ) ] )
}

ParseAlgn = function( algn, distMat, algnGenomes )
{
  # Make the alignment to
  outAlgn = vector( "character", length( algnGenomes ) )
  names( outAlgn ) = algnGenomes

  # Get the genome names of each sequence in the alignment
  genomeIds = sapply( names( algn ), GetGenomeId )

  # Get the length of the alignment
  algnLen = nchar( algn[1] )

  # For each genome, look up the best aligned sequence
  for ( i in 1:length( algnGenomes ) )
  {
    gIdx = which( genomeIds == algnGenomes[i] )

    # If there is one copy of this gene, return the aligned sequence
    if ( length( gIdx ) == 1 )
    {
      outAlgn[i] = algn[ gIdx ]

    # If this gene is missing, return gaps
    } else if ( length( gIdx ) == 0 ) {

      outAlgn[i] = paste( rep( '-', algnLen ), collapse = '' )

    # If there is more than one copy of this gene minimize the distance to
    # the reference allele
    } else {

      cat( length( gIdx ), '\n' )
      gIdx = gIdx[ which.min( rowSums( distMat[ gIdx,  ] ) ) ]
      outAlgn[i] = algn[ gIdx ]
    }
  }
  return( outAlgn )
}

ConcatAlignGenes = function( clAlgns, clDistMats, algnGenomes, clAlgnPath )
{
  # Concatenate the alignments by appending each alignment to the vector
  # of the concatenated alignemnt (passed by reference)
  concatAlgn = vector( "character", length( algnGenomes ) )
  for ( i in 1:length( clAlgns ) )
  {
    algn = ParseAlgn( clAlgns[[i]], clDistMats[[i]], algnGenomes )
    ConcatenateAlignments( concatAlgn, algn )
  }

  # Write the fasta file to the input path
  WriteFa( clAlgnPath, concatAlgn, algnGenomes )
}

SubmitFastTree = function( algnPath, trePath, nGs )
{
  faName  = ExtractGenomeNameFromPath( algnPath )
  projDir = "/nfs/turbo/umms-esnitkin/ryan/eip_hgt_analysis/"
  workDir = paste0( projDir, "analysis/2021_02_05_make_plasmid_trees/" )
  shDir   = paste0( workDir, "slurm_scripts/" )

  if ( nGs <= 10 )
  {
    memVal = 8
  } else if ( nGs <= 100 ) {
    memVal = 24
  } else {
    memVal = 64
  }

  partiton = "standard"

  sbat = c(
    "#!/bin/sh",
    paste0( "#SBATCH --job-name=", faName ),
    "#SBATCH --mail-user=rcrawfo@umich.edu",
    "#SBATCH --mail-type=BEGIN,END,FAIL",
    "#SBATCH --cpus-per-task=1",
    "#SBATCH --nodes=1",
    "#SBATCH --ntasks-per-node=1",
    paste0( "#SBATCH --mem-per-cpu=", memVal, "gb" ),
    "#SBATCH --time=1-00:00:00",
    "#SBATCH --account=esnitkin1",
    paste0( "#SBATCH --partition=", partiton ),
    paste0( "#SBATCH --output=", shDir, "%x-%j.log" ),
    paste0( "cd ", projDir, "lib/" ),
    paste0( "fasttree -nt ", algnPath, " > ", trePath )
    )

  sbatPath = paste0( shDir, faName, ".sh" )
  sink( sbatPath )
  cat( sbat, sep = '\n' )
  sink()

  system( paste( "sbatch", sbatPath ) )
}

source( "R/WriteFa.R" )

# ---- Import the data ---------------------------------------------------------

# Load the parsed nt alignments includes: "algns", "geneDists", "kpcClList"
load( "../data/2022_02_25_get_nt_cluster_algn_data.rData" )

# Load the data on the plasmid clusters: "eipClMat", "eipAcqEdge", "plMbrList",
# and "plMat"
load( "../data/2022_02_16_analyze_hgt_cluster_distribution.rData" )

# "plTrees", "featEdgeMats", "edgeColMats", "plTrMats"
load( "../data/2022_02_28_run_mpr_on_plasmid_trees.rData" )

# ----- For each gene cluster, generate an alignment and tree ------------------

algnDir =
  "../analysis/2021_02_05_make_plasmid_trees/concatenated_gene_alignments/"

# Create a list of the isolates which are positive for each cluster
algnGenomes =
  lapply( 1:ncol( plMat ), function(j) row.names( plMat )[ plMat[ , j ] == 1 ] )

# Generate the paths to the alignments to create
nCls        = length( plMbrList )
clAlgnPaths = paste0( algnDir, "cl_", seq( nCls ), "_cg_algn.fasta" )

for ( i in seq( nCls ) )
{
  clGenes = names( algns ) %in% plMbrList[[i]]

  ConcatAlignGenes(
    algns[ clGenes ],
    geneDists[ clGenes ],
    algnGenomes[[i]],
    clAlgnPaths[i]
    )
}

coreAlgnDir =
  "../analysis/2021_02_05_make_plasmid_trees/core_concatenated_gene_alignments/"
coreAlgnPaths = paste0( algnDir, "cl_", seq( nCls ), "_cg_algn.fasta" )

# Filter to core alignment positions
for ( i in 1:length( coreAlgnPaths ) )
  FilterAlgnPositions( clAlgnPaths[i], coreAlgnPaths[i], 0.1, 2 )

# Generate the paths to the trees
treeDir =
  "../analysis/2021_02_05_make_plasmid_trees/core_concatenated_gene_trees/"
cgTreePaths = paste0( treeDir, "cl_", seq( nCls ), "_cg_fastTree.tre" )

# Submit the jobs to make fast trees from the alignments
nGs = sapply( algnGenomes, length )
for ( i in seq( nCls ) )
  SubmitFastTree( coreAlgnPaths[i], cgTreePaths[i], nGs[i] )

# ---- Save the data -----------------------------------------------------------

print( pryr::mem_used() )
rData = "../data/2022_02_25_make_plasmid_trees.rData"
save( file = rData, list = c( "coreAlgnPaths", "clAlgnPaths", "cgTreePaths" ) )

# ------------------------------------------------------------------------------