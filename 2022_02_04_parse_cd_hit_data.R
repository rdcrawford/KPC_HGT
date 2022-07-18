# ------------------------------------------------------------------------------
#
# 2021/05/10
# Ryan D. Crawford
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( Rcpp )
library( cognac )
library( cluster )
library( future )
library( future.apply )
library( ape )

# ---- Constant declarations ---------------------------------------------------

# Column indices of the gff files
G_ID    = 1
G_ANNOT = 2
CONTIG  = 3
G_START = 4
G_END   = 5
G_CLUST = 7

# ---- Load the data -----------------------------------------------------------

# Load the data from the cognac run
load( "../data/2022_02_04_eip_cognac_analysis.rData" )
print( pryr::mem_used() )

# ---- Function definitions ----------------------------------------------------

sourceCpp( "R/ParseCdHit.cpp" )

# Add the cd-hit cluster name to the gff file
AddClustIdxToGff = function( rIdx )
{

  PrintPercComp( rIdx )

  # Initiaize a vector to store the cluster id for each gene
  clIdxs = vector( "integer", nrow( geneEnv$gfList[[ rIdx ]] ) )

  # For each cd-hit cluster...
  for ( i in  seq( nGenes ) )
  {
    # Find if the current genome has this gene
    isThisGenome = geneEnv$genomeIdList[[i]] == geneEnv$genomeNames[ rIdx ]

    # Match the gene id to the corresponding row in the gff file, and
    # mark the corresponding cluster as this gene index
    for ( j in which( isThisGenome ) )
    {
      isClGene =
        geneEnv$gfList[[ rIdx ]][ , G_ID ] == geneEnv$clustList[[i]][j]
      clIdxs[ isClGene ] = i
    }
  }

  # Add cluster indexs to the data frame
  gff = cbind.data.frame(
    geneEnv$gfList[[ rIdx ]],
    clIdxs,
    stringsAsFactors = FALSE
    )

  return( list( gff ) )
}

PrintPercComp = function( n )
{
  if ( !n %% 1000 ) cat( n, ":", round( n / nGenomes * 100, 2 ), "\n" )
}

# ---- Parse cd-hit results ----------------------------------------------------

# Parse the cd-hit cluster data and add it to the environment.
clstrFile = paste0( outDir, "temp_cognac_files/cdHitClusters.faa.clstr" )

# Parse the cd-hit results
minGeneNum = 2
ParseCdHitToList( clstrFile, minGeneNum, geneEnv )

# Add a cluster index to the gff files
nGenes         = length( geneEnv$clustList )
nGenomes       = length( geneEnv$genomeNames )
geneEnv$gfList = sapply( seq( nGenomes ), AddClustIdxToGff )

# ---- Save the data -----------------------------------------------------------

print( pryr::mem_used() )
rData = "../data/2022_02_04_2022_02_04_parse_cd_hit_data.rData"
save( file = rData, list = "geneEnv" )

# ------------------------------------------------------------------------------