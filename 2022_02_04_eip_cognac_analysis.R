# ------------------------------------------------------------------------------
# EIP cognac analysis
# 2022/02/04
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# Run cognac to generate the core gene alignment for the EIP genomes and
# matched public isolates.
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( cognac )

# ---- Load the data -----------------------------------------------------------

# Read in the relevant metas-data to the eip genomes
isolateMetaData = read.table(
  "../data/isolateMetaData.tsv",
  sep              = '\t',
  header           = TRUE,
  stringsAsFactors = FALSE
  )

# ---- Run cognac -------------------------------------------------------------

# Set the output directory to write the alignment and temp files
outDir = "../analysis/2022_02_04_eip_cognac_analysis/"

# Parse the genomic data
geneEnv = CreateGeneDataEnv(
  isolateMetaData$gffFiles,
  isolateMetaData$fastaFiles,
  isolateMetaData$genomeNames,
  outDir
  )

# Generate the alignment
algnEnv = cognac(
  geneEnv       = geneEnv,
  njTree        = TRUE,
  keepTempFiles = TRUE,
  mapNtToAa     = TRUE,
  outDir        = outDir
  )

# ---- Save the data -----------------------------------------------------------

print( pryr::mem_used() )
rData  = "../data/2022_02_04_eip_cognac_analysis.rData"
save( file = rData, list = ls() )

# ------------------------------------------------------------------------------
