# ------------------------------------------------------------------------------
# Make co-phylogeny plots
# 2022/02/27
# Ryan D. Crawford
# ------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( ape )
library( maps )
library( cluster )
library( phytools )
library( pheatmap )

# ---- Import the data ---------------------------------------------------------

# Load the paths to the alignments and the fastTrees: "clAlgnPaths", and
# "cgTreePaths"
load( "../data/2022_02_25_make_plasmid_trees.rData" )

# Load the annotations and colors: "annots" and "annotCols"
load( "../data/2022_02_07_make_annotations_and_colors.rData" )

# Read in the core gene fastTree
cgTree = read.tree( "../data/cgTree.tre" ) #midpoint.root()

# "plTrees", "featEdgeMats", "edgeColMats", "plTrMats"
load( "../data/2022_02_28_run_mpr_on_plasmid_trees.rData" )

# "arFeatureMats", "edgeColMat", "featEdgeMat", "edgeFeatList"
load( "../data/2022_02_19_run_mpr_for_location.rData" )

# ---- Function definitions ----------------------------------------------------

source( "R/PlotTree.R" )

#
PlotPlasmidTree = function( tree, annots, annotCols, soCols, facCols )
{
  # Set the layout for the plots
  layoutMat = matrix(c(
    1, 1, 2, 2,
    1, 1, 2, 2,
    1, 1, 2, 2,
    1, 1, 2, 2,
    4, 4, 3, 3,
    4, 4, 3, 3,
    5, 5, 3, 3,
    5, 5, 3, 3
    ), ncol = 4, byrow = TRUE )
  layout( layoutMat )

  # Ploty the tree
  PlotTree( tree, annots, annotCols, facCols$right, 2 )
  PlotSpeciesTree( tree, annots, annotCols, soCols$right, 2 )
  PlotLegend( annotCols )
  par( mar = c( 10, 6.5, 3, 1) )
  MakeTreeBarplot( tree, annots, 1 )
  MakeTreeBarplot( tree, annots, 2 )
}


i = which( names( plTrees ) == 49)
plTree = plTrees[[i]]; figId = trIds[i];
soCols = soEdgeCols[[i]]; strCols = strEdgeCols[[i]];
facCols = facEdgeCols[[i]]
MakeTreePlots = function(
  cgTree, plTree, annots, annotCols, figId, soCols, strCols, facCols
  )
{
  # Generate the plots to output
  figDir     = "../figures/2022_02_27_make_co_phylo_plots/"
  plTreePlot = paste0( figDir, "2022_02_27_", figId, "_pl_tree_plot.png" )
  souPhylo   = paste0( figDir, "2022_02_27_", figId, "_source_co_phylo.png" )
  strCoPhylo = paste0( figDir, "2022_02_27_", figId, "_strain_co_phylo.png" )

  png(
    filename = plTreePlot,
    units    = "in",
    width    = 5.5,
    height   = 9,
    res      = 300
    )
  PlotPlasmidTree( plTree, annots[ , 1:5 ], annotCols[ 1:5 ], soCols, facCols )
  dev.off()

  isEip = annots$Source != "Other"
  png(
    filename = souPhylo,
    units    = "in",
    width    = 5.5,
    height   = 4.75,
    res      = 300
    )
  # PlotCoPhylo( cgTree, plTree, annots[ isEip, ], 3, strCols )
  MakeCoPhyloPlot(
    cgTree, plTree, annots[ isEip, ], annotCols, 3, strCols$left, strCols$right
    )
  dev.off()

  png(
    filename = strCoPhylo,
    units    = "in",
    width    = 5.5,
    height   = 4.75,
    res      = 300
    )
  MakeCoPhyloPlot(
    cgTree, plTree, annots, annotCols, 6, soCols$left, soCols$right
    )
  # PlotCoPhylo( cgTree, plTree, annots, annotCols, 6, soCols, 2 )
  dev.off()


  return( c( plTreePlot, souPhylo, strCoPhylo ) )
}

GetStrainId = function( annots )
{
  strainId = sapply( 1:nrow( annots ), function(i)
  {
    if ( annots$ST[i] != "Other" ) return( annots$ST[i] )
    return( annots$Species[i] )
  })
  return( strainId )
}

# ---- Make the figures --------------------------------------------------------

annots     = cbind( annots, "strain" = GetStrainId( annots ) )
strainCols =
  list( "strain" = c( annotCols$ST[ 1:4 ], annotCols$species[ 1:5 ] ) )
annotCols  = append( annotCols, strainCols )

# plTrees = lapply( cgTreePaths, function(x) read.tree( x ) )

nTrees = length( plTrees )
trIds  = paste0( "cl_", names( plTrees ) )

strEdgeCols = lapply( seq( nTrees ),
  function(i)
    list( "left" = edgeColMat[ , 1 ], "right" = edgeColMats[[i]][ , 1 ] )
  )

soEdgeCols = lapply( seq( nTrees ),
  function(i)
    list( "left" = edgeColMat[ , 2 ], "right" =  edgeColMats[[i]][ , 2 ] )
  )

facEdgeCols = lapply( seq( nTrees ),
  function(i)
    list( "left" = edgeColMat[ , 3 ], "right" =  edgeColMats[[i]][ , 3 ] )
  )

# soCols = soEdgeCols[[i]]; strCols = strEdgeCols[[i]]; facCols = facEdgeCols[[i]]
trPlotPaths = lapply( seq( nTrees ), function(i)
{
  cat( i, '\n' )
  MakeTreePlots(
    cgTree, plTrees[[i]], annots, annotCols, trIds[i],
    soEdgeCols[[i]], strEdgeCols[[i]], facEdgeCols[[i]]
    )
})


# ---- Generate the Rmarkdown --------------------------------------------------

sink( "2022_02_27_make_co_phylo_plots.Rmd" )

cat(
"---
title: \"Plasmid co-phylogeny plots\"
author: \"Ryan D. Crawford\"
date: \"02/27/2022\"
output: pdf_document
---\n\n\n"
)

for ( i in seq( nTrees ) )
{
  cat( "## Gene cluster", i, "concatenated gene tree\n\n" )
  cat( paste0( "![](", trPlotPaths[[i]], ")\n\n", collapse = '' ) )
}
sink()

# ------------------------------------------------------------------------------