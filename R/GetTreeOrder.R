GetTreeOrder = function( tree, genomeIds )
{
  # Get the order of the tips
  isTip    = tree$edge[ , 2 ] <= length( tree$tip.label )
  tipOrder = tree$edge[ isTip, 2 ]

  orderedTips = tree$tip.label[ tipOrder ]
  orderedTips = orderedTips[ orderedTips %in% genomeIds ]

  treeOrder = sapply( orderedTips,
    function(x) which( genomeIds == x )
    )
  return( treeOrder )
}


