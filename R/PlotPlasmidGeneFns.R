PlotCongigGenes = function( geneMats, hmTitle, path )
{
  # 
  if ( missing( hmTitle ) ) hmTitle = NA
  if ( missing( path ) )    path    = NA
  
  # Find the number of genes on each contig
  geneCount =  sapply( geneMats, ncol )
  
  # Initialize the gene presence mat with the first matirx
  annotGeneMat = geneMats[[1]]
  
  # If there is more than one congig with the gene annotation, 
  # append it/them to the first gene presence absence matrix
  if ( length( geneMats ) > 1 )
  {
    for ( i in 2:length( geneMats) ) 
      annotGeneMat = cbind( annotGeneMat, geneMats[[i]] )
    
    # Find the gaps in each plasmid sequence 
    gaps = c(
      geneCount[1],
      sapply( 2:length( geneCount ), function(i) sum( geneCount[ 1:i ] ) )
    )
    
  } else {
    
    gaps = NULL
  }
  
  # Plot a heatmap
  pheatmap(
    annotGeneMat, 
    col               = c( "grey", "navyblue" ),
    gaps_col          = gaps,
    annotation_row    = annots,
    # annotation_colors = annotCols,
    show_rownames     = FALSE,
    show_colnames     = FALSE,
    cluster_cols      = FALSE, 
    cluster_rows      = FALSE,
    main              = hmTitle,
    filename          = path
  )
}



MakeGenePresenceHeatmap = function( contigData, annot )
{  
  genes = unique( unlist( contigData ) )
  geneMat = sapply( 1:length( contigData), 
                    function(i) as.numeric( genes %in% contigData[[i]] )
  )
  geneMat = geneMat[ rowSums( geneMat ) > 1, ]
  
}