## Give a circRNA count data, plot a certain row or plot all rows
#' @title Circ.ratioplot
#' 
#' @description
#' Plot the ratio of circRNAs to the total transcripts (circRNA plus liner RNA).
#' @param Circ CircRNACount file. A file of circRNA read count table. First three columns are circRNA coordinates, and followed by columns for circRNA read counts, each sample per column. 
#' @param Linear LinearCount file. A file of circRNA host gene expression count table. Same configuration as CircRNACount file.
#' @param plotrow The rownumber or rowname in the CircRNACount table corresponding the specific gene which you want to plot.
#' @param size the text size. Default 18.
#' @param ncol if groupindicator2 is provided, specify the panel layout. Default 2.
#' @param CircCoordinates BED format circRNA coordinates file. 
#' @param groupindicator1 A vector of group indicators. For example, ages. 
#' @param groupindicator2 An other vector of group indicators. For example, tissues. This indicator will be used to segement plots out.
#' @param x x axis lable. Default 'Conditions'.
#' @param y y axis lable. Default 'circRNA/(circRNA+Linear)'.
#' @param lab_legend legend title
#' @param circle_description Column indices which do not carry circle/linear read counts.
#' @param gene_column Column index of the column containing the gene name in CircCoordinates if available, otherwise its chosen from Circ.
#' @param y_axis_range Range of the y axis to keep all plots withint he same range
#' @examples data(Coordinates)
#' data(Circ)
#' data(Linear)
#' @examples data(Coordinates)
#' data(Circ)
#' data(Linear)
#' Circ.ratioplot(Circ,Linear,Coordinates,plotrow=10,groupindicator1=c(rep('1days',6),rep('4days',6),rep('20days',6)),groupindicator2=rep(c(rep('Female',4),rep('Male',2)),3),lab_legend='Ages', circle_description = c(1:3), gene_column = 4, y_axis_range = 1)
#' 
#' @export Circ.ratioplot
#' 
Circ.ratioplot <- function(Circ,Linear,CircCoordinates = None,plotrow='1',size=24,ncol=2,groupindicator1=NULL,groupindicator2=NULL,x='Conditions',y='circRNA/(circRNA+Linear)',lab_legend='groupindicator1', circle_description = c(1:3), gene_column = None, y_axis_range = 1, colour_mode = "colour"){

  if( !is.null(groupindicator1) & length(groupindicator1) != ncol(Circ)-length(circle_description) ){
    stop("If provided, the length of groupindicator1 should be equal to the number of samples.")
  }
  if( !is.null(groupindicator2) & length(groupindicator2) != ncol(Circ)-length(circle_description) ){
    stop("If provided, the length of groupindicator2 should be equal to the number of samples.")
  }
  if(is.null(groupindicator1)){
    stop("At least one grouping should be provided through groupindicator1.")
  }
  if(!is.null(groupindicator2)){
    twolevel <- TRUE
  }else{
    twolevel <- FALSE
  }
  
  rownames.circ <- rownames(Circ)
  Circ <- data.frame(lapply(Circ, as.character), stringsAsFactors=FALSE)
  rownames(Circ) <- rownames.circ
  
  rownames.linear <- rownames(Linear)
  Linear <- data.frame(lapply(Linear, as.character), stringsAsFactors=FALSE)
  rownames(Linear) <- rownames.linear
  
  if(!missing(CircCoordinates)){
    rownames.circCoordinates <- rownames(CircCoordinates)
    CircCoordinates <- data.frame(lapply(CircCoordinates, as.character), stringsAsFactors=FALSE)
    rownames(CircCoordinates) <- rownames.circCoordinates
  }else{
    CircCoordinates <- data.frame(Circ[,circle_description])
    rownames(CircCoordinates) <- rownames.circ
    rownames.circCoordinates <- rownames(CircCoordinates)
    CircCoordinates <- data.frame(lapply(CircCoordinates, as.character), stringsAsFactors=FALSE)
    rownames(CircCoordinates) <- rownames.circCoordinates  
  }
  
  groupindicator1 <- factor(groupindicator1,levels=unique(groupindicator1))
  groupindicator2 <- factor(groupindicator2,levels=unique(groupindicator2))
  
  # Get gene name, if no annotation, output NULL
  if (is.character(plotrow)){
    if ( ! plotrow %in% rownames(CircCoordinates) ){
      stop("Specified 'plotrow' not found.")
    }
  }else{
    if ( is.numeric(plotrow) ){
      if ( ! plotrow %in% 1:nrow(CircCoordinates) ){
        stop("Specified 'plotrow' not found.")
      }
    }else{
      stop("Specified plotrow should be ONE rowname or ONE rownumber.")
    }
  }
  # Choose your own column containing the gene name using gene_column. The genename will be displayed in the plot title if available
  if (missing(gene_column)){
    genename = NULL
  }else{
    genename <- as.character(CircCoordinates[plotrow,gene_column])
    if (genename == '.'){
      genename = NULL
    }
  }
  if(twolevel){
    plotdat <- summarySE( data.frame(Ratio=as.numeric(Circ[plotrow,-circle_description])/(as.numeric(Linear[plotrow,-circle_description])+as.numeric(Circ[plotrow,-circle_description])),
                                    groupindicator1,
                                    groupindicator2),
                         measurevar='Ratio',groupvars=c('groupindicator1','groupindicator2') )
  }else{
    plotdat <- summarySE( data.frame(Ratio=as.numeric(Circ[plotrow,-circle_description])/(as.numeric(Linear[plotrow,-circle_description])+as.numeric(Circ[plotrow,-circle_description])),
                                     groupindicator1),
                                     measurevar='Ratio',groupvars=c('groupindicator1') )
  }
# construct plot
  Q <- ggplot(plotdat, aes(x=groupindicator1, y=Ratio)) +
       geom_boxplot() + theme_classic() +
       theme(axis.text.x = element_blank())+
       theme(axis.text.y = element_text(size=size+4))+
       theme(axis.ticks = element_line(colour = 'black', size = 1)) +
       theme(axis.ticks.x = element_blank())+
       theme(legend.title=element_blank()) + 
       theme(text=element_text(size=size+4))+
       theme(legend.text=element_text(size=size)) +
       theme(plot.title = element_text(size=size)) + 
       theme(axis.text.y = element_text(margin=margin(5,5,10,5,"pt")))+
       #labs(list(title=paste("Annotation: ", genename, "\nChr ", toString(Circ[plotrow,circle_description]),sep=""),x=x,y=y)) +
       ggtitle(paste("Annotation: ", genename, "\nChr ", toString(Circ[plotrow,circle_description]),sep="")) +
       ylab("circRNA/(circRNA + Linear RNA)") + 
       xlab("Sample") +         
       geom_errorbar(aes(ymin=Ratio, ymax=Ratio+se), width=.2 , size=2) +
       geom_bar(stat="identity",aes(fill=groupindicator1), color = "black", size=2)

  if (colour_mode == "bw"){
      Q <- Q + scale_fill_grey(start = 0.0, end = 1)
  } else {
      Q <- Q + scale_fill_discrete(name=lab_legend)
  }

       Q <- Q +
       theme(legend.position="bottom") +
       theme(axis.ticks.length = unit(0.5, "cm")) +
       theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=3)) + 
        guides(fill=guide_legend(
                 keywidth=0.3,
                 keyheight=0.3,
                 default.unit="inch")
      ) + scale_y_continuous(expand=c(0,0), limits= c(0, y_axis_range))

  if(twolevel){
    Q <- Q + facet_wrap( ~ groupindicator2,ncol=ncol )
  }

  print(Q)
}

