#' @title Circ.test
#' 
#' @description
#' Test the independent variation of circRNAs in relevant to their host genes. 
#' @param Circ CircRNACount file. A file of circRNA read count table. First three columns are circRNA coordinates, and followed by columns for circRNA read counts, each sample per column. 
#' @param Linear LinearCount file. A file of circRNA host gene expression count table. Same configuration as CircRNACount file.
#' @param (Optional) CircCoordinates BED format circRNA coordinates file. 
#' @param group A vector of group indicators.
#' @param alpha p value cut off. Defaul 0.05.
#' @param plotsig If 'TRUE', significantly host-independently regulated circRNAs will be ploted.
#' @examples
#' data(Circ)
#' data(Linear)
#' test=Circ.test(Circ=Circ,Linear=Linear,group=c(rep(1,6),rep(2,6),rep(3,6)))
#' View(test$sig.dat)
#' # Plot one of them
#' Circ.ratioplot(Circ,Linear,Coordinates,plotrow=rownames(test$sig.dat)[1],groupindicator1=c(rep('1days',6),rep('4days',6),rep('20days',6)),groupindicator2=rep(c(rep('Female',4),rep('Male',2)),3),lab_legend='Ages')
#' @export Circ.test

Circ.test <- function(Circ,Linear,CircCoordinates=None,group,alpha=0.05,plotsig=T){
  # Requre packge
  require(aod)
  
  # check whether the input matrix are correct
  if ( nrow(Circ)!=nrow(Linear) | ncol(Circ) != ncol(Linear)){
    stop('Circ data and Linear data are not matched, dimention different.')
  }
  
  # A vector for pvalue
  p.val <- c()
  
  # groups
  if ( length(group) != ncol(Circ)-3 ){
    stop("length of 'group' must be equal to the number of samples of 'Circ' and 'Linear'. ")
  }
  group <- factor(group)

  ## test
  # constract test matrix for each circRNA  
  for ( i in rownames(Circ) ){
    #print (i)    
    # total read counts vector
    tot <- round( as.numeric(Linear[i,c(4:ncol(Circ))]) + as.numeric(Circ[i,c(4:ncol(Linear))]) )

    # circRNA read counts
    circ <- as.numeric(Circ[i,c(4:ncol(Circ))])

    
    # if there is 0 in the total count vector, the model will fail. So permute 0 to 1
    if ( 0 %in% tot ){
      tot[tot==0]=1
    }
    
    # Constract data frame
    testdat = data.frame(tot,circ,group)
    
    ## do test
    # Null model
    fitNull <- betabin(cbind(circ,tot-circ) ~ 1, ~ 1, data=testdat)
    # Alternative model
    fitAlt <- betabin(cbind(circ,tot-circ) ~ group, ~ 1, data=testdat)
    # test models
    a <- anova(fitNull,fitAlt)
    p.value <- a@anova.table[,11][2]
    print(p.value)
    p.val <- c( p.val, p.value ) 
  }
  
  p.adj <- p.adjust(p.val,n=sum(!is.na(p.val)),'BH')
  # select significant ones
  sig_dat <- Circ[p.adj<=alpha  & !is.na(p.adj),]
  sig_p <- p.adj[p.adj<=alpha  & !is.na(p.adj)]
  # sort
  sig_dat <- sig_dat[order(sig_p),]
  sig_p <- sort(sig_p)
  # A summary table
  if (missing(CircCoordinates)){
    summary_table <- cbind(sig_dat[,c(1:3)],sig_p)
  }else{
    summary_table <- cbind(CircCoordinates[rownames(sig_dat),],sig_p)
  }
  return(list(summary_table=summary_table,sig.dat=sig_dat,p.val=p.val,p.adj=p.adj,sig_p=sig_p))
}