#' @title Circ.test
#'
#' @description
#' Test the independent variation of circRNAs in relevant to their host genes.
#' @param Circ CircRNACount file. A file of circRNA read count table. First three columns are circRNA coordinates, and followed by columns for circRNA read counts, each sample per column.
#' @param Linear LinearCount file. A file of circRNA host gene expression count table. Same configuration as CircRNACount file.
#' @param (Optional) CircCoordinates BED format circRNA coordinates file.
#' @param group A vector of group indicators.
#' @param alpha p value cut off (FDR). Default 0.05.
#' @param plotsig If 'TRUE', significantly host-independently regulated circRNAs will be ploted.
#' @param circle_description Column indices which do not carry circle/linear read counts.
#' @examples
#' data(Circ)
#' data(Linear)
#' test=Circ.test(Circ=Circ,Linear=Linear,group=c(rep(1,6),rep(2,6),rep(3,6)))
#' View(test$sig.dat)
#' # Plot one of them
#' Circ.ratioplot(Circ,Linear,Coordinates,plotrow=rownames(test$sig.dat)[1],groupindicator1=c(rep('1days',6),rep('4days',6),rep('20days',6)),groupindicator2=rep(c(rep('Female',4),rep('Male',2)),3),lab_legend='Ages', circle_description = c(1:3))
#' @export Circ.test

Circ.test <- function(Circ, Linear, CircCoordinates=None, group, alpha=0.05, plotsig=T, circle_description = c(1:3)){

    # Requre packge
    require(aod)

    # check whether the input matrix are correct
    if ( nrow(Circ)!=nrow(Linear) | ncol(Circ) != ncol(Linear)){
        stop('Circ data and Linear data are not matched, dimention different.')
    }

    # A vector for pvalue and directions indicator
    p.val <- c()
    direction <- c()

    # groups
    if ( length(group) != ncol(Circ)-length(circle_description) ){
        stop("length of 'group' must be equal to the number of samples of 'Circ' and 'Linear'. ")
    }
    group <- factor(group)
    counter <- 0

    ## test
    # construct test matrix for each circRNA

    tmp_df = Circ[,FALSE]

    for (j in seq(1,length(unique(group)))){
        tmp_df[paste("group_",j,"_ratio_mean",sep="")] <- NA
    }

    for ( i in rownames(Circ) ){
        counter <- counter+1

        # total read counts vector
        tot <- round( as.numeric(Linear[i,-circle_description]) + as.numeric(Circ[i,-circle_description]) )

        # circRNA read counts
        circ <- as.numeric(Circ[i,-circle_description])

        # if there is 0 in the total count vector, the model will fail. So permute 0 to 1
        if ( 0 %in% tot ){
          tot[tot==0]=1
        }

        if (counter %% 1000 == 0){
            message(paste(counter, "candidates processed"))
        }

        tmp_rations <- data.frame(Ratio=as.numeric(Circ[i,-circle_description])/(as.numeric(Linear[i,-circle_description])+as.numeric(Circ[i,-circle_description])),
        group=group)
        for (rep_group in seq(1,max(as.numeric(levels(group))),1)){
            tmp_df[i, paste("group_",rep_group,"_ratio_mean",sep="")] <- mean(na.omit(unlist(tmp_rations[tmp_rations$group==rep_group,1])))
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
        # print(predict(fitAlt,testdat, se.fit=T))
        p.val <- c( p.val, p.value )
        # dir <- 1 # fitAlt@param[2][["group2"]]
        # direction <- c(direction, dir)
    }
    message(paste(counter, "candidates processed in total"))

    Circ$direction <- direction
    #names(Circ$direction ) <- c("direction")
    p.adj <- p.adjust(p.val,n=sum(!is.na(p.val)),'BH')
    # select significant ones
    sig_dat <- Circ[p.adj<=alpha  & !is.na(p.adj),]
    sig_ratios <- tmp_df[p.adj<=alpha  & !is.na(p.adj),]
    sig_p <- p.adj[p.adj<=alpha  & !is.na(p.adj)]
    # direction <- direction[p.adj<=alpha  & !is.na(p.adj)]

    # sort by p-val
    sig_dat <- sig_dat[order(sig_p),]
    sig_ratios <- sig_ratios[order(sig_p),]
    sig_p <- sort(sig_p)

    # A summary table
    if (missing(CircCoordinates)){
        summary_table <- data.frame(sig_dat[,circle_description],sig_p,sig_dat[,circle_description])

        rownames(summary_table) <- rownames(sig_dat)
        names(summary_table) <- c(names(sig_dat)[circle_description],"sig_p",names(sig_ratios)[circle_description])
    } else {
        # summary_table <- cbind(CircCoordinates[rownames(sig_dat),],sig_p,sig_dat$direction)
        # colnames(summary_table) <- c(colnames(CircCoordinates),"sig_p","direction")

        summary_table <- cbind(CircCoordinates[rownames(sig_dat),],sig_p,sig_ratios)
        colnames(summary_table) <- c(colnames(CircCoordinates),"sig_p",colnames(sig_ratios))
    }

    message(paste(nrow(summary_table), "candidates passed the specified thresholds"))

    # return all objects in a list
    return(list(summary_table=summary_table,
              sig.dat=sig_dat,
              p.val=p.val,
              p.adj=p.adj,
              sig_p=sig_p,
              ratios=sig_ratios
              # direction=direction
            )
        )
}
