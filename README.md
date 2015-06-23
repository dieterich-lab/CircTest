# CircTest: test the variation of circRNAs in respect to host genes.
# Installation

with `devtools`:

```S
> install.packages("devtools")
> require(devtools)
> install_github('dieterich-lab/CircTest',username='username')
```

# Usage

CircTest is strongly recommended to work with the output of DCC: https://github.com/dieterich-lab/DCC, otherwise, you need to organise your data in the same format as DCC outputs.


1) Read and load DCC output into R

```S

  library(CircTest)

  CircRNACount <- read.delim('CircRNACount',header=T)
  LinearCount <- read.delim('LinearCount',header=T)
  CircCoordinates <- read.delim('CircCoordinates',header=T)

  CircRNACount_filtered <- Circ.filter(circ = CircRNACount, linear = LinearCount, Nreplicates = 6, filter.sample = 6, filter.count = 5, percentage = 0.1)
  CircCoordinates_filtered <- CircCoordinates[rownames(CircRNACount_filtered),]
  LinearCount_filtered <- LinearCount[rownames(CircRNACount_filtered),]
```

Alternatively, load the processed Westholm et al. data from CircTest package.

```S
  
  library(CircTest)
  
  data(Circ)
  CircRNACount_filtered <- Circ
  data(Coordinates)
  CircCoordinates_filtered <- Coordinates
  data(Linear)
  LinearCount_filtered <- Linear
```

2) Test for host-independently regulated circRNAs

```S

 test=Circ.test(CircRNACount_filtered,LinearCount_filtered,CircCoordinates_filtered,group=c(rep(1,6),rep(2,6),rep(3,6)))
 # Significant result show in a summary table
 View(test$summary_table)
```

3) Visuallize the significantly host-independently regulated circRNAs

```S

 for (i in rownames(test$summary_table))  {
  Circ.ratioplot( CircRNACount_filtered, LinearCount_filtered, CircCoordinates_filtered, plotrow=i, 
                  groupindicator1=c(rep('1days',6),rep('4days',6),rep('20days',6)), 
                  lab_legend='Ages' )
 }
```
