# CircTest: test the variation of circRNAs in respect to host genes.
# Installation

with `devtools`:

```S
> install.packages("devtools")
> require(devtools)
> install_github('dieterich-lab/CircTest')
```
# Usage


# Usage with DCC output

CircTest is strongly recommended to work with the output of DCC: https://github.com/dieterich-lab/DCC, but can also run on your own tables. 
All you need is one table with circular counts and one with host-gene counts. These tables should have the same order, i.e. circ[i,j] and host[i,j] are read-counts for the same circle in the same sample.
Your table may have as many columns describing the circle or just one column with a circle ID followed by as many columns with counts.

Example tables:

```S
 # circRNA

   CircID	SRR1197279	SRR1197275	SRR1197273	SRR1197274
 1	2L:227372|231034	22	14	22	16
 2	2L:2704926|2712136	3	7	8	8
 3	2L:5274019|5275780	1	0	1	1
 4	2L:5850831|5853258	1	1	4	3
 5	2L:7107141|7107551	5	6	6	6
 6	2L:8377523|8377575	12	4	8	6

 # linear

```



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
