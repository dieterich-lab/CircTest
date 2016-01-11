# CircTest: test the variation of circRNAs in respect to host genes.
# Installation

with `devtools`:

```S
> install.packages("devtools")
> require(devtools)
> install_github('dieterich-lab/CircTest')
> library(CircTest)
```
CircTest is recommended to work with the output of DCC: https://github.com/dieterich-lab/DCC, but can also run on your own tables. 
All you need is one table with circular counts and one with host-gene counts. 
These tables should have the same order, i.e. circ[i,j] and linear[i,j] are read-counts for the same circle in the same sample.


# Usage

Your table may have many columns describing the circle or just one column containing the circle ID followed by many columns of read counts.

Example tables:  
  
Count table for back-spliced reads:
```S

           CircID Control_1 Control_2 Control_3 Treatment_1 Treatment_2 Treatment_3
     chr1:100|800         0         2         1           5           4           0
  chr1:1050|10080        20        22        21          10          13           0
   chr2: 600|1000         0         1         0          10           0           1
  chr10:4100|5400        55        54        52          56          53          50
   chr11:600|1500         3         0         1           2           2           3
 ```
   
Count table for host-gene reads:
```S
          CircID Control_1 Control_2 Control_3 Treatment_1 Treatment_2 Treatment_3
    chr1:100|800        10        11        12           9          10          10
 chr1:1050|10080        80        81        83          45          48          46
  chr2: 600|1000         5         5         2          12           8           7
 chr10:4100|5400       101       110       106         150         160         153
  chr11:600|1500        20        21        18          19          20          20

```

1) Read in tables

 ```S
  Circ <- read.delim('Circ.csv', header = T, as.is = T)
  Linear <- read.delim('Linear.csv', header = T, as.is = T)
```

2) Filter tables 
To model expression data using the beta binomial distribution and testing for differences in groups, it is benefical to only test well supported circles. You can use the package's function **Circ.filter** to filter your tables.

**Nreplicates** specifies the number of replicates in each condition.  
**filter.sample** specifies the number of samples the circle has to have enough circular reads in to be considered.  
**filter.count** specifies the circular read count threshold.  
**percentage** specifies the minimum circle to host-gene ratio.  
**circle_description** tells the function which columns are NOT filled with read counts but the circle's annotation.  

```S
  # filter circles by read counts
  Circ_filtered <- Circ.filter(circ = Circ, linear = Linear, Nreplicates = 3, filter.sample = 3, filter.count = 5, percentage = 0.1, circle_description = 1)
  #            CircID Control_1 Control_2 Control_3 Treatment_1 Treatment_2 Treatment_3
  # 2 chr1:1050|10080        20        22        21          10          13           0
  # 4 chr10:4100|5400        55        54        52          56          53          50

  
  # filter linear table by remaining circles
  Linear_filtered <- Linear[rownames(Circ_filtered),]
  #            CircID Control_1 Control_2 Control_3 Treatment_1 Treatment_2 Treatment_3
  # 2 chr1:1050|10080        80        81        83          45          48          46
  # 4 chr10:4100|5400       101       110       106         150         160         153

  ```

3) Test for changes

**Circ.test** uses the beta binomial distribution to model the data and performs an ANOVA to identify circles which differ in their relative expression between the groups.  
It is important that the grouping is correct (**group**) and the non-read-count columuns are specified (**circle_description**).

```S
  test <- Circ.test(Circ_filtered, Linear_filtered, group=c(rep(1,3),rep(2,3)), circle_description = 1)
  # $summary_table
  #            CircID      sig_p
  # 4 chr10:4100|5400 0.01747407
  # 
  # $sig.dat
  #            CircID Control_1 Control_2 Control_3 Treatment_1 Treatment_2 Treatment_3
  # 4 chr10:4100|5400        55        54        52          56          53          50
  # 
  # $p.val
  # [1] 0.153464107 0.008737037
  # 
  # $p.adj
  # [1] 0.15346411 0.01747407
  # 
  # $sig_p
  # [1] 0.01747407
  # 
  ```
4) Visualize data

The CircTest library has build in plot functions to view your significantly different genes.  
To visualize the ratio as barplot try:
```S
  for (i in rownames(test$summary_table))  { 
    Circ.ratioplot(Circ_filtered, Linear_filtered, plotrow=i, groupindicator1=c(rep('Control',3),rep('Treatment',3)), 
		   lab_legend='Condition', circle_description = 1 )
  }
```

To visualize the abundance of host-gene and circle separately in a line plot try:
```S
  for (i in rownames(test$summary_table))  {
    Circ.lineplot(Circ_filtered, Linear_filtered, plotrow=i, groupindicator1=c(rep('Control',3),rep('Treatment',3)),
		  circle_description = 1 )
 }
```

More examples on how to deal with different data can be found in the **example_data** folder.  

# Usage with DCC output



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
                  lab_legend='Ages', gene_column = 4 )
 }
```

# Plugin for CIRI data

If you have circle files generated by CIRI, you may but all sample files into one folder and run the python script: 
```S
python merge_ciri.py /path/to/all/ciri/files
```
This will procude three files. Please use the junction_reads.txt as Circ, the non_junction_reads.txt as Linear, and change **circle_description = c(1:6)**
