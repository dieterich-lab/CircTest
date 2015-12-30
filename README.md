# CircTest: test the variation of circRNAs in respect to host genes.
# Installation

with `devtools`:

```S
> install.packages("devtools")
> require(devtools)
> install_github('dieterich-lab/CircTest')
> library(CircTest)
```
CircTest is strongly recommended to work with the output of DCC: https://github.com/dieterich-lab/DCC, but can also run on your own tables. 
All you need is one table with circular counts and one with host-gene counts. These tables should have the same order, i.e. circ[i,j] and host[i,j] are read-counts for the same circle in the same sample.


# Usage

Your table may have as many columns describing the circle or just one column containing the circle ID followed by as many columns of read counts.

Example tables:

```S
 > Circ
 #		CircID SRR1197279 SRR1197275 SRR1197273 SRR1197274
 # 1 2L:8377523|8377575         22         14         22         16
 # 2 2L:8377523|8377575          3          7          8          8
 # 3 2L:8377523|8377575          1          0          1          1
 # 4 2L:8377523|8377575          1          1          4          3
 # 5 2L:8377523|8377575          5          6          6          6
 # 6 2L:8377523|8377575         12          4          8          6
 
 > Linear
 #		CircID SRR1197279 SRR1197275 SRR1197273 SRR1197274
 # 1 2L:8377523|8377575        122        135        242        165
 # 2 2L:8377523|8377575         34         25         31         39
 # 3 2L:8377523|8377575         29         33         42         55
 # 4 2L:8377523|8377575         28         21         34         21
 # 5 2L:8377523|8377575        110        104        109        110
 # 6 2L:8377523|8377575         31         28         36         28
```

1) Read in tables

2) Filter tables 
To model expression data using the beta binomial distribution and testing for differences in groups it is benefical to only test well supported circles. You can use the package's function Circ.filter to filter your tables.

Nreplicates specifies the number of replicates in each condition.
filter.sample specifies the number of samples the circle has to have enough circular reads in to be considered.
filter.count specifies the circular read count threshold.
percentage specifies the minimum circle to host-gene ratio
circle_description tells the function which columns are NOT filled with read counts but the circle's annotation.

```S
  # filter circles by read counts
  Circ_filtered <- Circ.filter(circ = Circ, linear = Linear, Nreplicates = 2, filter.sample = 3, filter.count = 5, percentage = 0.1, circle_description = 1)
  # filter linear table by remaining circles
  Linear_filtered <- Linear[rownames(Circ_filtered),]
  > Circ_filtered
  #		      CircID Young1 Young2 Young3 Young4 Young5 Young6 Adult1 Adult2 Adult3 Adult4 Adult5 Adult6 Old1 Old2 Old3 Old4 Old5 Old6
  # 53    3L:7655571|7656424     37     32     49     31     54     51     77     84     97     51    107     99  316  293  233  246  462  507
  # 65   3LHet:209156|209624    132    116    134     84    140    104    248    151    179    134    173    251  370  336  274  288  418  464
  # 18    2R:6681630|6681989     16     16     15     11     20     20     22     25     28     21     34     28   41   47   35   25   63   75
  # 11      2R:494269|494825     39     39     40     47     61     52     94     56     51     46     66    111   96   74   73   99  153  169
  # 93       4:958958|959743     10      8     18     16     10     13     16     14     20     10     19     25   31   24   31   24   37   63
  # 68 3LHet:2454979|2455331     99     88     79     94     75     91    130     80    108     85    117    198  197  187  169  184  280  308

  > Linear_filtered
  # 	 	      CircID Young1 Young2 Young3 Young4 Young5 Young6 Adult1 Adult2 Adult3 Adult4 Adult5 Adult6 Old1 Old2 Old3 Old4 Old5 Old6
  # 53    3L:7655571|7656424   1020   1138   1340   1049   1247   1108   1268   1237    926    863   1204   1181 1352 1260  773  846 1783 1525
  # 65   3LHet:209156|209624    482    427    508    367    484    383    611    449    463    375    516    687  559  555  436  480  696  856
  # 18    2R:6681630|6681989    266    260    309    241    279    245    231    243    241    186    280    257  208  235  151  144  287  270
  # 11      2R:494269|494825    857    661    792    623    827    714    878    577    587    448    811   1123  553  551  475  569  923 1133
  # 93       4:958958|959743    284    203    267    212    271    232    259    182    191    142    280    321  175  144  123  120  240  290
  # 68 3LHet:2454979|2455331    577    479    557    417    469    459    524    461    474    378    542    682  592  575  481  504  802  912
```

3) Test for changes

Circ.test uses the beta binomial distribution to model the data and perform and ANOVA to identify circles which differ in their relative expression between the groups.
It is important that the grouping is correct (group) and the non-read-count columuns are specified (circle_description).

```S
  test=Circ.test(Circ_filtered, Linear_filtered, group=c(rep(1,6),rep(2,6),rep(3,6)), circle_description = 1)

  > test
  # $summary_table
  #		sig_p
  # [1,] 3 1.399547e-12
  # [2,] 4 1.399547e-12
  # [3,] 2 7.694809e-10
  # [4,] 1 1.025458e-09
  # [5,] 6 6.384966e-09
  # [6,] 5 5.141057e-08

  # $sig.dat
  #                  CircID Young1 Young2 Young3 Young4 Young5 Young6 Adult1 Adult2 Adult3 Adult4 Adult5 Adult6 Old1 Old2 Old3 Old4 Old5 Old6
  # 53    3L:7655571|7656424     37     32     49     31     54     51     77     84     97     51    107     99  316  293  233  246  462  507
  # 65   3LHet:209156|209624    132    116    134     84    140    104    248    151    179    134    173    251  370  336  274  288  418  464
  # 18    2R:6681630|6681989     16     16     15     11     20     20     22     25     28     21     34     28   41   47   35   25   63   75
  # 11      2R:494269|494825     39     39     40     47     61     52     94     56     51     46     66    111   96   74   73   99  153  169
  # 93       4:958958|959743     10      8     18     16     10     13     16     14     20     10     19     25   31   24   31   24   37   63
  # 68 3LHet:2454979|2455331     99     88     79     94     75     91    130     80    108     85    117    198  197  187  169  184  280  308

  # $p.val
  # [1] 3.229639e-13 4.665157e-13 3.847405e-10 6.836387e-10 5.320805e-09 5.141057e-08

  # $p.adj
  # [1] 1.399547e-12 1.399547e-12 7.694809e-10 1.025458e-09 6.384966e-09 5.141057e-08

  # $sig_p
  # [1] 1.399547e-12 1.399547e-12 7.694809e-10 1.025458e-09 6.384966e-09 5.141057e-08
```
4) Visualize data


```S
 for (i in rownames(test$summary_table))  {
  Circ.ratioplot( Circ_filtered, Linear_filtered, plotrow=i, groupindicator1=c(rep('1days',6),rep('4days',6),rep('20days',6)), lab_legend='Ages', circle_description = 1 )
 }
```



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
                  lab_legend='Ages' )
 }
```
