library(CircTest)

Circ <- read.delim('Circ.csv', header = T, as.is = T)
Linear <- read.delim('Linear.csv', header = T, as.is = T)
CircCoord <- read.delim('CircCoordinates.csv', header = T, as.is = T)

Circ
#            CircID   Chr Start   End Control_1 Control_2 Control_3 Treatment_1 Treatment_2 Treatment_3
# 1    chr1:100|800  chr1   100   800         0         2         1           5           4           0
# 2 chr1:1050|10080  chr1  1050 10080        20        22        21          10          13           0
# 3   chr2:600|1000  chr2   600  1000         0         1         0          10           0           1
# 4 chr10:4100|5400 chr10  4100  5400        55        54        52          56          53          50
# 5  chr11:600|1500 chr11   600  1500         3         0         1           2           2           3

Linear
#            CircID   Chr Start   End Control_1 Control_2 Control_3 Treatment_1 Treatment_2 Treatment_3
# 1    chr1:100|800  chr1   100   800        10        11        12           9          10          10
# 2 chr1:1050|10080  chr1  1050 10080        80        81        83          45          48          46
# 3  chr2: 600|1000  chr2   600  1000         5         5         2          12           8           7
# 4 chr10:4100|5400 chr10  4100  5400       101       110       106         150         160         153
# 5  chr11:600|1500 chr11   600  1500        20        21        18          19          20          20

CircCoord
#            CircID   Chr Start   End   Gene
# 1    chr1:100|800  chr1   100   800 Gene_A
# 2 chr1:1050|10080  chr1  1050 10080 Gene_B
# 3  chr2: 600|1000  chr2   600  1000 Gene_C
# 4 chr10:4100|5400 chr10  4100  5400 Gene_D
# 5  chr11:600|1500 chr11   600  1500 Gene_E


Circ_filtered <- Circ.filter(circ = Circ, linear = Linear, Nreplicates = 3, filter.sample = 3, filter.count = 5, percentage = 0.1, circle_description = c(1:4))
Linear_filtered <- Linear[rownames(Circ_filtered),]
CircCoord_filtered <- CircCoord[rownames(Circ_filtered),]

Circ_filtered
#            CircID   Chr Start   End Control_1 Control_2 Control_3 Treatment_1 Treatment_2 Treatment_3
# 2 chr1:1050|10080  chr1  1050 10080        20        22        21          10          13           0
# 4 chr10:4100|5400 chr10  4100  5400        55        54        52          56          53          50

Linear_filtered
#            CircID   Chr Start   End Control_1 Control_2 Control_3 Treatment_1 Treatment_2 Treatment_3
# 2 chr1:1050|10080  chr1  1050 10080        80        81        83          45          48          46
# 4 chr10:4100|5400 chr10  4100  5400       101       110       106         150         160         153

CircCoord_filtered
#            CircID   Chr Start   End   Gene
# 2 chr1:1050|10080  chr1  1050 10080 Gene_B
# 4 chr10:4100|5400 chr10  4100  5400 Gene_D

test <- Circ.test(Circ_filtered, Linear_filtered, CircCoord_filtered, group=c(rep(1,3),rep(2,3)), circle_description = c(1:4))
# $summary_table
#            CircID   Chr Start  End   Gene      sig_p
# 4 chr10:4100|5400 chr10  4100 5400 Gene_D 0.01747407
# 
# $sig.dat
#            CircID   Chr Start  End Control_1 Control_2 Control_3 Treatment_1 Treatment_2 Treatment_3
# 4 chr10:4100|5400 chr10  4100 5400        55        54        52          56          53          50
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


pdf('example_pictures.pdf')
for (i in rownames(test$summary_table))  { 
  Circ.ratioplot(Circ_filtered, Linear_filtered, CircCoord_filtered, plotrow=i, groupindicator1=c(rep('Control',3),rep('Treatment',3)), 
		 lab_legend='Condition', circle_description = c(1:4), gene_column = 5 )
}

for (i in rownames(test$summary_table))  {
  Circ.lineplot(Circ_filtered, Linear_filtered, CircCoord_filtered, plotrow=i, groupindicator1=c(rep('Control',3),rep('Treatment',3)),
		circle_description = c(1:4), gene_column = 5 )
 }
dev.off()
      