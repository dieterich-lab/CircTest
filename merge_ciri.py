# merges all ciri files in a given folder outputs three files junction_reads.txt non_junction_reads.txt and junction_reads_ratio.txt


import os
import argparse


parser = argparse.ArgumentParser(description='')
parser.add_argument('in_folder', metavar = 'input', help = 'ciri output folder')

args = parser.parse_args()
ciri = args.in_folder

def read_file(infile, circles):
    circle_data = {}
    I = open(infile)
    I.readline()
    while True:
	Line = I.readline()
	if not Line:break
	L = Line.replace(' ', '_').split('\t')
	#if float(L[7]) < 1.0: # and float(L[7]) >= 0.05 :
	if not L[0] in circles:
	    circles[L[0]] = [L[1],L[2], L[3], L[8], L[9]] # chr start end circRNA_type gene_id 
	circle_data[L[0]] = [L[1],L[2], L[3], L[4], L[6], L[7], L[8], L[9]] # chr start end #junction_reads #non_junction_reads junction_reads_ratio circRNA_type gene_id 
    I.close()
    return(circles, circle_data)

def iterate_over_files(in_folder):
    Files = os.listdir(in_folder)
    sample_dict = {}
    circles = {}
    for lola in Files:
	if lola.split('.')[-1] == 'ciri':
	    print(lola)
	    circles, circle_dict = read_file('%s/%s'% (in_folder, lola), circles)
	    sample_dict['.'.join(lola.split('.')[:-1])] = circle_dict
	    print(len(circle_dict))
    print(len(circles))
    return(circles, sample_dict)

def write_circle_files(output_folder, circles, sample_dict):
    samples = sample_dict.keys()
    C_keys = circles.keys()
    SCK = sorted(C_keys)
    SS = sorted(samples)
    O_jun_reads = open('%s/junction_reads.txt' %(output_folder), 'w')
    O_no_jun_reads = open('%s/non_junction_reads.txt' %(output_folder), 'w')
    O_ratio = open('%s/junction_reads_ratio.txt' %(output_folder), 'w') 
    header = 'circle\tchr\tstart\tend\tcircRNA_type\tgene_id\t%s\n'%('\t'.join(SS))
    O_jun_reads.write(header)
    O_no_jun_reads.write(header)
    O_ratio.write(header)
    for lola in SCK:
	circle_info = '%s\t%s' %(lola, '\t'.join(circles[lola]))
	O_jun_reads.write(circle_info)
	O_no_jun_reads.write(circle_info)
	O_ratio.write(circle_info)
	for forrest in SS:
	    if lola in sample_dict[forrest]:
		O_jun_reads.write('\t%s' %(sample_dict[forrest][lola][3]))
		O_no_jun_reads.write('\t%s' %(sample_dict[forrest][lola][4]))
		O_ratio.write('\t%s' %(sample_dict[forrest][lola][5]))
	    else:
		O_jun_reads.write('\t0')
		O_no_jun_reads.write('\t0')
		O_ratio.write('\t0')
	O_jun_reads.write('\n')
	O_no_jun_reads.write('\n')
	O_ratio.write('\n')
    O_jun_reads.close()
    O_no_jun_reads.close()
    O_ratio.close()
    return




C, S = iterate_over_files(ciri)
write_circle_files(ciri, C, S)


