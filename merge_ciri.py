# merges all ciri files in a given folder outputs three files junction_reads.txt non_junction_reads.txt and junction_reads_ratio.txt


import os
import argparse


parser = argparse.ArgumentParser(description='Finds all ciri files in a given folder and summarizes all circles over all samples. A summary file of circle counts, linear counts and junction ratios are exported into separate files to use for CircTest or similar programs.')
parser.add_argument('in_folder', metavar = 'input_folder', help = 'ciri output folder')
parser.add_argument('-s', dest = 'suffix', default = 'ciri', help = 'Suffix of all ciri files after splitting filename by ".". Only circle files not log files. (default: ciri)')
parser.add_argument('-t', dest = 'threshold', default = 2, type = int, help = 'Only report circle if it is present in a least <t> samples. (default: 2)')

args = parser.parse_args()
ciri = args.in_folder
suffix = args.suffix
threshold = args.threshold


def read_file(infile, circles):
    """
    Parses ciri file into several dictionaries
    
    Args:
	infile: ciri file to be parsed
	circles: pre-initialized dictionary. May already contain circles from other files or empty dictionary
    
    Returns:
	circles: dictionary with all circles regardless in which sample they occur
	circle_data: sample specific circles
    """
    # sample circle dictionary
    circle_data = {}
    I = open(infile)
    # ignore first line. This should only be the header. If the files were produced using ciri, this is true, if files were produced using something else, one would need to double check this
    I.readline()
    # iterate over all lines
    while True:
	Line = I.readline()
	if not Line:break
	# split line by column separator
	L = Line.replace(' ', '_').split('\t')
	# if a circle appears for the first time, add it to the circles dictionary along with its annotation
	if not L[0] in circles:
	    circles[L[0]] = [L[1],L[2], L[3], L[8], L[9]] # chr start end circRNA_type gene_id 
	# add sample specific counts to a seperate circle dictionary
	circle_data[L[0]] = [L[1],L[2], L[3], L[4], L[6], L[7], L[8], L[9]] # chr start end #junction_reads #non_junction_reads junction_reads_ratio circRNA_type gene_id 
    I.close()
    return(circles, circle_data)



def iterate_over_files(in_folder, suffix):
    """
    Iterates over all files in a given folder which have the given suffix
    
    Args:
	in_folder: path to the folder containing all ciri output files which should be merged
	suffix: suffix of all ciri circle files, this is neccessary to avoid problems with the log files typically reported by ciri.
    
    Returns: 
	circles: dictionary of all circles. key: circle_id, value: circle annotation
	sample_dict: dictionary of all samples and circles. key: sample_name, value: { circle_id : junction/nonjunction counts, circle ratio }
    """
    Files = os.listdir(in_folder)
    sample_dict = {}
    circles = {}
    # iterate over all files
    for lola in Files:
	# find ciri files
	if lola.split('.')[-1] == suffix:
	    # print file so that the user may know which files were actually used
	    print(lola)
	    # parse ciri files
	    circles, circle_dict = read_file('%s/%s'% (in_folder, lola), circles)
	    # add circles to sample dict
	    sample_dict['.'.join(lola.split('.')[:-1])] = circle_dict
    return(circles, sample_dict)

def write_circle_files(output_folder, circles, sample_dict, threshold):
    """
    Creates separate files for junction_reads, non_junction_reads and junction_ratio.
    Writes out all circles that occur in at least t samples into the tab-separated files.
    
    Args:
	output_folder: folder in which output files will be deposited, same as input folder
	circles: Dictionary containing all circles along with their annotation
	sample_dict: Dictionary of all samples along with the junction, non-junction reads counts and junction ratio
	threshold: circle has to be present in at least this many samples to be reported
    """
    # Generate sorted keys for circles and samples
    samples = sample_dict.keys()
    C_keys = circles.keys()
    SCK = sorted(C_keys)
    SS = sorted(samples)
    
    # open junction, non_junction, and ratio file
    O_jun_reads = open('%s/junction_reads.txt' %(output_folder), 'w')
    O_no_jun_reads = open('%s/non_junction_reads.txt' %(output_folder), 'w')
    O_ratio = open('%s/junction_reads_ratio.txt' %(output_folder), 'w') 
    
    # generate header including all samples separated by tabs
    header = 'circle\tchr\tstart\tend\tcircRNA_type\tgene_id\t%s\n'%('\t'.join(SS))
    O_jun_reads.write(header)
    O_no_jun_reads.write(header)
    O_ratio.write(header)
    
    # iterate over all circles
    for lola in SCK:
	# for each circle count the number of samples it appears in:
	sample_threshold = 0
	for forrest in SS:
	    if lola in sample_dict[forrest]:
		sample_threshold += 1
	
	# if circles passes the threshold, write out circle information
	if sample_threshold >= threshold:
	    circle_info = '%s\t%s' %(lola, '\t'.join(circles[lola]))
	    O_jun_reads.write(circle_info)
	    O_no_jun_reads.write(circle_info)
	    O_ratio.write(circle_info)
	    
	    # and for each sample write out read counts
	    for forrest in SS:
		if lola in sample_dict[forrest]:
		    O_jun_reads.write('\t%s' %(sample_dict[forrest][lola][3]))
		    O_no_jun_reads.write('\t%s' %(sample_dict[forrest][lola][4]))
		    O_ratio.write('\t%s' %(sample_dict[forrest][lola][5]))
		else:
		    O_jun_reads.write('\t0')
		    O_no_jun_reads.write('\t0')
		    O_ratio.write('\t0')
	    
	    # finish up the line
	    O_jun_reads.write('\n')
	    O_no_jun_reads.write('\n')
	    O_ratio.write('\n')
    
    # close all files
    O_jun_reads.close()
    O_no_jun_reads.close()
    O_ratio.close()
    return




C, S = iterate_over_files(ciri, suffix)
write_circle_files(ciri, C, S, threshold)


