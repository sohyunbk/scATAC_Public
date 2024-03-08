import pysam
import argparse
import ast
import os


## From pablo # Sohyun edited for cellranger v2 and do not consider non nuclear configs

def parse_config_file(file_name):
    """
    Reads in the 10X config file used to generate the genome, parses this file
    for later coutning of the proportion of Nuclear and non-nuclear reads by
    bardcodes.
    """
    #Generate nested list. First list will be Nucelar contigs, second list will
    #be non nuclear contigs
    nuclear_non_nuclear_list = []
    try:
        with open(file_name, 'r') as f: ##Open .fai file
            list =[]
            for line in f:
                Chr = line.split("\t")[0]
                list.append(Chr)
            nuclear_non_nuclear_list.append(list)
        return(nuclear_non_nuclear_list)
    except OSError:
        print("File doest not exits %s, Exiting" % file_name)
        sys.exit(-1)

def read_bam_file(nuclear_non_nuclear_list, exp_name, bam_file):
    """
    First alter the CB tags as well as the other tags
    Next - count the tag and add tag to dictionary if not presenat and start at
    one. Add to total column, add to nuclear/nonnuclear column depending on the
    scaffold name and the list given above
    """
    dictionary_of_CB_counts = {}
    exp_name_tag = "-" + exp_name

    #Read the File
    #save = pysam.set_verbosity(0)
    #read_bam_file = pysam.AlignmentFile(bam_file,"rb", ignore_truncation=True)
    #pysam.set_verbosity(save)
    read_bam_file = pysam.AlignmentFile(bam_file,"rb")

    outfile = pysam.AlignmentFile("-", "w", template=read_bam_file)
    for read in read_bam_file :
        #For each read alter CB tag
        try:
            original_tag = (read.get_tag("CB"))
            #print("here")
            #print(original_tag)
            new_tag = original_tag.replace("-1", exp_name_tag)
            read.set_tag("CB", new_tag, replace=True)
            scaffold_name = (read.reference_name)

            final_tag = "CB:Z:" + read.get_tag("CB")

            if final_tag not in dictionary_of_CB_counts:
                if scaffold_name in nuclear_non_nuclear_list[0]:
                    dictionary_of_CB_counts[final_tag] = [1,1]
                #elif scaffold_name in nuclear_non_nuclear_list[1]:
                    #dictionary_of_CB_counts[final_tag] = [1,0]

            elif final_tag in dictionary_of_CB_counts:
                if scaffold_name in nuclear_non_nuclear_list[0]:
                    dictionary_of_CB_counts[final_tag][0] += 1
                    dictionary_of_CB_counts[final_tag][1] += 1
                #elif scaffold_name in nuclear_non_nuclear_list[1]:
                    #dictionary_of_CB_counts[final_tag][0] += 1
            outfile.write(read)
        except KeyError:
            pass
    return(dictionary_of_CB_counts)

def write_CB_file(exp_name, tissue, dict_to_write, output):
    try:
        os.remove(output)
    except OSError:
        columns = ["cellID", "total", "nuclear_library","tissue"]
        with open(output, 'a') as f:
            f.write('\t'.join(columns))
            f.write('\n')
            for key, value in dict_to_write.items():
                list_gen = [key, str(value[0]), str(value[1]), tissue]
                f.write('\t'.join(list_gen))
                f.write('\n')



def get_parser():
    parser = argparse.ArgumentParser(description='Pull our reads aligning to a\
        region from multiple list of BAM files, puts them into a BAM file\
        for later assembly.')
    parser.add_argument('-10x_config', "--10x_config", help="10x config file to \
            pull scaffold names of nulcear and non-nuclear scaffolds",
            required=True, dest='config')
    parser.add_argument('-exp_name', "--experiment_name", help="10x config file to \
            pull scaffold names of nulcear and non-nuclear scaffolds",
            required=True, dest='exp')
    parser.add_argument('-tissue', "--tissue", help="10x config file to \
        pull scaffold names of nulcear and non-nuclear scaffolds",
        required=True, dest='tis')
    parser.add_argument('-BAM', "--bam_file", help="Bam file to \
        pull reads from.", required=True, dest='bam_f')
    parser.add_argument('-output', "--output_count", help="Bam file to \
        pull reads from.", required=True, dest='out')


    args = vars(parser.parse_args())
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()

    #Load all Bed files
    nuclear_and_non_nuclear_scaffolds = parse_config_file(args.config)

    gathered_read_dict = read_bam_file(nuclear_and_non_nuclear_scaffolds, args.exp, args.bam_f)
    write_CB_file(args.exp, args.tis, gathered_read_dict, args.out)
