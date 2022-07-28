import argparse
import sys
import os.path
from sRNA_provider import sRNA_Provider


my_parser = argparse.ArgumentParser(description='Computes sRNAs for a given genome')

my_parser.add_argument('sequenceFile',
                       metavar='sequence_file',
                       type=str,
                       help='Annotated genome file (including absolute path)')

my_parser.add_argument('formatSequence',
                       metavar='format_sequence',
                       type=str,
                       help='Format of the annotated genome file')

my_parser.add_argument('shift',
                       metavar='shift_position',
                       type=int,
                       help='Offset position to compute the asRNAs')

my_parser.add_argument('length',
                       metavar='length',
                       type=int,
                       help='asRNAs length')

my_parser.add_argument('e_cutoff',
                       metavar='expected_cutoff',
                       type=float,
                       help='Expected cutoff when blasting sRNAs agains input genome')

my_parser.add_argument('identity_perc_cutoff',
                       metavar='Identity Percentage Cutoff',
                       type=float,
                       help='Percentage of identity used when blasting sRNAs agains input genome (a value between 0 and 1)')

my_parser.add_argument('-t',
                       '--tags',
                       action='store',
                       help='Excel file that includes the locus/gene tags to compute the sRNAS')

my_parser.add_argument('-r',
                       '--recompute',
                       action='store',
                       help='Shift position when recomputing sRNAS with hits')

my_parser.add_argument('-c',
                       '--csv',
                       action='store_true',
                       help='Save output to CSV format files (default Excel format only)')

my_parser.add_argument('-e',
                       '--excel',
                       action='store_true',
                       help='Save output to Excel format files (default Excel format only)')

# Execute the parse_args() method
args = my_parser.parse_args()

file_sequence_fullpath = args.sequenceFile
format = args.formatSequence
position = args.shift
length = args.length
e_cutoff = args.e_cutoff
identity_perc_cutoff = args.identity_perc_cutoff
file_tags = args.tags
recompute_position_st = args.recompute
excel_output = args.excel or not args.csv
csv_output = args.csv

if file_sequence_fullpath:
    split = os.path.split(file_sequence_fullpath)
    base_directory = split[0]
    seq_file = split[1]

if not os.path.isfile(file_sequence_fullpath):
    print('The specified input sequence file does not exist.')
    sys.exit()

if int(position)==0:
    print ('Shift position cannot be 0. It is either >0 or <0.')
    sys.exit(2)

if file_tags and not os.path.isfile(file_tags):
    print('The specified excel tags file does not exist.')
    sys.exit()

recompute_position = None
if recompute_position_st:
    try:
        recompute_position = int(recompute_position_st)
    except:
        print ('Recompute position must be an integer.')
        sys.exit()

    if recompute_position == 0:
        print('Shift position for recomputing cannot be 0. It is either >0 or <0.')
        sys.exit()
    else:
        if recompute_position == position:
            print('Shift positions cannot be the same.')
            sys.exit()


print ('Computing sRNAS with the following paramenters')
print('Input file:', file_sequence_fullpath)
print('Format: ', format)
print('Position for sRNA computation: ', position)
print('Length sRNA: ', length)
print('Expected cutoff (blast): ', e_cutoff)
print('Percentage Identity cutoff (blast): ', identity_perc_cutoff)
if file_tags:
    print('File tags: ', file_tags)
if recompute_position_st:
    print ('Position for recomputing: ', recompute_position)

sRNA_provider = sRNA_Provider()
sRNA_provider.compute_srnas(base_directory, seq_file, file_sequence_fullpath, format, position, length, e_cutoff, identity_perc_cutoff, file_tags, recompute_position, excel_output, csv_output)





