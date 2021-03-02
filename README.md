**Program to compute sRNAs for a given input genome in python**

**How to run the program in your local machine**
  
**Prerequisites**
 
- Blast:
If blast is not already installed, install blast in your local machine by following these instructions:

https://www.ncbi.nlm.nih.gov/books/NBK279671/

- Python:

1. Check if python is installed on your local machine by opening a command window and typing: *python --version*
2. If python is installed, the output should include the version of python that is installed in your machine. For example: *Python 2.7.16*
3. If python is not installed, you can follow these directions to install it. 
https://www.python.org/downloads/


**Run sRNA locally**
  
 1. Donwload the latest release of the code.
 2. Go to the directory where the sRNA code (srna directory) is located.
 4. Create a virtual enviroment for python by typing: 
    *python3 -m venv env*
    or by following these instructions: https://docs.python.org/3/library/venv.html
 4. Activate your virtual enviroment by typing:
    *source ./env/bin/activate*
 5. The program to be run is *main.py* and the parameters that it receives are the following:

Here is the order in which the program receives the parameters:

*python main.py  sequence_file format_sequence shift_position length expected_cutoff identity_percentage_cutoff [-t TAGS] [-r RECOMPUTE]*

Position Arguments | Meaning 
-------------------|---------
  sequence_file         | Sequence File that contains the genome (including absolute path)
  format_sequence       | Format of the sequence file
  shift_position        | Shift position to compute the sRNAs
  length                | sRNAs length
  expected_cutoff       | Expected cutoff when blasting sRNAs agains input genome
  identity_percentage_cutoff |Percentage of identity used when blasting sRNAs agains input genome (a value between 0 and 1)

Optional arguments | Meaning
------------------ | -------
  -h, --help       |    Shows help about how the program usage
  -t TAGS, --tags TAGS  | Excel file that includes the locus/gene tags to compute the sRNAS
  -r RECOMPUTE, --recompute RECOMPUTE | Shift position when recomputing sRNAS with hits

**Examples:

The directory *sequences* in the directory *srna* contains some sample genome sequences: K12.gb and JWGZ01.1.gbff. It also includes some samples of locus tags for these sequences.

- Example 1:
Suppose that you would like to compute all sRNAs for the sequence K12.gb. The position to compute the sRNAS would be -10 and the length of the sRNAS would be 19. Additionally, for the sRNAs that contain hits in the genome, the program should recompute these sRNAS with a position of -15. Assume that the program was installed in the following path: /home/srna. Therefore, the sequences are located at directory /home/srna/sequences. For this example, the program should be executed like this:








