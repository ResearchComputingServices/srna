<h1>Program to compute sRNAs for a given input genome in python<h1>

The algorithm to compute sRNAs is provided by Prof. Alex Wong at Carleton University. 

<h2>How to run the program in your local machine <h2>
  
<h3>Prerequisites<h3>
 
**Blast:**
1. If blast is not already installed, install blast in your local machine by following these instructions:

https://www.ncbi.nlm.nih.gov/books/NBK279671/

**Python:**

1. Check if python is installed on your local machine by opening a command window and typing: *python --version*
2. If python is installed, the output should include the version of python that is installed in your machine. For example: *Python 2.7.16*
3. If python is not installed, you can follow these directions to install it. 
https://www.python.org/downloads/


<h3> Run sRNA locally <h3>
  
 1. Donwload the latest release of the code.
 2. Go to the directory where the sRNA code is located.
 3. Create a virtual enviroment for python by typing: 
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





