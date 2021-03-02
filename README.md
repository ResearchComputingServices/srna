<h1>Program to compute sRNAs for a given input genome in python

The algorithm to compute sRNAs is provided by Prof. Alex Wong at Carleton University. 

<h2>How to run the program in your local machine
  
<h3>Prerequisites

**Blast:**
1. If blast is not already installed, install blast in your local machine by following these instructions:

https://www.ncbi.nlm.nih.gov/books/NBK279671/

**Python:**

1. Check if python is installed on your local machine by opening a command window and typing: *python --version*
2. If python is installed, the output should include the version of python that is installed in your machine. For example: *Python 2.7.16*
3. If python is not installed, you can follow these directions to install it. 
https://www.python.org/downloads/


<h3> Run sRNA locally
  
 1. Donwload the latest release of the code.
 2. Go to the directory where the sRNA code is located.
 3. Create a virtual enviroment for python by typing: 
    *python3 -m venv env*
    or by following these instructions: https://docs.python.org/3/library/venv.html
 4. Activate your virtual enviroment by typing:
    *source ./env/bin/activate*
 5. The program to be run is *main.py* and the parameters that it receives are the following:





usage: main.py [-h] [-t TAGS] [-r RECOMPUTE]
               Sequence file Format sequence Shift position Length Expected cutoff Identity
               Percentage Cutoff

Computes sRNAs for a given genome

positional arguments:
  Sequence file         Sequence File that contains the genome (including absolute path)
  Format sequence       Format of the sequence file
  Shift position        Shift position to compute the sRNAs
  Length                sRNAs length
  Expected cutoff       Expected cutoff when blasting sRNAs agains input genome
  Identity Percentage Cutoff
                        Percentage of identity used when blasting sRNAs agains input genome
                        (a value between 0 and 1)

optional arguments:
  -h, --help            show this help message and exit
  -t TAGS, --tags TAGS  Excel file that includes the locus/gene tags to compute the sRNAS
  -r RECOMPUTE, --recompute RECOMPUTE
                        Shift position when recomputing sRNAS with hits





