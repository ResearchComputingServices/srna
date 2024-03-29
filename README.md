
### Algorithm for computing antisense RNAs

CAREn is described in Romero et al. (submitted). It generates antisense RNAs according to the following procedure.

Along with other parameters, the algorithm receives as input a genome sequence in a a GenBank or EMBL formatted file. These file contain two parts:  (1) annotations which give the locations of different genomic features (such as genes, CDS, or tRNAs), and (2) the DNA sequence itself.

The algorithm will compute asRNAs from the CDS features. For simplicity, there are two varieties of CDS:  (1) the CDS whose position is written like: CDS *i..m*. This indicates that the gene in this CDS starts at the position *i* and ends at position *m*. We will refer to this CDS as "encoded forward".

The second type of CDS are those whose position information includes the word "complement" CDS  *complement(i..m)*.  The word "complement" means that the gene is encoded backwards. For complement CDS, *m* would be the start of the gene and *i* is the end of the gene. 

Given a CDS *i…m*, for a gene encoded forward, the steps to compute an antisense RNA are as follows:
 
1. Get the position of the start of the gene. This will be the first number *i*. 
2. For a given shift position *p*,  
  - If *p<0*,  go *p* positions immediately before the start of the gene without including the *i* position. Let us denote this position as  *s = i + p*.   
  - If *p>0*,  go *p* positions at the right of the start of the gene including the start of the gene. That is, *s =i +p-1*.
 **Notice that there is no zero position p=0**. For a value of  *p=-1*, *s* would be the position immediately before the start of the gene (*s=i-1*), and for *p=1*, *s* would be the first position in the start of the gene (*s=i*).

3. For a given parameter *l*, capture a *l* sequence with the *s* position in the middle. Let *h = (l-1)/2*. The sequence would be captured at positions: *s-h…s+h*. 

4. Find the reverse complement to this sequence. This could be done by replacing every A with T, G with C, T with A, and C with G, then reverse it.  **This is the antisense sequence for the gene**.

Given a complement *CDS(i..m)* for a gene encoded backward, the steps to compute the asRNA are as follows:

1. Get the position of the start of the gene. This will be the second number: *m* 
2. For a given shift position p, 
  - If *p<0*,  go *p* positions before the start position without including the position *m*. Notice that since the gene is encoded backward, going before the start position corresponds to going to the right of the start position. Also notice that since *p<0*, for going to the right of the sequence, we should subtract *p*. That is: *s = m – p*.   
  - If *p>0*, go *p* positions after the start position including the start position. This would be equivalent to go *p* bases (letters) to the left of *m* including *m*, i.e., *s = m – p +1*.
Once again, there is no zero position *p=0*. For *p=-1*, *s* would be the position immediately before the start of the gene (*s=m+1*), and for *p=1*, *s=m* would be the first position in the start of the gene.

3. For a given parameter *l*, capture a *l* sequence with the *s* position in the middle. Let *h = (l-1)/2*. The sequence would be captured at the positions: *s-m…s+m*. **This is the antisense sequence for the gene**.
 

### Algorithm for computing a set of antisense RNAs from a genome

We describe here the general steps behind CAREn. For a given input genome, CAREn can compute either the asRNAs for **all** CDS in the genome, or for **a subset** of these CDS. The user can specify for which CDS the antisense RNAs (asRNAs) will be computed through a file of tags. Each tag  will correspond to a gene tag or a locus tag. If the tag is present in the genome, then an asRNA will be obtained for such a tag as described in the previous section.

After the set of asRNAs is obtained, the next step will be to look for similar sequences to these RNAs in the genome. This is to try to ensure that an asRNA only occurs in the gene or CDS where the asRNA was obtained from. The similarity information will be obtained using BLAST [1].   
For each asRNA, BLAST will output (among other information) a set of subsequences with an E value and a percentage of identity. Let us call these subsequences “hits”. For each pair (asRNA, hit), there is an E value and a percentage of identity returned by BLAST.  The E value describes how many times you would “expect” a match (asRNA, hit) to occur in the genome by chance, the closer to zero the value of E, the less likely is that the match (asRNA ,hit) will occur more than once. The percentage of identity describes how similar is the hit to the asRNA sequence. The higher the value of the percentage, the more similar are the hit and the asRNA [1,2].  The user will specify the E value (expected_cutoff) and the percentage of identity (identity_percentage_cutoff) for which BLAST will do a cutoff. That is, BLAST will return only the hits that pass those thresholds. 

Let us define as *R = [R1, R2,…,Rn]* the set of asRNAs which have more than one hit with *E<=expected_cutoff* and *P>=identity_percentage_cutoff*. As an option, the user could specify if a second computation would take place. That is, if for each asRNA *Ri* in *R*, the program will identify from which CDS *Ri* was obtained from and the program will re-obtain a new asRNA for that CDS at a different shift position *p’*, where *p<>p’*. After that, the program will BLAST the new set of reobtained asRNAs. 
 
The following summarizes CAREn’s operation:

**Input:** 
A genome sequence *S*,  
a set of gene/locus tags *T* (optional),
a shift position  *p*,
a length value *l*, 
a cutoff expected value *expected_cutoff*,
a cutoff percentage of identity *identity_percentage_cutoff*, 
a second shift position *p’* *p’<>p* (optional),

**Output:**  a set of asRNAs R and a set of asRNAs R' (optional)

1.	If *T* is not provided (T=none), make T=all CDS in the genome S.
2.	For each CDS in *T*, compute the asRNAs following the steps as described in Section 1 using parameters *p* and *l*.
3.	For each computed asRNA, run BLAST with parameters *expected_cutoff*, and *percentage_of_identity*.
4.	Let *R=[R1,…,Rn]* be the subset of asRNAs which have more than one hit in the genome at the expected value and percentage of identity cutoffs.
5.	If *p’* is none, return R.
6.	Otherwise, for each asRNA *Ri* in *R*, determine from which CDS the *Ri* was obtained from and compute a new asRNA as described in Section 1 using *p’* and *l*. Make this *R'* this new set of asRNAs.
7.	Return *R'*

## Program to compute asRNAs for a given input genome in python

### How to run the program in your local machine
  
**Prerequisites**
 
- Blast:
If blast is not already installed, install blast in your local machine by following these instructions:

https://www.ncbi.nlm.nih.gov/books/NBK279671/

- Python:

1. Check if python is installed on your local machine by opening a command window and typing: *python --version*
2. If python is installed, the output should include the version of python that is installed in your machine. For example: *Python 2.7.16*
3. If python is not installed, you can follow these directions to install it. 
https://www.python.org/downloads/


**Run CAREn locally**
  
 1. Download the latest release of the code.
 2. Go to the directory where the CAREn code (srna directory) is located.
 4. Create a virtual environment for python by typing: 
    *python3 -m venv env*
    or by following these instructions: https://docs.python.org/3/library/venv.html
 5. Activate your virtual environment by typing:
    *source ./env/bin/activate*
 6. Install the requirements for the environment:
    *pip install -r requirements.txt*
 7. The program to be run is *main.py* and the parameters that it receives are below.
 
**Note: Steps 1, 4 and 6 will be executed the first time the CAREn will be run. Afterwards, only execute steps 2, 5 and 7.**

Here is the order in which the program receives the parameters:

*python main.py  sequence_file   format_sequence    shift_position    length    expected_cutoff    identity_percentage_cutoff    [-t TAGS]   [-r RECOMPUTE]*

- Position Arguments:
  - sequence_file: Sequence file that contains the genome (including absolute path)
  - format_sequence: Format of the sequence file (e.g., genbank or embl)
  - shift_position: Shift position to compute the asRNAs
  - length: asRNAs length
  - expected_cutoff: Expected cutoff when blasting asRNAs against input genome
  - identity_percentage_cutoff: Percentage of identity cutoff used when blasting asRNAs against input genome (a value between 0 and 1)

- Optional arguments:
  - -h, --help: Shows help about how the program usage
  - -t TAGS, --tags TAGS: Excel file that includes the locus/gene tags to compute the asRNAS
  - -r RECOMPUTE, --recompute RECOMPUTE Shift position when recomputing asRNAS with hits
  
**Important:**
The set of tags should follow the format as illustrated in the sample file tags_k12.xlsx. That is, the tags file should contain the following headers and in that order:

Gene_Tag	Locus_Tag![image](https://user-images.githubusercontent.com/72103416/110243569-faea2380-7f28-11eb-8f78-f9b12e555c56.png)

**Output:**
The program will export the computed asRNAs into an excel file which will be located in the subdirectory *sequences* under the *srna* directory. The name of the output file has this format sequence_file_datetime_asrna.xlsx. In addition, CAREn also exports the set of gene and locus tags of the asRNAs that contain offset-hits in the genome. Notice that this file could be used as input to the program to recompute asRNAs (parameter -T). The name of the tag file has this format sequence_file_datetime_tags.xlsx.


**Examples:**

The directory *sequences* in the directory *srna* contains some sample genome sequences: K12.gb and JWGZ01.1.gbff. It also includes some samples of locus/gene tags for these sequences.

- Example 1:
Suppose that you would like to compute **all** asRNAs for the sequence K12.gb. The format for this file is genbank. The offset position to compute the asRNAS would be -10 and the length of the asRNAS would be 19. The expected_cutoff and identity_percentage_cutoff  are 0.01 and 0.8, respectively. Additionally, for the asRNAs that contain offset-hits in the genome, the program should recompute these asRNAS with a position of -15. Assume that the program was installed in the following path: /home/srna. Therefore, the sequences are located at directory /home/srna/sequences. For this example, the program should be executed like this:

python main.py /home/srna/sequences/K12.gb genbank -10 19 0.01 0.8 -r -15

The output of the program could look like this: K12.gb_03-08-2021 15:55:32_srna.xlsx and K12.gb_03-08-2021 15:55:34_tags.xlsx.

- Example 2:
Suppose that you would like to compute the asRNAs for *a set* of gene/locus tags given in a file. The input file is JWGZ01.1.gbff and the set of tags is in tags_jwz.xlsx. The format for this file is genbank. The position to compute the asRNAs would be -8 and the length of the asRNAs would be 21. The expected_cutoff and identity_percentage_cutoff are 0.01 and 0.8, respectively. Additionally, for the asRNAs that contain hits in the genome, the program should recompute these asRNAs with a position of -10. Assume that the program was installed in the following path: /home/srna. Therefore, the sequences are located at directory /home/srna/sequences. For this example, the program should be executed like this:


python main.py /home/srna/sequences/JWGZ01.1.gbff genbank -8 21 0.01 0.8 -r -10 -t home/srna/sequences/tags_jwz.xlsx

The output of the program could look like this: JWGZ01.1.gbff_03-08-2021 16:03:55_srna.xlsx and JWGZ01.1.gbff_03-08-2021 16:03:55_tags.xlsx



### References

[1] https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=FAQ

[2]https://ase.tufts.edu/chemistry/walt/sepa/Activities/BLASTpractice.pdf









