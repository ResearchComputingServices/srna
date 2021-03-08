
## Algorithm

### **Algorithm to compute a sRNAs for a gene**

The sRNA computation follows the algorithm provided by Prof. Alex Wong at Carleton University [1]. 

Along with other parameters, the algorithm receives as input a genome sequence in a file. This file should contain two parts:  (1) annotations which gives the locations of different genomic features (such as genes, CDS, or tRNAs), and (2) the DNA sequence itself. The most common formats for this file are genbank and fasta.

The algorithm will compute sRNAs from the CDS features. For simplicity, there are two varieties of CDS:  (1) the CDS whose position is written like: CDS *i..m*. This indicates that the gene in this CDS starts at the position *i* and ends at position *m*. We will refer to this CDS as encoded forward.

The second type of CDS are those whose position information includes the word "complement" CDS  *complemen(i..m)*.  The word "complement" means that the genome is encoded backwards. For complement CDS, *m* would be the start of the gene and *i* is the end of the gene. 

Given a CDS *i…m*, for a gene encoded forward, the steps to compute a sRNA are the following.
 
1. Get the position of the start of the gene. This will be the first number *i*. 
2. For a given shift position *p*,  
  - If *p<0*,  go *p* positions immediately before the start of the gene without including the i position. Let us denote this position as  *s = i + p*.   
  - If *p>0*,  go *p* positions at the right of the start of the gene including the start of the gene. That is, *s =i +p-1*.
 **Notice that there is no zero position p=0**. For a value of  *p=-1*, s would be the position immediately before the start of the gene (*s=i-1*), and for *p=1*, *s* would be the first position in the start of the gene (*s=i*).

3. For a given parameter *l*, capture a *l* sequence with the *s* position in the middle. Let *h = (l-1)/2*. The sequence would be captured at the positions: *s-h…s+h*. 

4. Find the reverse complement to this sequence. This could be done by replacing every A with T, G with C, T with A, and C with G, then reverse it.  **This is the sRNA sequence for the gene**.

Given a complement *CDS(i..m)* for a gene encoded backward, the steps to compute the sRNA are as follows:

1. Get the position of the start of the gene. This will be the second number: *m* 
2. For a given shift position p, 
  - If *p<0*,  go *p* positions before the start position without including the position *m*. Notice that since the gene is encoded backward, going before the start position corresponds to going to the right of the start position. Also notice that since *p<0*, for going to the right of the sequence, we should subtract *p*. That is: *s = m – p*.   
  - If *p>0*, go *p* positions after the start position including the start position. This would be equivalent to go *p* bases to the left of *m* including *m*, i.e., *s = m – p +1*.
Once again, there is no zero position *p=0*. For *p=-1*, *s* would be the position immediately before the start of the gene (*s=m+1*), and for *p=1*, *s=m* would be the first position in the start of the gene.

3. For a given parameter *l*, capture a *l* sequence with the *s* position in the middle. Let *h = (l-1)/2*. The sequence would be captured at the positions: *s-m…s+m*. **This is the sRNA sequence for the gene**.
 

### **Algorithm to compute a set of sRNAs from a genome**

We describe here the general steps behind the sRNA computation program. For a given input genome, the sRNA program can compute either the sRNAs for all CDS in the genome, or for a subset of these CDS. The user can specify for which CDS the sRNAs will be computed through a file of tags. Each tag  will correspond to a gene tag or a locus tag. If the tag is present in the genome, then a sRNA will be obtained for such tag as described in the previous section.

After the set of sRNAs is obtained, the next step will be to look for similar sequences to these sRNAs in the genome. This is to try to ensure that an sRNA only occur in the gene or CDS where the sRNA was obtained from. The similarity information will be obtained using BLAST [2].   
For each sRNA, BLAST will output (among other information ) a set of subsequences with an E value and a percentage of identity. Let us call these subsequences “hits”. For each pair (sRNA, hit), there is an E value and a percentage of identity returned by BLAST.  The E value describes how many times you would “expect” a match (sRNA, hit) to occur in the genome by chance, the closer to zero the value of E, the less likely is that the match (sRNA ,hit) will occur more than once. The percentage of identity describes how similar is the hit to the sRNA sequence. The higher the value of the percentage, the more similar are the hit and the sRNA [2].  The user will specify the E values (expected_cutoff) and the percentage of identity (identity_percentage_cutoff) for which BLAST will do a cutoff. That is, BLAST will return only the hits that pass those thresholds. 

Let us define as *R = [R1, R2,…,Rn]* the set of sRNAs which have at least more than hit with *E<=expected_cutoff* and *P>=identity_percentage_cutoff*. As an option, the user could specify if a second computation would take place. That is, if for each sRNA *Ri* in *R*, the program will identify from which CDS, *Ri* was obtained from and it will re-obtain a new sRNA for that CDS at a different shift position *p’*, where *p<>p’*. After that the program, will BLAST the new set of reobtained sRNAs. 
 
We summarize next the steps of the sRNAs computation program.

**Input:** 

A genome sequence | S 
A set of gene/locus tags T (optional)
A shift position | p 
A length value | l 
A cutoff expected value | expected_cutoff
A cutoff percentage of identity | identity_percentage_cutoff  
A second shift position | p’ p’<>p (optional)

**Output:**  a set of sRNAs R and a set of sRNAs R' (optional)

1.	If T is not provided (T=none), make T=all CDS in the genome S.
2.	For each CDS in T, compute the sRNAs following the steps as described in Section 1 using parameters p and l.
3.	For each computed sRNA, run BLAST with parameters expected_cutoff, and percentage_of_identity.
4.	Let R=[R1,…,Rn] be the subset of sRNAs which have more than one hit in the genome at the expected value and percentage of identity cutoffs.
5.	If p’ is none, return R.
6.	Otherwise, for each sRNA in Ri in R, determine from which CDS the Ri was obtained from and compute a new sRNA as described in Section 1 using p’ and l. Make this R' this new set of sRNAs.
7.	Return R'

## Program to compute sRNAs for a given input genome in python**

### How to run the program in your local machine**
  
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
 5. Install the requirements for the enviroment:
    *pip install -r requirements.txt*
 7. The program to be run is *main.py* and the parameters that it receives are the following:

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
  
**Important:
The set of tags should follow the format as illustrated in the sample file tags_k12.xlsx. That is, the tags file should contain the following headers and in that order:

Gene_Tag	Locus_Tag![image](https://user-images.githubusercontent.com/72103416/110243569-faea2380-7f28-11eb-8f78-f9b12e555c56.png)
 
 
**Examples:**

The directory *sequences* in the directory *srna* contains some sample genome sequences: K12.gb and JWGZ01.1.gbff. It also includes some samples of locus tags for these sequences.

- Example 1:
Suppose that you would like to compute all sRNAs for the sequence K12.gb. The format for this file is genbank. The position to compute the sRNAS would be -10 and the length of the sRNAS would be 19. The expected_cutoff and identity_percentage_cutoff  are 0.01 and 0.8, respectively. Additionally, for the sRNAs that contain hits in the genome, the program should recompute these sRNAS with a position of -15. Assume that the program was installed in the following path: /home/srna. Therefore, the sequences are located at directory /home/srna/sequences. For this example, the program should be executed like this:

python main.py /home/srna/sequences/K12.gb genbank -10 19 0.01 0.8 -r -15


- Example 2:
Suppose that you would like to compute the sRNAs for the set of tags given in a file. The input file is JWGZ01.1.gbff and the set of tags if in tags_jwz.xlsx. The format for this file is again genbank. The position to compute the sRNAS would be -8 and the length of the sRNAS would be 21. The expected_cutoff and identity_percentage_cutoff  are 0.01 and 0.8, respectively. Additionally, for the sRNAs that contain hits in the genome, the program should recompute these sRNAS with a position of -10. Assume that the program was installed in the following path: /home/srna. Therefore, the sequences are located at directory /home/srna/sequences. For this example, the program should be executed like this:


python main.py /home/srna/sequences/JWGZ01.1.gbff genbank -8 21 0.01 0.8 -r -10 -t home/srna/sequences/tags_jwz.xlsx









