# srna
Run syntax:

print('python main.py -i <sequence_input_file> -f <sequence_format e.g. genbank, fasta> (-a <alphabet> optional) (-s <shift> optional (-l <length> optional)')

-s: if not provided is 8
-l: if not provided is 21


Example of runs:

-i /Users/jazminromero/development/sequences/K12.gb -f genbank -s 8 -l 21

python main.py -i /Users/jazminromero/development/sequences/K12.gb -f genbank -s 8 -l 21


sRNA Documentation.

Basic algorithm for computing a sRNA

First, some background on the data file itself. A .gb file has two main parts: (1) annotations, giving the locations of different genomic features (such as genes, CDS, or tRNAs), and (2) the DNA sequence itself, which starts after the word "ORIGIN". We're interested in all of the "CDS" features; "CDS" stands for "coding sequence", which is a piece of DNA that codes for a protein.
 
For each CDS, we'll design a short RNA, or sRNA, sequence (for now; we may not want to design sRNAs for every CDS in the end, but it's probably easiest to start this way). Note that, for each CDS, the .gb file gives a location. For example, the thrA CDS is located from positions 337 to 2799.
 
CDS             337..2799
                     /gene="thrA"
                     /locus_tag="b0002"
                     ...
 
For current purposes, there are two varieties of CDS: those whose position is written like thrA, and those whose position information includes the word "complement", e.g., yaaA:
 
CDS             complement(5683..6459)
                     /gene="yaaA"
                     /locus_tag="b0006"
                     ...
 
"complement" means that the protein is encoded backwards. If "complement" is not indicated, then the protein is encoded forwards.
 
For genes encoded backwards ("complement"), do the following:
 
(1) Get the position of the start of the gene. This will be the second number - 6459 in the case of yaaA.
(2) In the DNA sequence, go to the position 8 letters (or "bases") to the right of the start position - 6467 for yaA.
(3) Capture a 21 base sequence with this position in the middle, e.g., positions 6457-6477 inclusive. This is the sRNA sequence for this gene.
 
For genes encoded forwards (not complement), do the following:
 
(1) Get the position of the start of the gene. This will be the first number - 337 for thrA.
(2) In the DNA sequence, go the the position 8 bases to the left of the start position - 329 for thrA.
(3) Capture a 21 base sequence with this position in the middle, here, 319-339 inclusive.
(4) Find the reverse complement to this sequence. You could do this with the reverse_complement method in BioPython. Alternatively, replace every A with T, G with C, T with A, and C with G, then reverse it. This is the sRNA sequence for this gene.
 


Computation of short RNA sequences in a input genome
