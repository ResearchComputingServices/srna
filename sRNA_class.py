from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

class sRNA_Class:
    input_sequence = Seq("", IUPAC.unambiguous_dna)
    input_record = 0
    start_position_sRNA = 0
    end_position_sRNA = 0
    length_sRNA = 0
    sequence_sRNA = Seq("", IUPAC.unambiguous_dna)
    start_position_CDS = 0
    end_position_CDS = 0
    shift = 0
    strand = 0
    gene = ''
    locus_tag=''
    list_hits = []



    def __init__(self, start_position_sRNA=0, end_position_sRNA=0, length_sRNA=0, sequence_sRNA=Seq(""), start_position_CDS=0, end_position_CDS=0, shift=0, strand=0, gene='', locus_tag=''):
        self.start_position_sRNA = start_position_sRNA
        self.end_position_sRNA = end_position_sRNA
        self.length_sRNA = length_sRNA
        self.sequence_sRNA = sequence_sRNA
        self.start_position_CDS = start_position_CDS
        self.end_position_CDS = end_position_CDS
        self.shift = shift
        self.strand = strand
        self.gene = gene
        self.locus_tag = locus_tag


    def to_dict(self):
        dict={}
        dict['sRNA'] = str(self.sequence_sRNA)
        dict["Strand"] = self.strand
        dict["Gene"] = self.gene
        dict["Locus Tag"] = self.locus_tag
        dict["sRNA Start"] = self.start_position_sRNA + 1
        dict["sRNA End"] = self.end_position_sRNA
        dict["sRNA Length"] = self.length_sRNA
        dict["Shift/Position"] = self.shift
        dict["CDS Start"] = self.start_position_CDS+1
        dict["CDS End"] = self.end_position_CDS
        list_hits =[]
        for hit in self.list_hits:
            dict_hit={}
            dict_hit["Hit Start"] = int(hit.sbjct_start)
            dict_hit["Hit End"] = hit.sbjct_end
            dict_hit["Expected Value"] = hit.expect
            dict_hit["Align Length"] = hit.align_length
            dict_hit["Perc. Identity"] = (hit.score / self.length_sRNA) * 100
            list_hits.append(dict_hit)
        dict["Hits"] = list_hits

        return dict



















