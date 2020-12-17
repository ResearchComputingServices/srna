from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import sys
import traceback
from sRNA_class import sRNA_Class
from blast import Blast
import pandas as pd
from io import BytesIO
import json

DEBUG=100


class sRNA_Provider:

    blastProvider = Blast()

    #Reads a data file and return a list of seq_record
    def read_input_sequence(self, input_file, format, alphabet):
        try:
            if alphabet!='NULL':
                seq_iterator = SeqIO.parse(input_file,format,alphabet)
            else:
                seq_iterator = SeqIO.parse(input_file,format)
        except:
            print("Unexpected error at read_sequence:", sys.exc_info()[0])
            traceback.print_exc()

        seq_record_list = list(seq_iterator)
        return seq_record_list


    def __print_seq_record(self, seq_record):
        print('Record ID: ', seq_record.id)
        # repr is a "representative" output of the sequence
        print('Sequence: ', repr(seq_record.seq))
        print('Sequence Length:', "{:,}".format(len(seq_record)))
        # print ('Annotations: ', seq_record.annotations)
        #print('\n')


    #Print a list of sequences
    def print_input_sequence(self, seq_record_list):
        try:
            count = 0
            for seq_record in seq_record_list:
                self. __print_seq_record(seq_record)
                count = count + 1

            print('\nTotal Records: ', count)
        except:
            print("Unexpected error at print_sequence:", sys.exc_info()[0])
            traceback.print_exc()



    def __sign(self, number):
        if number>0:
            return '+'
        else:
            if number <0:
                return '-'
            else:
                return ' '

        # Python sequence start at 0 to end-1
        # However when we obtain a substring from a string the method is [start...end)
        # That is end is not inclusive
        # Therefore a start position in the file is position+1
        # But end position remains the same, because an end position actually is end-1 in python
        # python sequence
        # [start...end)
        # file sequence
        # [start+1...end]

    def _start_in_file(self, position):
        return position + 1

    def _end_in_file(self, position):
        return position



    def total_srnas(self, list_srna):
        total = 0
        for list in list_srna:
            total = total + len(list)
        return total



    def __print_srna(self, srna):
        print ('sRNA: ',srna.sequence_sRNA)


        print ('sRNA location: ['+str("{:,}".format(self._start_in_file(srna.start_position_sRNA))) + ' - ' + str("{:,}".format(srna.end_position_sRNA)) + ']')
        print ('sRNA Length: ',srna.length_sRNA)
        print ('CDS: [' + str("{:,}".format(self._start_in_file(srna.start_position_CDS))) + '-' + str("{:,}".format(srna.end_position_CDS)) + '] (' + self.__sign(srna.strand) + ')')
        if srna.gene!='':
            print ('Info: ', srna.gene)

        print('\n')

        if len(srna.list_hits)>0:
            print ('Total Hits: ',len(srna.list_hits))
            print ("Hit# \t\t Expect Value \t Per. Id. \t Align Length \t Hit Start - Hit End")
            for index, hit in enumerate(srna.list_hits):
                if srna.strand>0:
                    print(str(index+1) + "\t\t\t" + str(hit.expect) + "\t \t" + str((hit.score / srna.length_sRNA) * 100) + "\t \t" + str(
                    hit.align_length) + "\t \t" + "[" + str("{:,}".format(hit.sbjct_end)) + "-" + str(
                    "{:,}".format(hit.sbjct_start)) + "]")
                else:
                    print(str(index+1)  + "\t\t\t" + str(hit.expect) + "\t\t \t" + str((hit.score/srna.length_sRNA)*100) + "\t \t" + str(hit.align_length) + "\t \t \t" + "[" + str("{:,}".format(hit.sbjct_start)) + "-" + str("{:,}".format(hit.sbjct_end)) + "]")
        else:
            print ('No hits found.')

        print('\n')


    def __print_list_srna(self, list_sRNA):
        for index, srna in enumerate(list_sRNA):
                print ('\n')
                print("sRNA #" + str(index + 1))
                self.__print_srna(srna)



    def print_list_srna(self, list_sRNA, seq_record_list):
        for record_index, list_sRNA_per_record in enumerate(list_sRNA):
                print('----------------------------------------------------------------------------')
                print ("Record: # ", record_index+1)
                self.__print_seq_record(seq_record_list[record_index])
                print ('Total sRNAs: ', len(list_sRNA_per_record))
                self.__print_list_srna(list_sRNA_per_record)




    def sRNA_Forward(self, sequence, start_CDS, end_CDS, strand, gene, locus_tag, position, length):
        sRNA = sRNA_Class()
        try:

            shift = position
            start_gene = start_CDS

            #There is no 0 position
            #-1 is the position immediately before the start of the gene
            #+1 is the first position in the gene
            #Therefore if +1 the position is actually start_gene
            if position > 0:
                shift = shift - 1

            mid_position = start_gene + shift
            half_length =  int((length-1)//2)    #We substract -1 due to the mid_position
            start_sRNA_position = mid_position - half_length
            end_sRNA_position = start_sRNA_position + length

            #Note that end_SRNA_position is not inclusive. So the real sequence goes up to end_SRNA_position-1
            sub_sequence = sequence[start_sRNA_position:end_sRNA_position]
            sub_sequence = sub_sequence.reverse_complement()
            sRNA = sRNA_Class(start_sRNA_position, end_sRNA_position, len(sub_sequence), sub_sequence, start_CDS, end_CDS, position, strand, gene, locus_tag)
        except:
            print("Unexpected error at get_sub_sequence:", sys.exc_info()[0])
            traceback.print_exc()
        return sRNA




    def sRNA_Complement(self, sequence, start_CDS, end_CDS, strand, gene, locus_tag, position, length):
        sRNA = sRNA_Class()
        try:
            if position < 0:
                shift = - position
            else:
                shift = -(position-1)

            # We need to remove -1 because the .end position is not inclusive in python
            # That is, a sequence starts with 1 in file, but 0 in python
            # but feature.location.end gives [starts,end] where end is NOT an inclusive
            # character.

            start_gene = end_CDS - 1
            mid_position = start_gene + shift
            half_length =  int((length-1)//2)    #We substract -1 due to the mid_position
            start_sRNA_position = mid_position - half_length
            end_sRNA_position = start_sRNA_position + length
            sub_sequence = sequence[start_sRNA_position:end_sRNA_position]
            sRNA = sRNA_Class(start_sRNA_position, end_sRNA_position, len(sub_sequence), sub_sequence, start_CDS, end_CDS, position, strand,
                              gene, locus_tag)
        except:
            print("Unexpected error at get_sub_sequence:", sys.exc_info()[0])
            traceback.print_exc()
        return sRNA




    def __get_sRNA_from_input_all_CDS(self, record_index, seq_record, position, length):

        list_sRNA = []
        for feature in seq_record.features:
            if feature.type == 'CDS':
                gene = ''
                locus_tag=''
                if 'gene' in feature.qualifiers.keys() and len(feature.qualifiers['gene'])>0:
                    gene = feature.qualifiers['gene']

                if 'locus_tag' in feature.qualifiers.keys() and len(feature.qualifiers['locus_tag'][0])>0:
                    locus_tag = feature.qualifiers['locus_tag']

                sequence = seq_record.seq
                if feature.location.strand == 1:  # Forward
                    sRNA = self.sRNA_Forward(sequence, feature.location.start, feature.location.end, feature.location.strand, gene, locus_tag, position, length)

                else:
                    if feature.location.strand == -1:  # Reverse
                        sRNA = self.sRNA_Complement(sequence, feature.location.start, feature.location.end, feature.location.strand, gene, locus_tag, position, length)

                list_sRNA.append(sRNA)
                sRNA.input_sequence = seq_record
                sRNA.input_record = record_index
        return(list_sRNA)



    def __get_sRNA_from_input_listCDS(self, record_index, seq_record, position, length, gene_tags, locus_tags):

        list_sRNA = []
        for feature in seq_record.features:
            if feature.type == 'CDS':
                gene = ''
                locus_tag = ''
                if 'gene' in feature.qualifiers.keys() and len(feature.qualifiers['gene'])>0:
                    gene = feature.qualifiers['gene']

                if 'locus_tag' in feature.qualifiers.keys() and len(feature.qualifiers['locus_tag'][0])>0:
                    locus_tag = feature.qualifiers['locus_tag']

                if gene[0] in gene_tags or locus_tag[0] in locus_tags:

                    sequence = seq_record.seq
                    if feature.location.strand == 1:  # Forward
                        sRNA = self.sRNA_Forward(sequence, feature.location.start, feature.location.end,
                                             feature.location.strand, gene, locus_tag, position, length)

                    else:
                        if feature.location.strand == -1:  # Reverse
                            sRNA = self.sRNA_Complement(sequence, feature.location.start, feature.location.end,
                                                    feature.location.strand, gene, locus_tag, position, length)

                    sRNA.input_sequence = seq_record
                    sRNA.input_record = record_index
                    list_sRNA.append(sRNA)

        return (list_sRNA)



    def compute_sRNAs_from_genome(self, seq_record_list, position, length, gene_tags=None, locus_tags=None):

        list_sRNA = []
        for record_index, seq_record in enumerate(seq_record_list):
            if gene_tags or locus_tags:
                list_sRNA_perRecord = self.__get_sRNA_from_input_listCDS(record_index, seq_record, position, length, gene_tags, locus_tags)
            else:
                list_sRNA_perRecord = self.__get_sRNA_from_input_all_CDS(record_index, seq_record, position, length)

            list_sRNA.append(list_sRNA_perRecord)

        return(list_sRNA)    #This is a list of lists. [i][j] The jth sRNA for the ith seqRecord




    def __recompute_sRNA(self, sRNA_original, position, length):
        if sRNA_original.strand == 1:  # Forward
            sRNA = self.sRNA_Forward(sRNA_original.input_sequence.seq, sRNA_original.start_position_CDS,
                                     sRNA_original.end_position_CDS, sRNA_original.strand, sRNA_original.gene, sRNA_original.locus_tag, position,
                                     length)
        else:
            if sRNA_original.strand == -1:  # Reverse
                sRNA = self.sRNA_Complement(sRNA_original.input_sequence.seq, sRNA_original.start_position_CDS,
                                            sRNA_original.end_position_CDS, sRNA_original.strand, sRNA_original.gene, sRNA_original.locus_tag,
                                            position, length)

        sRNA.input_sequence = sRNA_original.input_sequence
        sRNA.input_record = sRNA.input_record

        return (sRNA)



    def __blast_sRNA_against_genome(self, list_sRNA, e_cutoff, identity_perc_cutoff):

        query_file = "seq_sRNA.fasta"
        subject_file = "seq_Source.fasta"

        for index, srna in enumerate(list_sRNA):
            SeqIO.write(srna.input_sequence, subject_file, "fasta")
            seq_sRNA = SeqRecord(srna.sequence_sRNA, id="sRNA")
            SeqIO.write(seq_sRNA, query_file, "fasta")
            if len(srna.sequence_sRNA)>0 and len(srna.input_sequence)>0 and len(srna.list_hits)==0:
                    list_hits = self.blastProvider.blast(query_file, subject_file, str(srna.sequence_sRNA),float(e_cutoff), float(identity_perc_cutoff))
                    srna.list_hits = list_hits



    def blast_sRNAs_against_genome_for_DEBUG(self, list_sRNA, e_cutoff, identity_perc_cutoff):
        for list_sRNA_per_record in list_sRNA:
            list_sRNA_per_record_ = list_sRNA_per_record[0:DEBUG]    ####For debuging only!
            #list_sRNA_per_record_ = list_sRNA_per_record
            self.__blast_sRNA_against_genome(list_sRNA_per_record_, e_cutoff, identity_perc_cutoff)
        return (list_sRNA)




    def blast_sRNAs_against_genome(self, list_sRNA, e_cutoff, identity_perc_cutoff):
        for list_sRNA_per_record in list_sRNA:
            self.__blast_sRNA_against_genome(list_sRNA_per_record, e_cutoff, identity_perc_cutoff)
        return (list_sRNA)



    def get_sRNAs_with_hits(self, list_sRNA):
        list_sRNA_with_hits =[]

        for seq_record in list_sRNA:
            list = []
            for sRNA in seq_record:
                #Blasting always returns the sRNA subsequence as a hit
                if len(sRNA.list_hits)>1:
                    list.append(sRNA)
            list_sRNA_with_hits.append(list)

        return list_sRNA_with_hits


    def recompute_sRNAs(self, list_sRNA, times, position, length):
        list_recomputed_sRNAs = []

        for seq_record in list_sRNA:
            list = []
            for sRNA in seq_record:
                list.append(sRNA)  # Original sRNA
                i = 1
                while i <= times:
                    i = i + 1
                    sRNA_rec = self.__recompute_sRNA(sRNA, position * 1, length)
                    list.append(sRNA_rec)
            list_recomputed_sRNAs.append(list)

        return list_recomputed_sRNAs


    def sRNA_hit_to_dict(self, srna, hit):
        dict={}
        dict["sRNA"] = str(srna.sequence_sRNA)
        dict["Strand"] = srna.strand
        dict["Gene"] = srna.gene
        dict["Locus Tag"] = srna.locus_tag
        dict["sRNA Start"] = self._start_in_file(srna.start_position_sRNA)
        dict["sRNA End"] = self._end_in_file(srna.end_position_sRNA)
        dict["sRNA Length"]= srna.length_sRNA
        dict["Shift/Position"] = srna.shift
        dict["CDS Start"]= self._start_in_file(srna.start_position_CDS)
        dict["CDS End"]= self._end_in_file(srna.end_position_CDS)
        dict["Hit Start"]= int(hit.sbjct_start)
        dict["Hit End"]= hit.sbjct_end
        dict["Expected Value"]= hit.expect
        dict["Align Length"]= hit.align_length
        dict["Perc. Identity"]= (hit.score/srna.length_sRNA)*100
        return dict


    def write_tags_to_file(self, base_directory, list_sRNA):
        file_name = base_directory + '/tags.xlsx'

        gene_tags = []
        locus_tags = []

        for record in list_sRNA:
            for sRNA in record:

                for gene in sRNA.gene:
                    if gene not in gene_tags:
                        gene_tags.append(gene)


                for locus in sRNA.locus_tag:
                    if locus not in locus_tags:
                        locus_tags.append(locus)


        dict = {'Gene_Tag': gene_tags, 'Locus_Tag': locus_tags}
        df = pd.DataFrame(dict)
        df.to_excel(file_name, header=True, index=False)



    def load_tags(self, file_name):
        df = pd.read_excel(file_name)
        gene_tags = df['Gene_Tag'].values.tolist()
        locus_tags = df['Locus_Tag'].values.tolist()
        return gene_tags, locus_tags



    def export_sRNAs(self, list_sRNA, base_directory, seq_file, format, position, length, e_cutoff, perc_identity):

       file_name = base_directory + '/srna.xlsx'
       name= seq_file
       general_info = {
                    "User Parameters": '',
                    "Input File": seq_file,
                    "Format": format,
                    "Position": position,
                    "Length": length,
                    "Expected cutoff": e_cutoff,
                    "Percentage Indentity": perc_identity
                    }

       #df = pd.DataFrame([general_info])
       df_ginfo = pd.DataFrame.from_dict(general_info, orient='index')
       df_ginfo.rename(columns={0: ' '}, inplace=True)

       headings = {
           "sRNA": '',
           "Strand" :'',
           "Gene":'',
           "Locus Tag":'',
           "sRNA Start": '',
           "sRNA End": '',
           "sRNA Length": '',
           "Shift/Position":'',
           "CDS Start": '',
           "CDS End": '',
           "Hit Start": '',
           "Hit End" : '',
           "Expected Value" : '',
           "Align Length": '',
           "Perc. Identity": ''
       }

       df_headings = pd.DataFrame([headings])

       with pd.ExcelWriter(file_name, engine='xlsxwriter') as writer:
           df_ginfo.to_excel(writer, sheet_name="{} summary".format(name + 'sRNA'), index=True)
           df_headings.to_excel(writer, startrow=9, sheet_name="{} summary".format(name + 'sRNA'), index=False)

           row = 10
           for record in list_sRNA:
                for srna in record:
                    for hit in srna.list_hits:
                        dict = self.sRNA_hit_to_dict(srna, hit)
                        pd.DataFrame([dict]).to_excel(writer,startrow=row,sheet_name="{} summary".format(name + 'sRNA'), index=False, header=False)
                        row = row+1
                    row=row+1

       #output = BytesIO()
           workbook = writer.book
           worksheet = writer.sheets["{} summary".format(name + 'sRNA')]
           format = workbook.add_format()
           format.set_align('left')
           format.set_align('vcenter')
           #worksheet.set_column('A:A', 16, format)
           worksheet.set_column('A:N', 25, format)
           writer.save()



    def write_json(self, list, filename):


        list_dict = []
        for item in list:
            ini_list =[]
            for sRNA in item:
                if len(sRNA.list_hits)>0:
                    #ini_list.append(sRNA.__dict__)
                    ini_list.append(sRNA.to_dict())
            list_dict.append(ini_list)


        with open(filename, 'w') as outfile:
            json.dump(list_dict, outfile)

