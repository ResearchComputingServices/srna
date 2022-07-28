from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
from Bio.SeqFeature import SeqFeature, FeatureLocation
import sys
import traceback
from sRNA_class import sRNA_Class
from blast import Blast
import pandas as pd
from io import BytesIO
import json
import os.path
from datetime import datetime
import uuid
import os
import time



DEBUG=3


class sRNA_Provider:

    blastProvider = Blast()

    def fetch_and_save_input_sequence(self, accession, base_directory):
        filename = base_directory + '/' + accession + ".gbk"
        if not os.path.isfile(filename):
            # Downloading...
            net_handle = Entrez.efetch(
                db="nucleotide", id=accession, rettype="gb", retmode="text"
            )
            out_handle = open(filename, "w")
            out_handle.write(net_handle.read())
            out_handle.close()
            net_handle.close()
            print("Saved")

        print("Parsing...")
        seq_iterator = SeqIO.parse(filename, "genbank")
        seq_record_list = list(seq_iterator)
        return seq_record_list


    def fetch_input_sequence(self, accession):
        try:
            # Downloading sequence
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
            seq_iterator = SeqIO.parse(handle, "genbank")
            seq_record_list = list(seq_iterator)
        except:
            seq_record_list = []
            print("Unexpected error at fetch_input_sequence:", sys.exc_info()[0])
            traceback.print_exc()

        return seq_record_list


    #Reads a data file and return a list of seq_record
    def read_input_sequence(self, input_file, format):

        try:
            seq_iterator = SeqIO.parse(input_file,format)
            seq_record_list = list(seq_iterator)
        except:
            seq_record_list=[]
            print("Unexpected error at read_sequence:", sys.exc_info()[0])
            traceback.print_exc()

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

            if start_sRNA_position<0:
                start_sRNA_position = 0

            if end_sRNA_position>len(sequence):
                end_sRNA_position = len(sequence)

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

            if start_sRNA_position<0:
                start_sRNA_position = 0

            if end_sRNA_position > len(sequence):
                end_sRNA_position = len(sequence)

            # Note that end_SRNA_position is not inclusive. So the real sequence goes up to end_SRNA_position-1
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



    def __get_sRNA_from_input_listCDS(self, record_index, seq_record, position, length, gene_tags_input, locus_tags_input):

        list_sRNA = []

        for feature in seq_record.features:
            if feature.type == 'CDS':
                gene = ''
                locus_tag = ''

                if 'gene' in feature.qualifiers.keys() and len(feature.qualifiers['gene'])>0:
                    gene = feature.qualifiers['gene']

                if 'locus_tag' in feature.qualifiers.keys() and len(feature.qualifiers['locus_tag'][0]) > 0:
                    locus_tag = feature.qualifiers['locus_tag']

                compute  = False
                for g in gene:
                    if g.lower() in gene_tags_input:
                        compute = True
                        break

                if not compute:
                    for l in locus_tag:
                        if l.lower() in locus_tags_input:
                            compute = True
                            break

                if compute:

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
        sRNA.input_record = sRNA_original.input_record

        return (sRNA)


    def __blast_sRNA_against_genome(self, list_sRNA, e_cutoff, identity_perc_cutoff):

        #Creates temporary files
        query_file = str(uuid.uuid4()) + '.fasta'
        subject_file = str(uuid.uuid4()) + '.fasta'

        #query_file = "seq_sRNA.fasta"
        #subject_file = "seq_Source.fasta"

        for index, srna in enumerate(list_sRNA):
            SeqIO.write(srna.input_sequence, subject_file, "fasta")
            seq_sRNA = SeqRecord(srna.sequence_sRNA, id="sRNA")
            SeqIO.write(seq_sRNA, query_file, "fasta")
            if len(srna.sequence_sRNA)>0 and len(srna.input_sequence)>0 and len(srna.list_hits)==0:
                    list_hits = self.blastProvider.blast(query_file, subject_file, str(srna.sequence_sRNA),float(e_cutoff), float(identity_perc_cutoff))
                    srna.list_hits = list_hits

        #Remove temporary files
        if os.path.exists(query_file):
            os.remove(query_file)

        if os.path.exists(subject_file):
            os.remove(subject_file)




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

    def _get_sRNAs_without_blast_info(self,list):
        list_sRNA_without_blast = []

        for seq_record in list:
            list = []
            for sRNA in seq_record:
                if len(sRNA.list_hits) == 0:
                    list.append(sRNA)

            if len(list) > 0:
                list_sRNA_without_blast.append(list)

        return list_sRNA_without_blast

    def _get_sRNAs_without_hits(self, list_sRNA):
        list_sRNA_without_hits =[]

        for seq_record in list_sRNA:
            list = []
            for sRNA in seq_record:
                if len(sRNA.list_hits)==1:
                    list.append(sRNA)

            if len(list)>0:
                list_sRNA_without_hits.append(list)

        return list_sRNA_without_hits


    def get_sRNAs_without_hits(self, list_sRNA, list_sRNA_recomputed):
        srnas_without_hits = []
        list = self._get_sRNAs_without_hits(list_sRNA)
        for srna in list:
            if srna != []:
                srnas_without_hits.append(srna)
        list = self._get_sRNAs_without_hits(list_sRNA_recomputed)
        for srna in list:
            if srna != []:
                srnas_without_hits.append(srna)

        return srnas_without_hits



    def get_not_processed_tags(self, list_sRNA, gene_tags_input, locus_tags_input):
        gene_tags_found = []
        locus_tags_found = []

        gene_tags_not_found = []
        locus_tags_not_found = []

        for record in list_sRNA:
            for srna in record:
                for g in srna.gene:
                    if g not in gene_tags_found:
                        g = g.strip()
                        gene_tags_found.append(g.lower())

                for l in srna.locus_tag:
                    if l not in locus_tags_found:
                        l = l.strip()
                        locus_tags_found.append(l.lower())

        for tag in gene_tags_input:
            if tag not in gene_tags_found:
                gene_tags_not_found.append(tag)

        for tag in locus_tags_input:
            if tag not in locus_tags_found:
                locus_tags_not_found.append(tag)

        return gene_tags_not_found, locus_tags_not_found


    def get_gene_locus_tags_correlated(self, list_sRNA):
        list_dict = []
        for record in list_sRNA:
            for srna in record:
                dict = {}
                genes = ''
                locus = ''
                for g in srna.gene:
                    genes = genes + g + ' '

                for l in srna.locus_tag:
                    locus = locus + l + ' '

                dict['Gene_Tag'] = genes
                dict['Locus_Tag'] = locus
                list_dict.append(dict)

        return list_dict


    def get_hit_tags_correlated(self, list_srna):

       list_dict = []
       for srnas in list_srna:
           index = 0
           total = len(srnas)
           while (index<total):
               hits_recomputed_srna = len(srnas[index+1].list_hits)
               #Safety check
               if (srnas[index].start_position_CDS==srnas[index+1].start_position_CDS) and hits_recomputed_srna>1:
                   dict = {}
                   genes = ''
                   locus = ''
                   for g in srnas[index].gene:
                       genes = genes + g + ' '

                   for l in srnas[index].locus_tag:
                       locus = locus + l + ' '

                   dict['Gene_Tag'] = genes
                   dict['Locus_Tag'] = locus
                   list_dict.append(dict)
               index = index + 2

       return list_dict



    def sRNA_hit_to_dict(self, srna,  display_if_hit, hit=None):
        dict={}
        display = True

        if hit:
            is_a_hit = True
            if srna.strand == -1 and srna.start_position_sRNA + 1 == int(hit.sbjct_start):
                is_a_hit = False
            if srna.strand == 1 and srna.start_position_sRNA + 1 == int(hit.sbjct_end):
                is_a_hit = False

            display = True
            if display_if_hit:
                if is_a_hit:
                    display = True
                else:
                    display = False

        if display:
            dict["Record"] = str(srna.input_record+1)
            dict["asRNA"] = str(srna.sequence_sRNA)
            dict["Strand"] = srna.strand
            dict["Gene"] = srna.gene
            dict["Locus Tag"] = srna.locus_tag
            dict["asRNA Start"] = self._start_in_file(srna.start_position_sRNA)
            dict["asRNA End"] = self._end_in_file(srna.end_position_sRNA)
            dict["asRNA Length"]= srna.length_sRNA
            dict["Offset Position"] = srna.shift
            dict["CDS Start"]= self._start_in_file(srna.start_position_CDS)
            dict["CDS End"]= self._end_in_file(srna.end_position_CDS)
            if hit:
                if display:
                    dict["Hit Start"]= int(hit.sbjct_start)
                    dict["Hit End"]= hit.sbjct_end
                    dict["Expected Value"]= hit.expect
                    dict["Align Length"]= hit.align_length
                    dict["Perc. Identity"]= (hit.score/srna.length_sRNA)*100
                    if is_a_hit:
                        dict["Description"] = "Offset-Hit"
                    else:
                        dict["Description"] = "asRNA"
            else:
                dict["Hit Start"] = 'Blast did not return any output for the selected threshold and percentage of identity parameters.'
                dict["Hit End"] = ''
                dict["Expected Value"] = ''
                dict["Align Length"] = ''
                dict["Perc. Identity"] = ''

        return dict


    def write_tags_to_file(self, base_directory, seq_file, list_sRNA, excel_output=True, csv_output=False):
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

        if len(gene_tags)!=len(locus_tags):
            if len(gene_tags)>len(locus_tags):
                i = len(locus_tags)
                while (i<len(gene_tags)):
                    locus_tags.append('')
                    i= i +1

            if len(locus_tags) > len(gene_tags):
                i = len(gene_tags)
                while (i < len(locus_tags)):
                    gene_tags.append('')
                    i = i + 1

        current = datetime.now()
        date_time = current.strftime("%m-%d-%Y %H:%M:%S")

        name = seq_file
        file_name = base_directory + '/' + name + '_' + date_time + '_tags'
        dict = {'Gene_Tag': gene_tags, 'Locus_Tag': locus_tags}
        df = pd.DataFrame(dict)
        if excel_output:
            df.to_excel(file_name+'.xlsx', header=True, index=False)
            print('Exported tags to: {0}'.format(file_name+'.xlsx'))
        if csv_output:
            df.to_csv(file_name+'.csv', header=True, index=False)
            print('Exported tags to: {0}'.format(file_name+'.csv'))
            



    def load_locus_gene_tags(self, file_name):
        try:
            gene_tags = []
            locus_tags = []

            df = pd.read_excel(file_name)
            gene_tags_input = df['Gene_Tag'].values.tolist()
            locus_tags_input = df['Locus_Tag'].values.tolist()


            for tag in gene_tags_input:
                if str(tag)!='nan':
                    tag = tag.lower()
                    tag = tag.strip()
                    gene_tags.append(tag)


            for tag in locus_tags_input:
                if str(tag) != 'nan':
                    tag = tag.lower()
                    tag = tag.strip()
                    locus_tags.append(tag)

        except:
            print("Unexpected error at load_locus_gene_tags:", sys.exc_info()[0])
            traceback.print_exc()

        return gene_tags, locus_tags


    def sRNAs_to_data_frames(self, list_sRNA, title, seq_file, name_seq, description_seq, format, position, position_rec, length, e_cutoff, perc_identity, display_only_hits, msg):

        name = seq_file

        # Creates panda data frame with general info
        general_info = {
            "Result " : title,
            "User Parameters": '',
            "Input File": seq_file,
            "Sequence Name": name_seq,
            "Sequence Description" : description_seq,
            "Format": format,
            "Position": position,
            "Length": length,
            "Expected cutoff": e_cutoff,
            "Percentage Indentity": perc_identity,
            "Position for Recomputing": position_rec
        }

        df_ginfo = pd.DataFrame.from_dict(general_info, orient='index')
        df_ginfo.rename(columns={0: ' '}, inplace=True)

        # Creates panda date frame with column headings
        headings = {
            "Seq. Record": '',
            "asRNA": '',
            "Strand": '',
            "Gene": '',
            "Locus Tag": '',
            "asRNA Start": '',
            "asRNA End": '',
            "asRNA Length": '',
            "Offset Position": '',
            "CDS Start": '',
            "CDS End": '',
            "Hit Start": '',
            "Hit End": '',
            "Expected Value": '',
            "Align Length": '',
            "Perc. Identity": '',
            "Description ": '',
        }

        df_headings = pd.DataFrame([headings])

        # Creates panda data frame with sRNA information
        list_rows = []
        index_record = 1
        for record in list_sRNA:
            for srna in record:
                if len(srna.list_hits) == 0:
                    dict = self.sRNA_hit_to_dict(srna, display_only_hits)
                    if dict!={}:
                        list_rows.append(dict)

                for hit in srna.list_hits:
                    dict = self.sRNA_hit_to_dict(srna, display_only_hits, hit)
                    if dict!={}:
                        list_rows.append(dict)
            index_record = index_record + 1

        if len(list_rows)==0:
            list_rows.append(msg)
        df_rows = pd.DataFrame(list_rows)

        return df_ginfo, df_headings, df_rows

    def export_output_to_file(self, base_directory, seq_file, name_seq, description_seq, format, position, position_rec, length, e_cutoff, perc_identity, list_sRNA_recomputed=None, list_sRNA=None, gene_tags_input=None, locus_tags_input=None, input_tags_filename=None, excel_output=True, csv_output=False):

        current = datetime.now()
        date_time = current.strftime("%m-%d-%Y %H:%M:%S")

        name = seq_file
        base_filename = os.path.join(base_directory, f"{name}_{date_time}_asrna")
        excel_filename = base_filename + '.xlsx'

        # Writes frames to excel
        writer = pd.ExcelWriter(excel_filename, engine='xlsxwriter')


        #Check if there is blast info
        if list_sRNA_recomputed and len(list_sRNA_recomputed) > 0:
            l = self._get_sRNAs_without_hits(list_sRNA)
            all_srnas = l + list_sRNA_recomputed

            srnas_with_blast_info = 0
            for record in all_srnas:
                srnas_with_blast_info= srnas_with_blast_info + len(record)

        if srnas_with_blast_info > 0:
            if list_sRNA:
                #Export good sRNAs: sRNAs without hits
                srnas_without_hits= self.get_sRNAs_without_hits(list_sRNA, list_sRNA_recomputed)
                title = 'Good asRNAs: asRNAs without off-target hits in the genome'
                msg = 'There are no good asRNAs. '
                df_ginfo, df_headings, df_rows = self.sRNAs_to_data_frames(srnas_without_hits, title, seq_file, name_seq, description_seq, format, position, position_rec, length, e_cutoff, perc_identity, False,msg)
                df_ginfo.to_excel(writer, sheet_name="Good_srnas", index=True)
                df_headings.to_excel(writer, startrow=13, sheet_name="Good_srnas", index=False)
                df_rows.to_excel(writer, sheet_name="Good_srnas", startrow=14, index=False, header=False)

            #Export hits
            if list_sRNA_recomputed:
                title = 'asRNAs with off-target hits in the genome'
                msg = 'There are no asRNAs with off-target hits.'
                df_ginfo, df_headings, df_rows = self.sRNAs_to_data_frames(list_sRNA_recomputed, title, seq_file, name_seq, description_seq, format, position, position_rec, length, e_cutoff, perc_identity, True,msg)
                df_ginfo.to_excel(writer, sheet_name="Hits", index=True)
                df_headings.to_excel(writer, startrow=13, sheet_name="Hits", index=False)
                df_rows.to_excel(writer, sheet_name="Hits", startrow=14, index=False, header=False)



            #Export all sRNAS
            if list_sRNA:
                title = 'All computed asRNAs and its hits.'
                msg = 'There are no computed asRNAs. Verify input sequence and/or input tags file (if provided).'
                if list_sRNA_recomputed and len(list_sRNA_recomputed) > 0:
                    l1 = self._get_sRNAs_without_hits(list_sRNA)
                    l2 = self._get_sRNAs_without_blast_info(list_sRNA)
                    all_srnas = l1 + l2 + list_sRNA_recomputed
                else:
                    all_srnas = list_sRNA
                df_ginfo, df_headings, df_rows = self.sRNAs_to_data_frames(all_srnas, title, seq_file,  name_seq, description_seq, format, position, position_rec, length, e_cutoff, perc_identity, False,msg)
                df_ginfo.to_excel(writer, sheet_name="All", index=True)
                df_headings.to_excel(writer, startrow=13, sheet_name="All", index=False)
                df_rows.to_excel(writer, sheet_name="All", startrow=14, index=False, header=False)

            #Export gene tags and locs tags of sRNAS without hits
            if srnas_without_hits:
                title = 'Pairs of gene and locus tags of good asRNAs (asRNAs without off-target hits)'
                description = 'The columns are correlated.'
                general_info=[]
                general_info.append(title)
                general_info.append(description)
                df_ginfo = pd.DataFrame(general_info)
                df_ginfo.to_excel(writer, sheet_name="Good Tags",  header=False, index=False)

                tag_list = self.get_gene_locus_tags_correlated(srnas_without_hits)
                headerF = True
                if len(tag_list)==0:
                    tag_list.append('There are no good asRNAs. ')
                    headerF = False

                df = pd.DataFrame(tag_list)
                df.to_excel(writer, sheet_name="Good Tags", startrow=4, header=headerF, index=False)

            # Export gene tags and locs tags of sRNAS with hits
            if list_sRNA_recomputed:
                title = 'Pairs of gene and locus tags of asRNAs WITH off-target hits'
                description = 'The columns are correlated.'

                general_info = []
                general_info.append(title)
                general_info.append(description)

                df_ginfo = pd.DataFrame(general_info)
                df_ginfo.to_excel(writer, sheet_name="Bad Tags",  header=False, index=False)

                tag_list = self.get_hit_tags_correlated(list_sRNA_recomputed)
                headerF= True
                if len(tag_list)==0:
                    tag_list.append('There are no asRNAs with off-target hits.')
                    headerF = False

                df = pd.DataFrame(tag_list)
                df.to_excel(writer, sheet_name="Bad Tags", startrow=4, header=headerF, index=False)
        else:
            title = 'All computed asRNAs.'
            msg = 'There are no computed asRNAs. Verify input sequence and/or input tags file (if provided).'
            all_srnas = list_sRNA
            df_ginfo, df_headings, df_rows = self.sRNAs_to_data_frames(all_srnas, title, seq_file, name_seq,
                                                                       description_seq, format, position, position_rec,
                                                                       length, e_cutoff, perc_identity, False, msg)
            df_ginfo.to_excel(writer, sheet_name="All", index=True)

            df_headings.to_excel(writer, startrow=13, sheet_name="All", index=False)
            df_rows.to_excel(writer, sheet_name="All", startrow=14, index=False, header=False)



        # Export which genes were not found from input tags
        if gene_tags_input and len(gene_tags_input)>0 and locus_tags_input and len(locus_tags_input)>0:
            gene_tags_not_found, locus_tags_not_found = self.get_not_processed_tags(list_sRNA,gene_tags_input,locus_tags_input)
            title = 'Input tags file: ' + input_tags_filename

            if len(gene_tags_not_found)==0 and len(locus_tags_not_found)==0:
                description = 'All gene and/or locus tags were found in the genome, and therefore at least one asRNA was computed for those tags.'
            else:
                description = 'The following genes and/or locus tags were NOT found in the genome and therefore asRNAs were NOT computed for those tags. Notice: the columns are NOT correlated. '

            general_info = []
            general_info.append(title)
            general_info.append(description)

            df_ginfo = pd.DataFrame(general_info)
            df_ginfo.to_excel(writer, sheet_name="Input Tags Result", header=False, index=False)

            headings = {
                "Gene_Tag": '',
                " ": ' ',
                "Locus_Tag ": '',
            }

            df_headings = pd.DataFrame([headings])
            df_headings.to_excel(writer, startrow=4, sheet_name="Input Tags Result", index=False)
            df = pd.DataFrame(gene_tags_not_found)
            df.to_excel(writer, sheet_name="Input Tags Result", startrow=5, startcol=0, header=False, index=False)
            df = pd.DataFrame(locus_tags_not_found)
            df.to_excel(writer, sheet_name="Input Tags Result", startrow=5, startcol=2, header=False, index=False)

        workbook = writer.book
        workbook.formats[0].set_font_size(14)

        # Close the Pandas Excel writer and output the Excel file.
        writer.save()


        # Convert Excel to CSV, if CSV output was requested
        if csv_output:
            sheets = pd.read_excel(excel_filename, sheet_name=None)
            for sheet_name in sheets.keys():
                csv_filename = f"{base_filename}_{sheet_name}.csv"
                sheets[sheet_name].to_csv(csv_filename, header=False, index=False)
            print (f"Exported sRNA output to: {base_filename}_*.csv")

        if excel_output:
            print ('Exported sRNA output to: {0}'.format(excel_filename))
        else:

            # If user only asked for CSV output, the Excel file was
            # only temporarily needed for conversion to CSV.  In this
            # case, delete the Excel file
            os.unlink(excel_filename)


    def compute_srnas(self, base_directory, seq_file, file_sequence_fullpath, format, position, length, e_cutoff, identity_perc_cutoff, file_tags, recompute_position, excel_output, csv_output):

        start_program = time.time()
        # Read Sequence
        print('Reading input sequence\n')
        sequence_record_list = self.read_input_sequence(file_sequence_fullpath, format)

        if sequence_record_list and not len(sequence_record_list) > 0:
            print('The specified format does not correspond to the input sequence.')
            sys.exit()

        name_seq = sequence_record_list[0].name
        description_seq = sequence_record_list[0].description
        print(name_seq)
        print(description_seq)

        # Print sequence
        #self.print_input_sequence(sequence_record_list)
        print('Total Records: ', len(sequence_record_list))

        # Compute sRNAS
        gene_tags=[]
        locus_tags=[]
        if file_tags:
            if os.path.isfile(file_tags):
                print('Computing sRNAs for locus/gene tags \n')
                gene_tags, locus_tags = self.load_locus_gene_tags(file_tags)
                if len(gene_tags)==0 and len(locus_tags)==0:
                    print ('The file with gene tags is invalid. Please verify the format.')
                    sys.exit(2)
                list_sRNA = self.compute_sRNAs_from_genome(sequence_record_list, int(position), int(length),gene_tags, locus_tags)
            else:
                print("Tag file could not be found!")
                sys.exit(2)
        else:
            print('Computing all sRNAs\n')
            list_sRNA = self.compute_sRNAs_from_genome(sequence_record_list, int(position), int(length))

        print('Total of computed sRNAs', self.total_srnas(list_sRNA))
        print('\n')

        print('Blast each sRNA against input sequence\n')
        self.blast_sRNAs_against_genome(list_sRNA, e_cutoff, identity_perc_cutoff)

        #print ('Printing sRNA info \n')
        #self.print_list_srna(list_sRNA,sequence_record_list)

        print('Get sRNA with hits \n')
        list_sRNA_with_hits = self.get_sRNAs_with_hits(list_sRNA)

        if recompute_position:
            print('Recompute sRNAs for sRNAs with hits')
            list_sRNA_recomputed = self.recompute_sRNAs(list_sRNA_with_hits, 1, int(recompute_position), int(length))
            self.blast_sRNAs_against_genome(list_sRNA_recomputed, e_cutoff, identity_perc_cutoff)
        else:
            list_sRNA_recomputed = []

        print('Export Info')
        self.export_output_to_file(base_directory, seq_file, name_seq, description_seq, format, position, recompute_position, length, e_cutoff, identity_perc_cutoff, list_sRNA_recomputed, list_sRNA, gene_tags, locus_tags, file_tags, excel_output, csv_output)

        print('Write tags of sRNAs with hits to file')
        self.write_tags_to_file(base_directory,seq_file, list_sRNA_recomputed,excel_output, csv_output)

        end_program = time.time()
        print('Elapsed minutes program in mins: ')
        print((end_program - start_program) / 60)


