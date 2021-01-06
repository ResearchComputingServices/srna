import sys, getopt
import os.path
import time


from sRNA_provider import sRNA_Provider

def help():
   print('python main.py -i <sequence_input_file> -f <sequence_format e.g. genbank, fasta> (-a <alphabet> optional) (-s <shift> optional (-l <length> optional)')


def main(argv):

   start_program = time.time()
   alphabet = 'NULL'
   position = -8
   length = 21
   e_cutoff = 0.01
   identity_perc_cutoff = 0.8
   tags_filename = None

   try:
      opts, args = getopt.getopt(argv,"d:i:f:p:l:e:m:t:r:")
   except getopt.GetoptError:
      help()
      sys.exit(2)
   for opt, arg in opts:
      if opt in ("-d", "--base directory"):
         base_directory = arg
      else:
         if opt in ("-i", "--seq_file name"):
            seq_file = arg
         else:
            if opt in ("-f", "--format"):
               format = arg
            else:
               if opt in ("-a", "--alphabet"):
                  alphabet = arg
               else:
                  if opt in ("-p", "--position"):
                     position = arg
                  else:
                     if opt in ("-l", "--length"):
                        length = arg
                     else:
                        if opt in ("-e", "--expected cutoff"):
                           e_cutoff = arg
                        else:
                           if opt in ("-m", "--identity percentage cutoff"):
                              identity_perc_cutoff = arg
                           else:
                              if opt in ("-t", "--gene/locus tag filename"):
                                 tags_filename = arg
                              else:
                                 if opt in ("-r", "--recompute position"):
                                    recompute_position = arg

   #sRNA_provider = sRNA_Provider()
   #list = sRNA_provider.fetch_input_sequence('CP060121', base_directory)

   file_sequence = base_directory+'/'+seq_file

   if os.path.isfile(file_sequence):
      print('\nInput parameters: \n')
      print('Input file:', file_sequence)
      print('Format: ', format)
      print ('Position for sRNA computation: ',position)
      print('Length sRNA: ', length)
      print ('Expected cutoff (blast): ', e_cutoff)
      print ('Percentage Identity cutoff (blast): ', identity_perc_cutoff)

      print ('\n\n')

      if alphabet != 'NULL':
         print('Alphabet: ', alphabet)

      if int(position)==0:
         print ('Position cannot be 0. It is either >0 or <0.')
         sys.exit(2)

      #Read Sequence
      print('1. Reading input sequence\n')
      sRNA_provider = sRNA_Provider()
      sequence_record_list = sRNA_provider.read_input_sequence(file_sequence, format, alphabet)

      #Print sequence
      #sRNA_provider.print_input_sequence(sequence_record_list)
      print ('Total Records: ', len(sequence_record_list))

      #Compute sRNAS
      if tags_filename:
         file_tags = base_directory + '/' + tags_filename
         if os.path.isfile(file_tags):
             print('2. Computing sRNAs - locus/gene tags \n')
             gene_tags, locus_tags = sRNA_provider.load_tags(file_tags)
             list_sRNA = sRNA_provider.compute_sRNAs_from_genome(sequence_record_list, int(position), int(length),gene_tags, locus_tags)
         else:
             print ("Tag file could not be found!")
             sys.exit(2)
      else:
         print('2. Computing sRNAs\n')
         list_sRNA = sRNA_provider.compute_sRNAs_from_genome(sequence_record_list, int(position), int(length))

      print ('Total of computed sRNAs', sRNA_provider.total_srnas(list_sRNA))
      print('\n')

      start = time.time()
      print ('3. Blast each sRNA against input sequence\n')
      sRNA_provider.blast_sRNAs_against_genome(list_sRNA, e_cutoff, identity_perc_cutoff)
      end = time.time()
      print ('Elapsed minutes blast 1: ')
      print ((end-start)/60)

      #print ('4. Printing sRNA info \n')
      #sRNA_provider.print_list_srna(list_sRNA,sequence_record_list)


      print ('5. Get sRNA with hits \n')
      list_sRNA_with_hits = sRNA_provider.get_sRNAs_with_hits(list_sRNA)


      print ('6. Recompute sRNAs for sRNAs with hits')
      list_sRNA_recomputed = sRNA_provider.recompute_sRNAs(list_sRNA_with_hits, 1, int(recompute_position), int(length))

      start = time.time()
      print ('7. Blast the re-computed sequence')
      sRNA_provider.blast_sRNAs_against_genome(list_sRNA_recomputed, e_cutoff, identity_perc_cutoff)
      end = time.time()
      print('Elapsed minutes blast 2: ')
      print((end - start) / 60)

      #print ('. Printing sRNA info \n')
      #sRNA_provider.print_list_srna(list_sRNA_recomputed,sequence_record_list)

      print('8. Export Info')
      #sRNA_provider.export_sRNAs(list_sRNA, base_directory, seq_file, format, position, length, e_cutoff, identity_perc_cutoff)
      sRNA_provider.export_sRNAs(list_sRNA_recomputed, base_directory, seq_file, format, position, length, e_cutoff,identity_perc_cutoff)

      print ('9. Write tags to file')
      sRNA_provider.write_tags_to_file(base_directory,list_sRNA_recomputed)

      end_program = time.time()
      print('Elapsed minutes program in mins: ')
      print((end_program - start_program) / 60)


   else:
      print ('\n\n File:', seq_file, 'not found! \n')


if __name__ == "__main__":
   if len(sys.argv) < 5:
      help()
   else:
      main(sys.argv[1:])

