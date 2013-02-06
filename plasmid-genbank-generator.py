#import BioPython tools
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


#import csv tools and related module
import csv
import itertools

with open('construction-master-list.csv', 'rU') as inputfile:
    construction_dict = csv.DictReader(inputfile)

    construction_list = []
    
    construct_number = 1

    plasmid_sequence_list = []


    def get_construct_number(row):
        return row["construct"]

    def get_strategy(row):
        return row["strategy"]

    #convert the dictionary into a list of dictionary entries
    for row in construction_dict:
        construction_list.append(row)

##    for row in construction_list:
##        print row            

    for key, items in itertools.groupby(construction_list, key=get_construct_number):
        for subitems in items:
##            print subitems
            plasmid_sequence_list.append(subitems['sequence'])
            plasmid_name = subitems['construct name']
            plasmid_number = subitems['construct']
            plasmid_filename = str(plasmid_number) + '. ' + str(plasmid_name) + '.gb'
        plasmid_sequence_text = str(' '.join(plasmid_sequence_list))
##        print(plasmid_sequence_text)
##        plasmid_sequence = Seq(str(plasmid_sequence_text), IUPAC.unambiguous_dna)
        outputfile = open(plasmid_filename, 'w+')
        #convert plasmid_sequence into an actual sequence object
        plasmid_sequence_obj = Seq(plasmid_sequence_text, IUPAC.unambiguous_dna)
        plasmid_sequence_record = SeqRecord(plasmid_sequence_obj)
        SeqIO.write(plasmid_sequence_record, outputfile, 'genbank')
##        print '\n'
        #clear the plasmid sequence list
        plasmid_sequence_list = []
