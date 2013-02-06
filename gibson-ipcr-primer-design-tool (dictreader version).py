#import BioPython Tools
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

#import csv tools
import csv
import sys
import os
import itertools

#import other tools
import operator
from operator import itemgetter
import time

##from tempita import looper

#import Google Data API
import gdata.docs.service
import gdata.spreadsheet.service
import gdata.docs
import re, os

with open('constructs-to-make-shortened2.csv', 'rU') as constructs:
    construct_list = csv.DictReader(constructs)
    def get_construct_number(row):
        return row["Construct"]
    def get_strategy(row):
        return row["Strategy"]
    primer_list = []
    gibson_primer_temp_list = []
    counter = 1
    groups = []

##the case where Gibson (regular) is the strategy##
##===============================================##
    for key, items in itertools.groupby(construct_list, key=get_construct_number):
##        print(list(items))
##        print(key, items)
##        print(items_list)
##        print(zip(items))
##        print(len(list(zip(items))))
        counter = 1
        string_number = 1
        for subitems in items:
            if subitems['Strategy'] == 'Gibson' and subitems['Construct'] == key:
##                print(key)
                fp_anneal = Seq(subitems['Sequence'][0:40], IUPAC.unambiguous_dna)
                gibson_primer_temp_list.append([fp_anneal, key, counter, 'fp_anneal', string_number])
                string_number += 1
##                print(fp_anneal)
                tp_anneal = Seq(subitems['Sequence'][-40:], IUPAC.unambiguous_dna).reverse_complement()
                gibson_primer_temp_list.append([tp_anneal, key, counter, 'tp_anneal', string_number])
                string_number += 1
##                print(tp_anneal)
                fp_overhang = Seq(subitems['Sequence'][0:20], IUPAC.unambiguous_dna).reverse_complement()
                gibson_primer_temp_list.append([fp_overhang, key, counter, 'fp_overhang', string_number])
                string_number += 1
##                print(fp_overhang)
                tp_overhang = Seq(subitems['Sequence'][-20:], IUPAC.unambiguous_dna)
                gibson_primer_temp_list.append([tp_overhang, key, counter, 'tp_overhang', string_number])
                string_number += 1
##                print(tp_overhang)
                counter += 1
    ##the case where Gibson30 is the strategy##
    ##---------------------------------------##
            elif subitems['Strategy'] == 'Gibson30' and subitems['Construct'] == key:
##                print(key)
                fp_anneal = Seq(subitems['Sequence'][0:30], IUPAC.unambiguous_dna)
                gibson_primer_temp_list.append([fp_anneal, key, counter, 'fp_anneal', string_number])
                string_number += 1
##                print(fp_anneal)
                tp_anneal = Seq(subitems['Sequence'][-30:], IUPAC.unambiguous_dna).reverse_complement()
                gibson_primer_temp_list.append([tp_anneal, key, counter, 'tp_anneal', string_number])
                string_number += 1
##                print(tp_anneal)
                fp_overhang = Seq(subitems['Sequence'][0:30], IUPAC.unambiguous_dna).reverse_complement()
                gibson_primer_temp_list.append([fp_overhang, key, counter, 'fp_overhang', string_number])
                string_number += 1
##                print(fp_overhang)
                tp_overhang = Seq(subitems['Sequence'][-30:], IUPAC.unambiguous_dna)
                gibson_primer_temp_list.append([tp_overhang, key, counter, 'tp_overhang', string_number])
                string_number += 1
##                print(tp_overhang)
                counter += 1
    ##the case where iPCR is the strategy##
    ##-----------------------------------##
        ##NOTE: here, I make the primers directly.
            if subitems['Strategy'] == 'iPCR' and subitems['Construct'] == key:
                fw_primer = Seq(subitems['Sequence'][0:60], IUPAC.unambiguous_dna)
                primer_list.append([fw_primer, subitems['Construct'], 1, 'fw primer'])
                re_primer = Seq(subitems['Sequence'][-60:], IUPAC.unambiguous_dna).reverse_complement()
                primer_list.append([re_primer, subitems['Construct'], 1, 're primer'])

##here, i process all the gibson primers to get the final list of primers##
##=======================================================================##
    construct_num = 1
    temp = []
    part_num = 1
    temp_row_num = 1
    max_seq_num = 0

    for row in gibson_primer_temp_list:
##        print row

        max_seq_num = 0

        for x in gibson_primer_temp_list:
            if int(x[1]) > construct_num:
                pass
            if int(x[1]) == construct_num:
                max_seq_num += 1
##        print('Const. number counter is at ' + str(construct_num) + ' and current maximum known number of sequences is ' + str(max_seq_num))

##        print(row[1])

##        if int(row[1]) < construct_num:
##            while construct_num < int(row[1]):
##        print(max_seq_num)
##        for row in gibson_primer_temp_list:
##            if int(row[1]) == construct_num:
##                max_seq_num += 1
##            if int(row[1]) > construct_num:
##                break

        #print('Construct number is ' + str(row[1]) + ' and seq. number is ' + str(row[4]))
        #print('Const. number counter is ' + str(construct_num) + ' and max. seq. number is ' + str(max_seq_num) + '.')

        if int(row[1]) > construct_num:
            part_num = 1
            while construct_num < int(row[1]):
                #print('Construct number is ' + str(construct_num))
                construct_num += 1
##                temp_row_num += 1 #do not uncomment
            #continue - not to be added back again!

        if int(row[1]) == construct_num:

            if int(row[4]) == max_seq_num:

##                print(row)
                temp.append(row)
                temp_row_num += 1
##                print('We are going to make primers that join the first and last part in construct ' + str(construct_num))
##                print('Grabbing overhang portion from part ' + str(part_num) + ', which is sequence ' + str(row[4]) + '. It has the sequence ' + str(row[0]))

                for w in gibson_primer_temp_list: #grabs overhang from last part
                    if w[1] == row[1] and w[4] == max_seq_num and w[3] == 'tp_overhang':
                        overhang = w
                        break

##                print('Grabbing the first sequence...')
                for x in gibson_primer_temp_list:
##                    print(row[1] == x[1] and x[4] == 1)
                    if x[1] == row[1] and x[4] == 1:
##                        print(x)
                        anneal = x
##                        print('The first sequence is ' + str(anneal))
                        fw_primer = overhang[0] + anneal[0]
##                        print('The forward primer on the first part is: ' + str(fw_primer))
                        primer_list.append([fw_primer, construct_num, x[2], 'fw primer'])
                        break

##                print('Grabbing the third sequence...')
                for y in gibson_primer_temp_list:
##                    print(row[1] == y[1] and y[4] == 3)
                    if row[1] == y[1] and y[4] == 3 and y[3] == 'fp_overhang':
##                        print(y[0])
                        overhang = y
##                        print('The third sequence is ' + str(overhang))
                        break

##                print('Grabbing the (n-2)th sequence...')
                steps_backward = 2
                target_seq_num = max_seq_num - steps_backward
                for z in gibson_primer_temp_list:
##                    print(row[1] == z[1] and z[4] == target_seq_num)
                    if row[1] == z[1] and z[4] == target_seq_num:
##                        print(z[0])
                        anneal = z
##                        print('The n-2th sequence is ' + str(anneal))
                        break

                re_primer = overhang[0] + anneal[0]
                primer_list.append([re_primer, construct_num, z[2], 're primer'])
                continue

            if part_num == int(row[2]) and part_num == 1: #if the part number counter = part number and is the first part number
                #print(row)
                temp.append(row)
                temp_row_num += 1
                continue #do NOT delete this continue

            if part_num < int(row[2]):
                #print('Current part is: ' + str(part_num) + '. Upping part number.' + '\n')
                part_num += 1
                #do NOT add in a "continue" here
                                
##            print('Current part number is...')
##            print(part_num)
##            print(row)

            if part_num == int(row[2]) and row[3] == 'tp_overhang':
##                print('Appending tp_overhang from part ' + str(row[2]))
##                print(row)
                temp.append(row)
                temp_row_num += 1
                part_num += 1
                continue


            if part_num == int(row[2]) and row[3] == 'fp_anneal':
##                print('\n')
##                print(row)
                temp.append(row)
                temp_row_num += 1
##                print('Current part is: ' + str(part_num))
##                print('Grabbing tp_overhang from part ' + str(part_num - 1) + '...')
##                print('Currnt construct number is ' + str(construct_num))
                for y in gibson_primer_temp_list:
##                    print(y)
                    if int(y[1]) == construct_num and y[2] == (row[2] - 1) and y[3] == 'tp_overhang':
                        prev_tp_overhang = y
##                        print prev_tp_overhang
##                print('Sequence of tp_overhang from part ' + str(part_num - 1) + ' is: ' + prev_tp_overhang[0])
                fw_primer_current = prev_tp_overhang[0] + row[0]
##                print('Appending to master primer list...')
                primer_list.append([fw_primer_current, construct_num, part_num, 'fw primer'])
##                print('Forward primer is: ' + str(fw_primer_current) + '\n')
                continue

            if part_num == int(row[2]) and row[3] == 'tp_anneal':
##                print(row)
                temp.append(row)
                temp_row_num += 1
                continue


            if part_num == int(row[2]) and row[3] == 'fp_overhang':
##                print(row)
                temp.append(row)
                temp_row_num += 1
##                print('Current temp_row_num is ' + str(temp_row_num))
##                print('Current part is: ' + str(part_num))
##                print('Grabbing tp_anneal from part ' + str(part_num - 1) + '...')
##                x = 1
                for y in temp:
##                    print(y)
##                    print(int(y[1]) == construct_num and y[2] == (part_num - 1) and y[3] == 'tp_anneal')
                    if int(y[1]) == construct_num and int(y[2]) == (part_num - 1) and y[3] == 'tp_anneal':
                        prev_tp_anneal = y
##                        print(row)
##                        print(y)
                        pass
##                print('Sequence of tp_anneal from part ' + str(part_num - 1) + ' is: ' + prev_tp_anneal[0])
                re_primer_prev = row[0] + prev_tp_anneal[0]
##                print('Appending to master primer list...')
                primer_list.append([re_primer_prev, construct_num, part_num - 1, 're primer'])
##                print('Reverse primer for previous part is: ' + str(re_primer_prev) + '\n')
##                part_num += 1 #do not uncomment
                continue


            continue


##here, i create a list containing the names for each primer##
##==========================================================##
primers_names_list = []

with open('constructs-to-make-shortened2.csv', 'rU') as constructs:
    construct_list = csv.DictReader(constructs)

    part_num = 1
    construct_counter = 1
    part_counter = 1

    for row in construct_list:

##        print(construct_counter > int(row['Construct']))
        if construct_counter > int(row['Construct']):
            construct_counter = 1
            part_counter = 1
##        print(construct_counter < int(row['Construct']))
        if construct_counter < int(row['Construct']):
            construct_counter += 1
            part_counter = 1

        if construct_counter == int(row['Construct']):
            construct_number = row['Construct']
            source = row['Source']
            content = row['Content']
            strategy = row['Strategy']
            fw_primer_description = str('Fw ' + strategy + ' primer on ' + source + ' to extract ' + content)
            re_primer_description = str('Re ' + strategy + ' primer on ' + source + ' to extract ' + content)
##            print(fw_primer_description)
##            print(re_primer_description)
            primers_names_list.append({'notes':fw_primer_description,
                                       'construct number':construct_number,
                                       'direction':'fw primer',
                                       'part number':part_counter,
                                       'source':source})
            primers_names_list.append({'notes':re_primer_description,
                                       'construct number':construct_number,
                                       'direction':'re primer',
                                       'part number':part_counter,
                                       'source':source})
            part_counter += 1
            continue

##for row in primers_names_list:
##    print(row)

##here i print out all the primers generated##
##==========================================##
##    print('\n')
##    for row in primer_list:
##        print(row)

##here, i write all the contents of the primers list to a csv file##
##================================================================##

formatted_primer_list = []
for row in primer_list:
    formatted_primer_list.append([str(row[0]), row[1], row[2], row[3]])

primers_output = open('primers-list.csv', 'w+')
wr = csv.writer(primers_output, dialect='excel')
primers_output.write(str('Primer Sequence,Construct Number,Part Number,Direction'))
primers_output.write('\n')
for row in formatted_primer_list:
    wr.writerow(row)

primers_output.close()



##here, i sort all the primers##
##============================##
primers_list = open('primers-list.csv', 'rU')
primers_unsorted = csv.DictReader(primers_list)
primers_sorted = sorted(primers_unsorted, key=operator.itemgetter('Construct Number', 'Part Number'))

primers_output = open('primers-sorted.csv', 'w+')
fields = ['Part Number', 'Construct Number', 'Direction', 'Primer Sequence']
wr = csv.DictWriter(primers_output, delimiter=',', fieldnames=fields)
primers_output.write(str('Part Number,Construct Number,Direction,Primer Sequence'))
primers_output.write('\n')
for row in primers_sorted:
    wr.writerow(row)

primers_output.close()

##here, i attach the primer names to the primer sequences##
##=======================================================##
    #identifying information is the construct number, part number AND direction

with open('primers-sorted.csv', 'rU') as f:
    primers_without_names_dict = csv.DictReader(f)

    primers_without_names = []

    for row in primers_without_names_dict:
        primers_without_names.append(row)

    
##    for row in primers_names_list:
##        print row

##    for row in primers_without_names:
##        print row

    primers_with_notes = []


    for row in primers_names_list: #recall that primers_names_list is a list of dictionaries
        for x in primers_without_names:
            if (
                int(x['Part Number']) == row['part number'] and
                x['Construct Number'] == row['construct number'] and
                x['Direction'] == row['direction']
            ):
                primers_with_notes.append(
                    {
                        'part number': row['part number'], 
                        'construct number': row['construct number'], 
                        'notes': row['notes'], 
                        'primer sequence':x['Primer Sequence'],
                        'length':len(x['Primer Sequence']),
                        'direction':x['Direction']
                    }
                )

                break
                # If you are only expecting one match from the primers_without_names
                # collection, or wish to enforce that, you can add a break statement after
                # the insertion here to break out of the inner comparison loop and move on
                # to the next row item

##    for p in primers_with_notes:
##        print p
##
##    print
##    print len(primers_with_notes)

with open('primers-with-notes.csv', 'w+') as primers_output:

    fields = ['construct number','part number','length','notes','primer sequence','direction']
    wr = csv.DictWriter(primers_output, delimiter = ',', fieldnames = fields)
    primers_output.write('construct number,part number,length,notes,primer sequence,direction')
    primers_output.write('\n')

    for row in primers_with_notes:
        wr.writerow(row)

##here, i try to access to Google Docs spreadsheet##
##================================================##
email = 'youremailaddress'
password = 'yourpassword'
# Find this value in the url with 'key=XXX' and copy XXX below
##spreadsheet_key = 'get_spreadsheet_key'
# All spreadsheets have worksheets. I think worksheet #1 by default always
# has a value of 'gid=0'
##worksheet_id = 'od6'

spr_client = gdata.spreadsheet.service.SpreadsheetsService()
spr_client.email = email
spr_client.password = password
spr_client.source = 'Test Python Spreadsheet'
spr_client.ProgrammaticLogin()


##from the csv file primers-with-notes.csv, get the entire csv file as a list of dictionary items##
##===============================================================================================##

f = open('primers-with-notes.csv', 'rU')
primers_with_notes_dict = csv.DictReader(f)
primers_with_notes = []

for row in primers_with_notes_dict:
    primers_with_notes.append(row)

primers_id = 'od7' #get the spreadsheet id for the primers list. first is od6, second is od7, third is od8 etc...

doc_name = 'Test Python Spreadsheet'

gd_client = gdata.spreadsheet.service.SpreadsheetsService()
gd_client.email = email
gd_client.password = password
gd_client.source = 'Test Python Spreadsheet'
gd_client.ProgrammaticLogin()


q = gdata.spreadsheet.service.DocumentQuery()
q['title'] = doc_name
q['title-exact'] = 'true'
feed = gd_client.GetSpreadsheetsFeed(query=q)
spreadsheet_id = feed.entry[0].id.text.rsplit('/',1)[1]
feed = gd_client.GetWorksheetsFeed(spreadsheet_id)
worksheet_id = feed.entry[0].id.text.rsplit('/',1)[1]

current_primers_online = gd_client.GetListFeed(spreadsheet_id, primers_id).entry

current_primers_copy = [] #a list of dictionary items

#find out the highest primer number#
#----------------------------------#
primers_counter = 1
for row in current_primers_online:
    
    primers_counter += 1

primers_with_notes_names = []

for row in primers_with_notes:
    
    #Write new primer to the primers worksheet on the test spreadsheet while also checking to make sure it's not already there#
    #-------------------------------------------------------------------------------------------------------------------------#
    primers_dict = {}
    primers_dict['name'] = 'EMP' + str(primers_counter)
    primers_dict['sequence'] = row['primer sequence']
    primers_dict['notes'] = row['notes']
    primers_dict['length'] = row['length']
##    primers_dict['construct number'] = row['construct number']
##    primers_dict['part number'] = row['part number']
##    print(primers_dict)

    primers_dict2 = {}
    primers_dict2['name'] = 'EMP' + str(primers_counter)
    primers_dict2['sequence'] = row['primer sequence']
    primers_dict2['notes'] = row['notes']
    primers_dict2['length'] = row['length']
    primers_dict2['construct number'] = row['construct number']
    primers_dict2['part number'] = row['part number']
    primers_dict2['direction'] = row['direction']

    for x in current_primers_online:
        length = x.custom['length'].text
        notes = x.custom['notes'].text
        name = x.custom['name'].text
        sequence = x.custom['sequence'].text

        current_primers_copy.append({'length':length,
                                     'notes':notes,
                                     'name':name,
                                     'sequence':sequence})
    primer_exists_already = any(row['sequence'] == primers_dict['sequence'] for row in current_primers_copy)
    
    if not primer_exists_already:
        current_primers_copy.append(primers_dict) #appends the list to current_primers_copy, so that duplicates don't get added into the list.
        new_entry = spr_client.InsertRow(primers_dict, spreadsheet_id, primers_id)
        
        primers_with_notes_names.append(primers_dict2)
        
        primers_counter += 1
##        if isinstance(new_entry, gdata.spreadsheet.SpreadsheetsList):
##          print "Insert row succeeded."
##          
##        else:
##          print "Insert row failed."
          
    if primer_exists_already: #print: the primer already in list; add the existing primer from the 
##        print "The primer is already in the list."
        for row in current_primers_copy:
            if row['sequence'] == primers_dict['sequence']:
                primers_existing_dict = {}
                primers_existing_dict['name'] = row['name']
                primers_existing_dict['sequence'] = row['sequence']
                primers_existing_dict['notes'] = row['notes']
                primers_existing_dict['length'] = row['length']
                primers_existing_dict['construct number'] = primers_dict2['construct number']
                primers_existing_dict['part number'] = primers_dict2['part number']
                primers_existing_dict['direction'] = primers_dict2['direction']
                primers_with_notes_names.append(primers_existing_dict)
                break

##write primers_with_notes_names to a csv file##
##============================================##
with open('primers-with-notes-names.csv', 'w+') as primers_output:

    fields = ['name','construct number','part number','length','notes','sequence','direction']
    wr = csv.DictWriter(primers_output, delimiter = ',', fieldnames = fields)
    primers_output.write('name,construct number,part number,length,notes,primer sequence,direction')
    primers_output.write('\n')

    for row in primers_with_notes_names:
        wr.writerow(row)

##now, assign primers to construct and part numbers, spit out master spreadsheet containing the following headers:
##Construct	Strategy	Construct Name	Sequence	Source	Content	Fw 	Fw sequence	Description	Length (bp)	Re	Re sequence	Description	Length (bp)	PCR (bp)	#	PCR ###
##===========================================================================================================================================================================================================================##
master_list = []

construct_counter = 1
part_number = 1


##f = open('constructs-to-make-shortened2.csv', 'rU')
##reader = csv.DictReader(f)
##for row in reader:
##    print(row)
##
##f.close()

##f = open('primers-with-notes-names.csv', 'rU')
##reader = csv.DictReader(f)
##for row in reader:
##    print(row)

##f.close()

with open('constructs-to-make-shortened2.csv', 'rU') as constructs:
    construct_list = csv.DictReader(constructs)

    with open('primers-with-notes-names.csv', 'rU') as primers:
        primers_list = list(csv.DictReader(primers))

        #make list of constructs for checking later on#
##        construct_numbers_list = []
##        for row in primers_list:
##            construct_numbers_list.append(row['construct number'])
##
##        print(construct_numbers_list)
        

        for construct in construct_list:
##            print('Currently at construct number ' + construct['Construct'])
##            print('Construct counter at ' + str(construct_counter))
##            print('Part number counter is at ' + str(part_number))
            master_row = {}
            master_row['construct'] = construct['Construct']
            master_row['strategy'] = construct['Strategy']
            master_row['construct name'] = construct['Construct Name']
            master_row['sequence'] = construct['Sequence']
            master_row['source'] = construct['Source']
            master_row['content'] = construct['Content']


##            print('We are at construct number ' + str(construct['Construct']))
##            print('Construct counter is at ' + str(construct_counter))
            is_next_construct = (int(construct['Construct']) > construct_counter)
##            print('Are we at the next construct?')
##            print(is_next_construct)
            
            if is_next_construct:
                part_number = 1
                construct_counter = int(construct['Construct'])
##            print('Part number is now ' + str(part_number))

            for primer in primers_list:
##                print(primer)


##                    print('Is primer ' + str(primer['name']) + ' associated with the construct?')
                is_associated_with_construct = bool(primer['construct number'] == construct['Construct'] and str(primer['part number']) == str(part_number))
##                    print(is_associated_with_construct)
                if(is_associated_with_construct == False):
                    continue

                is_forward = bool(primer['construct number'] == construct['Construct'] and str(primer['part number']) == str(part_number) and primer['direction'] == 'fw primer')
                
##                print('Primer ' + str(primer['name']) + ' is a forward primer?')
##                print(is_forward)

                is_reverse = bool(primer['construct number'] == construct['Construct'] and str(primer['part number']) == str(part_number) and primer['direction'] == 're primer')

##                print('Primer ' + str(primer['name']) + ' is a reverse primer?')
##                print(is_reverse)
                
                if is_forward:
                    master_row['primer1'] = primer['name']
                    master_row['primer1 sequence'] = primer['primer sequence']
                    master_row['primer1 description'] = primer['notes']
                    master_row['primer1 length'] = primer['length']
                    master_row['primer1 part number'] = primer['part number']
##                        print(master_row)
                    continue

                elif is_reverse:
                    master_row['primer2'] = primer['name']
                    master_row['primer2 sequence'] = primer['primer sequence']
                    master_row['primer2 description'] = primer['notes']
                    master_row['primer2 length'] = primer['length']
                    master_row['primer2 part number'] = primer['part number']
##                        print(master_row)
                    part_number += 1
##                    print('Part number now = ' + str(part_number) + '\n')
                    master_list.append(master_row)
                    break
##            break
##                master_list.append(master_row)

                    

##for row in master_list:
##    print(row)
##

#write to csv flie#
with open('construction-master-list.csv', 'w+') as output:

    fields = ['construct',
              'strategy',
              'construct name',
              'sequence',
              'source',
              'content',
              'primer1',
              'primer1 sequence',
              'primer1 description',
              'primer1 length',
              'primer1 part number',
              'primer2',
              'primer2 sequence',
              'primer2 description',
              'primer2 length',
              'primer2 part number']
    wr = csv.DictWriter(output, delimiter = ',', fieldnames = fields)
    output.write('construct,strategy,construct name,sequence,source,content,primer1,primer1 sequence,primer1 description,primer1 length,primer1 part number,primer2,primer2 sequence,primer2 description,primer2 length,primer2 part number')
    output.write('\n')
    
    for row in master_list:
        wr.writerow(row)
