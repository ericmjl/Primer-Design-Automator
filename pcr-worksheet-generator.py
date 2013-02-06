import csv

with open('construction-master-list.csv', 'rU') as f:
    construction_master_dict = csv.DictReader(f)

    construction_master_list = []
    for row in construction_master_dict:
        construction_master_list.append(row)

##for row in construction_master_list:
##    print(row)

pcr_worksheet = []

for row in construction_master_list:

##    print row
    pcr_worksheet_row = {}

    pcr_worksheet_row['construct'] = row['construct']
    pcr_worksheet_row['strategy'] = row['strategy']
    pcr_worksheet_row['construct name'] = row['construct name']
    pcr_worksheet_row['source'] = row['source']
    
    pcr_worksheet_row['primer1'] = row['primer1']
##    pcr_worksheet_row['primer1 sequence'] = row['primer1 sequence']
##    pcr_worksheet_row['primer1 length'] = row['primer1 length']
    
    pcr_worksheet_row['primer2'] = row['primer2']
##    pcr_worksheet_row['primer2 sequence'] = row['primer2 sequence']
##    pcr_worksheet_row['primer2 length'] = row['primer2 length']

    #in the following code, we determine the length of the PCR product.
    if row['strategy'] == 'Gibson':
        pcr_worksheet_row['pcr (bp)'] = len(row['sequence']) + 40

    elif row['strategy'] == 'Gibson30':
        pcr_worksheet_row['pcr (bp)'] = len(row['sequence']) + 60

    elif row['strategy'] == 'iPCR':
        pcr_worksheet_row['pcr (bp)'] = len(row['sequence'])


    pcr_worksheet_row['part number'] = row['primer1 part number']
    pcr_worksheet_row['pcr number'] = str(pcr_worksheet_row['construct']) + '.' + str(pcr_worksheet_row['part number'])
        
    pcr_worksheet.append(pcr_worksheet_row)

    

##for row in pcr_worksheet:
##    print row

with open('pcr-worksheet.csv', 'w+') as output:
    fields = ['construct',
              'part number',
              'strategy',
              'construct name',
              'source',
              'primer1',
              'primer2',
              'pcr (bp)',
              'pcr number']
    wr = csv.DictWriter(output, delimiter = ',', fieldnames = fields)
    output.write('construct,part number,strategy,construct name,source,primer1,primer2,pcr (bp),pcr number')
    output.write('\n')

    for x in pcr_worksheet:
        wr.writerow(x)
                                                                        
