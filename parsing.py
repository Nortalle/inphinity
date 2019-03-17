import xml.etree.ElementTree as etree


primary_call = []
secondary_call = []
proteins = []



def parseXML(xml):
  
  """
  Parse XML files of the website : http://locate.imb.uq.edu.au

  Those files contains proteins and we extract those different informations :
    accn
    protein_sequence
    transcript_sequence
    tier1
    tier2

  It occurs that literature, experimental_data and externalannot tags could 
  contain tier1 and tiers2, but only the literature one seems to be reliable

  :param xml: xml file

  :type xml: xml file
  """

  # TODO a class would be better
  tier1 = {'design': None, 'goid': None}
  tier2 = {'design': None, 'goid': None}

  # To know if we are in the "literature" tag
  inLiterature = False
  inProtein = False

  # This is all the element we need to fill the database
  current_protein = {'acc': None, 'seq_prot': None, 'seq_adn': None, 'tier1' : None, 'tier2': None, 'goid': None}
  
  # We Browse in the xml document
  for event, elem in etree.iterparse(xml, events=('start', 'end')):
    
    # if it's a start tag
    if event == 'start':  

      if elem.tag == 'protein_sequence':
        current_protein['seq_prot'] = elem.text

      if elem.tag == 'transcript_sequence':
        current_protein['seq_adn'] = elem.text

      if elem.tag == 'literature':
        inLiterature = True

      if elem.tag == 'protein':
        inProtein = True

      if inLiterature:
        if elem.tag == 'tier1' and current_protein['tier1'] == None:
          current_protein['tier1'] = elem.text
        if elem.tag == 'tier2' and current_protein['tier2'] == None:
          current_protein['tier2'] = elem.text
        if elem.tag == 'location' and current_protein['goid'] == None:
          current_protein['goid'] = elem.attrib['goid']

      if inProtein:  
        if elem.tag == 'accn':
          current_protein['acc'] = elem.text
          
    # if it's an end tag
    if event == 'end':
      if elem.tag == 'LOCATE_protein':

        # TODO look what is the correct way to test multiple conditions that can be readable
        hasAcc = current_protein['acc'] != None
        hasSequ_prot = current_protein['seq_prot'] != None
        hasSeq_adn = current_protein['seq_adn'] != None
        hasTier1 = current_protein['tier1'] != None
        hasGoid = current_protein['goid'] != None

        # if we have the basic information on that protein
        if hasAcc and hasSequ_prot and hasSeq_adn and hasTier1 and hasGoid:
          goids = current_protein['goid'].split(';')
          tier1['design'] = current_protein['tier1']
          tier1['goid'] = goids[0]
          primary_call.append(tier1)
          if len(goids) > 1 and current_protein['tier2'] != None: # Some goid looks like "smth;   "
            tier2['design'] = current_protein['tier2']
            tier2['goid'] = goids[1]
            secondary_call.append(tier2)
          proteins.append(current_protein)

        # reset 
        current_protein = {'acc': None, 'seq_prot': None, 'seq_adn': None, 'tier1' : None, 'tier2': None, 'goid': None}
        tier1 = {'design': None, 'goid': None}
        tier2 = {'design': None, 'goid': None}
        elem.clear() # IMPORTANT LINE, if commented ==> memory heap space

      if elem.tag == 'literature':
        inLiterature = False
        
      if elem.tag == 'protein':
        inProtein = False

parseXML("./dataSets/LOCATE_human_v6_20081121.xml")
parseXML("./dataSets/LOCATE_mouse_v6_20081121.xml")

#----------------------------------------------------------------------------------------------
# TODO find a better way to make a set out of an array of dict
primary_call = list({v['design']:v for v in primary_call}.values())
secondary_call = list({v['design']:v for v in secondary_call}.values())

#----------------------------------------------------------------------------------------------
# just a way to display all the tiers of that file
print('Tier1 (', len(primary_call), ') : ')
for elem in sorted(primary_call, key=lambda k: k['design']):
  nbr_tab = 6 - (len(elem['design'])/8)
  
  tab = ''
  for i in range(int(nbr_tab)):
    tab = tab + '\t'

  print(elem['design'], tab, elem['goid'])

print('------- ------------ ---------- ----------- ')
print('Tier1 (', len(secondary_call), ') : ')
for elem in sorted(secondary_call, key=lambda k: k['design']):
  nbr_tab = 6 - (len(elem['design'])/8)
  
  tab = ''
  for i in range(int(nbr_tab)):
    tab = tab + '\t'

  print(elem['design'], tab, elem['goid'])

print('Number of proteins : ', len(proteins))


#----------------------------------------------------------------------------------------------
# Giving all the tiers an id
for (id, tier1) in zip(range(len(primary_call)), primary_call):
  tier1['id'] = id
  # print(tier1)
  
for (id, tier2) in zip(range(len(secondary_call)), secondary_call):
  tier2['id'] = id
  # print(tier2)

#----------------------------------------------------------------------------------------------
# Preparing the proteins to be filled in the DB, with foreign key of the primary call
db_proteins = []
db_protein = {'acc': None, 'seq_prot': None, 'seq_adn': None, 'Fk_tier1' : None, 'Fk_tier2': None}

for protein in proteins:
  db_protein['acc'] = protein['acc']
  db_protein['seq_prot'] = protein['seq_prot']
  db_protein['seq_adn'] = protein['seq_adn']

  tier1 = protein['tier1']
  tier2 = protein['tier2']

  #Corresponding foreign key
  for call in primary_call:
    if call['design'] == tier1:
      db_protein['Fk_tier1'] = call['id']


  #Corresponding foreign key
  for call in secondary_call:
    if call['design'] == tier2:
      db_protein['Fk_tier2'] = call['id']

  db_proteins.append(db_protein)    
  db_protein = {'acc': None, 'seq_prot': None, 'seq_adn': None, 'Fk_tier1' : None, 'Fk_tier2': None}

#beta testing
for i in range(10):
  current = db_proteins[i]
  print(current['acc'], current['Fk_tier1'], current['Fk_tier2'])


