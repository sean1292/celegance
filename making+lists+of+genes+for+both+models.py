
# coding: utf-8

# In[30]:

import xml.etree.ElementTree as ET
from string import punctuation
import re 

kaletaModel = ET.parse('C:\Users\user\Documents\Lewis_LAB\SBML files\Celegans_ecoli.xml').getroot()
icellMod = ET.parse('C:\Users\user\Documents\Lewis_LAB\SBML files\iCEL1273.xml').getroot()

root = kaletaModel.getchildren()
kaletaListOfRxns = root[0].getchildren()[2]
icellListOfRxns = icellMod.getchildren()[0].getchildren()[4]

genesK = []
genesI = []
iGeneFormat = re.compile('(?:(?:^GENE_LIST: )|(?:; ))([^\s;]+)')
kGeneFormat = re.compile('(?:GENE ASSOCIATION: \(| or \(|GENE_ASSOCIATION: \(| and | and \(| or |GENE ASSOCIATION: |GENE_ASSOCIATION: )((\w+\.\w+)|(\w+-\w+\.\w+)|(\w+.\w+))')

for rxnI in icellListOfRxns:
    if len(rxnI)>1:
        rxnNotesI = rxnI[0]
        for x in range(len(rxnNotesI)):
            geneList = re.findall(iGeneFormat,str(rxnNotesI[x].text))
            if geneList:
                break
        if len(geneList)>0:
            for gene in geneList:
                
                if gene not in genesI:
                    genesI.append(gene)
                    
iGenesTxt = open('iGenesList.txt', 'w')
for gene in genesI:
    iGenesTxt.write(str(gene) +'\n')
                    
for rxnK in kaletaListOfRxns:
    rxnNotesK = rxnK[0][0]
    length =  len(rxnNotesK)
    for x in range(length):
        string = str(rxnNotesK[x].text)
        if re.search(r'\bgene association\b', string, re.I) or re.search(r'\bgene_association\b', string, re.I):
            #print string
            geneListK = re.findall(kGeneFormat,string)
            if geneListK:
                break
    if len(geneListK)>0:
        for gene in geneListK:
            #print str(gene)    
            
            if gene not in genesK:
                genesK.append(gene)   

kGenesTxt = open('kGenesList.txt', 'w')
for gene in genesK:
    
    kGenesTxt.write(str(gene[0]) +'\n')


# In[29]:

print len(genesI)
print len(genesK)


# In[15]:

rxnK = kaletaListOfRxns[0]
print len(rxnK[0][0])


# In[58]:

kGenesTxt = open('kGenesList.txt', 'w')
genesConverted = []
with open('kGenesListConverted.txt') as f:
    for line in f:
        if re.search(r'(Cel\.\d+)', line):          
            geneId = re.search(r'(Cel\.\d+)', line)
            #rint str(geneId.group())
            if str(geneId.group()) not in genesConverted:
                genesConverted.append(str(geneId.group()))
for gene in genesConverted:
    kGenesTxt.write(gene + '\n')


# In[45]:

print len(genesConverted)


# In[57]:

iGenesTxt = open('iGenesList.txt', 'w')
genesConverted = []
with open('iGenesListConverted.txt') as f:
    for line in f:
        if re.search(r'(Cel\.\d+)', line):          
            geneId = re.search(r'(Cel\.\d+)', line)
            #rint str(geneId.group())
            if str(geneId.group()) not in genesConverted:
                genesConverted.append(str(geneId.group()))
for gene in genesConverted:
    iGenesTxt.write(gene + '\n')


# In[48]:

print len(genesConverted)


# In[51]:

iGenesTxt = open('IcellNotKaletaUNIGENE.txt', 'w')
genesConverted = []
with open('IcellNotKaletaUNIGENE_RAW.txt') as f:
    for line in f:
        if re.search(r'(Cel\.\d+)', line):          
            geneId = re.search(r'(Cel\.\d+)', line)
            if str(geneId.group()) not in genesConverted:
                genesConverted.append(str(geneId.group()))
for gene in genesConverted:
    iGenesTxt.write(gene + '\n')


# In[54]:

kGenesTxt = open('KaletaNotIcellUNIGENE.txt', 'w')
genesConverted = []
with open('KaletaNotIcellUNIGENE_RAW.txt') as f:
    for line in f:
        if re.search(r'(Cel\.\d+)', line):          
            geneId = re.search(r'(Cel\.\d+)', line)
            if str(geneId.group()) not in genesConverted:
                genesConverted.append(str(geneId.group()))
for gene in genesConverted:
    kGenesTxt.write(gene + '\n')


# In[ ]:



