# find orthogroups and perform alignments
#this script takes a specific list of genes, finds the orthogroup for each gene, groups all the sequences into one fasta, then aligns them
#might need biopython installed? But I don't think so
# cursor.execute('DROP TABLE singleton_trees')
# import requests, json, time
#from ete3 import PhyloTree
from subprocess import Popen, PIPE
from Bio.Align.Applications._Muscle import MuscleCommandline
import os
import time
import sqlite3
db = sqlite3.connect('drosophilaDatabase_neop')
db.isolation_level = None
cursor = db.cursor()
realGroupDict = {}
cursor.execute('SELECT * FROM finalOrthogroups')
for a,b in cursor.fetchall():
    realGroupDict[a] = b.split(',')
sing_list = []
with open('dsuz_singletons.txt','r') as file:
    for line in file:
        line = line.strip('\n')
        sing_list.append(line)
count, ti = 0, time.time()
noTree, noTreeList = 0, []
# geneTreeDict, geneTreePrunedDict = {},{}
count = 0
#the list here is obviously D. suzukii singleton specific, I think any list of gene ids can be feed through the loop to perform alignments though
#you will need to use the Results_Aug27 orthogroups
print('Fetching orthogroups and performing alignments...')
for g in sing_list:
    #check if the current gene is not present in any orthogroups
    searchCmd3 = ['grep',g,'/neopTranslations_forOrthofinder/OrthoFinder/Results_Jun25/Orthogroups/Orthogroups_UnassignedGenes.tsv']
    searchP3 = Popen(searchCmd3,stdout = PIPE, stderr = PIPE)
    out3, err3 = searchP3.communicate()
    #finding originial orthogroup in Orthofinder output
    if out3 != bytes('',encoding='UTF-8'):
        orthoGroup = None
#         print(g,'unassigned to group')
        noTree += 1
        continue
    else:
        #genes that are not 'unassigned'
        t = '/home/zoe/neopTranslations_forOrthofinder/OrthoFinder/Results_Jun25/Orthogroups/Orthogroups.txt'
        searchCmd1 = ['grep', g,t]
        searchP1 = Popen(searchCmd1,stdout = PIPE, stderr = PIPE)
        out, err = searchP1.communicate()
        if out != bytes('',encoding='UTF-8'):
            orthoGroup = out.split(bytes(':','utf-8'))[0].decode('utf-8')
        # else:
            # print(g,'is a mystery cat')
            # print(out, err)
    #find actual orthoGroup (the finalOrthogroups table should have all the groups, not just the singletons) and make fasta for alignment
    try:
        realGroup = [key for key in realGroupDict if g in realGroupDict[key]][0] #CHECK!
    except IndexError:
        noTree += 1
        realGroup = None
        continue
    unfilteredFasta = '/home/zoe/neopTranslations_forOrthofinder/OrthoFinder/Results_Jun25/Orthogroup_Sequences/' + orthoGroup + '.fa'
    inFasta = 'neopTranslations_by_orthogroup/group'+str(realGroup)+'.fa'
    speciesList = realGroupDict[realGroup]
    with open(unfilteredFasta, 'r') as file, open(inFasta,'w') as out:
        incSeq = False
        for line in file:
            if line.startswith('>') and line.strip('\n').strip('>') in speciesList:
                incSeq = True
                out.write(line)
            elif line.startswith('>'):
                incSeq = False
            elif incSeq == True:
                out.write(line)
# run MUSCLE on the fasta file to create an aligned fasta, muscle will need to be installed
    outFasta = 'home/zoe/neopTranslations_alignments/group' + str(realGroup) + '_alignment.fa'
    cline = ['muscle','-in', inFasta, '-out', outFasta]
    p = Popen(cline, stdout=PIPE, stderr=PIPE)
    out, err = p.communicate()
    count += 1
    if count % 200 == 0:
        print(count, 'done', round(time.time()-ti, 2), 'seconds')
