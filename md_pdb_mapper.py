import os, sys, gzip
from pymol import cmd

for line in gzip.open('../covid19/data/interactions/find_interactions/P0DTC2_mdi.txt.gz', 'rt'):
    if line[0] != '#' and line.split('\t')[-1].replace('\n', '') == '3DID':
        motif = line.split('\t')[2]
        domain = line.split('\t')[9]
        startm = line.split('\t')[34]
        endm = line.split('\t')[35]
        flag = 1
        for line2 in open('../covid19/scripts/temp_blastp/P0DTC2+'+domain+'+'+motif+'_blastp.txt', 'r'):
            if 'No hits found' in line2:
                flag = 0
                break
            elif line2[0] == '>':
                print (motif, domain, startm, endm)
                print (line2)
                pdb = line2.split('>')[1].split('_')[0]
                chain = line2.split('_')[1]
                print (chain)
                break
        if flag == 0:
            continue
        else:
            break

#pdb = '6mnl'
cmd.load('../covid19/data/pdbs/'+pdb.upper()+'.pdb')
cmd.hide('everything', 'all')
cmd.show('cartoon', 'chain '+chain)
#cmd.color('red', 'chain A')
#cmd.color('yellow', 'chain B')
#cmd.color('blue', 'resi 7-10')
#cmd.label('i. 7 and n. CA', "'L'")
