#!/usr/bin/env python3

## pdb_mapper.py takes a CoV2 protein accession along with its motif's match, start and end positions, and interacting domain and its source as an input
## It then superimpoees the given peptide (motif'match) onto a chain along with its interacting domain onto another chain within same complex

## Import built-in libraries
import os, sys, gzip
from pymol import cmd
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB import *
import urllib.request

## Header
__author__ = "Gurdeep Singh"
__license__ = "GPL"
__email__ = "gurdeep330@gmail.com"

def download_pdb_cif(pdb):
    if os.path.isfile('pdbs/'+pdb+'.cif') == False:
        os.system('wget https://files.rcsb.org/view/'+pdb+'.pdb -O '+'pdbs/'+pdb+'.pdb')
        os.system('wget https://files.rcsb.org/view/'+pdb+'.cif -O '+'pdbs/'+pdb+'.cif')

def main(d1, d2):
    ## Make an output folder
    outdir = 'output_dd/'+input_protein+'+'+d1+'+'+d2+'/'
    #print (outdir)
    if os.path.isdir(outdir) == False:
        os.system('mkdir '+outdir)

    domains = d1 + '+' + d2
    for pdb in dd[domains]:
        download_pdb_cif(pdb)
        for coordinates in dd[domains][pdb]:
            #print (pdb, domains, dd[domains][pdb][coordinates])
            chainA = coordinates.split('+')[0].split(':')[0]
            chainA_res = coordinates.split('+')[0].split(':')[1]
            chainB = coordinates.split('+')[1].split(':')[0]
            chainB_res = coordinates.split('+')[1].split(':')[1]
            cmd.load('pdbs/'+pdb.upper()+'.pdb')
            cmd.hide('everything')
            cmd.show('cartoon', 'chain '+chainA+'+'+chainB)
            cmd.set('cartoon_fancy_helices', 1)
            print (pdb, chainA, chainA_res, chainB, chainB_res)
            try:
                x = cmd.centerofmass('chain '+chainA)
                print ('Success chainA')
                cmd.pseudoatom('chainA_label', pos=x)
                global nameA
                nameA = id_to_name[d1] + '(' + d1+')'
                cmd.label('chainA_label', 'nameA')
            except:
                print ('Failed chainA')
            try:
                x = cmd.centerofmass('chain '+chainB)
                print ('Success chainB')
                cmd.pseudoatom('chainB_label', pos=x)
                global nameB
                nameB = id_to_name[d2] + '(' + d2+')'
                cmd.label('chainB_label', 'nameB')
                x = cmd.centerofmass('chain '+chainA+'+'+chainB)
                cmd.origin(position=x)
                cmd.center('origin')
            except:
                print ('Failed chainB')
            cmd.set('label_size', 7.5)
            cmd.set('cartoon_fancy_helices', 1)
            cmd.color('red', 'chain '+chainA)
            cmd.color('orange', 'chain '+chainB)
            for row in dd[domains][pdb][coordinates]:
                res1 = row[2]
                res2 = row[3]
                #cmd.distance('chain '+chainA+' and i. '+res1+' and n. CB', 'chain '+chainB+' and i. '+res2+' and n. CB')
                cutoff = 6.5
                m = cmd.distance('dist', 'chain '+chainA+' and i. '+res1+' and n. CB', 'chain '+chainB+' and i. '+res2+' and n. CB', cutoff, 0)
                #m = cmd.get_distance(atom1='chain '+chainA+' and i. '+res1+' and n. CB',atom2='chain '+chainB+' and i. '+res2+' and n. CB',state=0)
                #print (pdb, m, chainA, res1, chainB, res2)
                cmd.select('res1', 'chain '+chainA+' and resi '+res1)
                cmd.select('res2', 'chain '+chainB+' and resi '+res2)

                if float(m) != 0.0:
                    cmd.show('sticks', 'res1')
                    cmd.show('sticks', 'res2')
                    cmd.color('cyan', 'res1')
                    cmd.color('yellow', 'res2')

            cmd.save(outdir+pdb+'_'+input_protein+'_'+d1+'_'+d2+'_'+coordinates.replace(':', '_').replace('+', '_').replace('-', '_')+'.pse')
            cmd.save(outdir+pdb+'_'+input_protein+'_'+d1+'_'+d2+'_'+coordinates.replace(':', '_').replace('+', '_').replace('-', '_')+'.pdb')
            cmd.delete('all')
        #break
    '''
    sys.exit()

    ## Pymol starts from here
    for pdb in pdbs:
        if pdbs[pdb]['motif_chains'] != [] and pdbs[pdb]['domain_chains'] != '':
            #pdb = '6mnl'
            print (pdb)
            #print (pdbs)
            print (pdbs[pdb]['domain_chains'].split(';'))
            #sys.exit()
            print (pdbs[pdb]['motif_chains'][0])
            print (pdbs[pdb]['motif_chains'])
            print (pdbs[pdb]['domain_chains'])
            mot_chain = pdbs[pdb]['motif_chains'][0].split(':')[0]
            mot_chain_res = pdbs[pdb]['motif_chains'][0].split(':')[1]
            cmd.load('pdbs/'+pdb.upper()+'.pdb')
            cmd.fab('Z/'+str(startm)+'/ '+match, 'input_protein_motif')
            cmd.hide('everything')
            cmd.show('cartoon', 'chain '+str(mot_chain))
            cmd.show('cartoon', 'input_protein_motif')
            cmd.color('red', 'chain '+mot_chain)
            cmd.color('yellow', 'input_protein_motif')
            cmd.select('struc_motif', 'resi '+mot_chain_res+' and chain '+mot_chain)
            cmd.color('cyan', 'struc_motif')
            cmd.super('input_protein_motif', 'struc_motif')
            cmd.center('struc_motif')
            cmd.label('mot, i. '+str(mot_chain_res.split('-')[0])+' and n. CA and chain '+mot_chain, 'motif')
            cmd.set('label_color', 'cyan', 'struc_motif')
            for dom in pdbs[pdb]['domain_chains'].split(';')[:-1]:
                dom_chain = dom.split(':')[0]
                dom_chain_res = dom.split(':')[1]
                cmd.show('cartoon', 'chain '+str(dom_chain))
                cmd.color('grey', 'chain '+str(dom_chain))
                cmd.color('orange', 'resi '+dom_chain_res+' and chain '+dom_chain)
                cmd.label('dom, i. '+str(dom_chain_res.split('-')[0])+' and n. CA and chain '+dom_chain, 'domain')
            cmd.save(outdir+pdb+'_'+input_protein+'_'+motif+'_'+domain+'_'+str(startm)+'_'+str(endm)+'_'+source+'.pse')
            cmd.save(outdir+pdb+'_'+input_protein+'_'+motif+'_'+domain+'_'+str(startm)+'_'+str(endm)+'_'+source+'.pdb')
            cmd.delete('all')
            #break
    '''

def read_3did_ddi():
    dd = {}
    id_to_name = {}
    for line in gzip.open('/home/gurdeep/projects/DB/3did/3did_flat.gz', 'rt'):
        if line.split()[0] == '#=ID':
            domain1 = line.split()[3].split('@')[0].replace('(', '').replace(')', '').split('.')[0]
            domain2 = line.split()[4].split('@')[0].replace('(', '').replace(')', '').split('.')[0]
            id_to_name[domain1] = line.split()[1]
            id_to_name[domain2] = line.split()[2]
            domains = domain1+'+'+domain2
            if domains not in dd:
                dd[domains] = {}
        elif line.split()[0] == '#=3D':
            pdb = line.split()[1].upper()
            coordinates = line.split()[2]+'+'+line.split()[3]
            if pdb not in dd[domains]:
                dd[domains][pdb] = {}
            dd[domains][pdb][coordinates] = []
        elif line[:2] != '//':
            dd[domains][pdb][coordinates].append(line.split()[:4])
    return dd, id_to_name

dd, id_to_name = read_3did_ddi()
#print (dd)
print ('starting...')
nameA=''; nameB=''
done = []
for files in os.listdir('/home/gurdeep/projects/covid19/data/interactions/find_interactions/'):
    if files.endswith('_ddi.txt.gz'):
        for line in gzip.open('/home/gurdeep/projects/covid19/data/interactions/find_interactions/'+files, 'rt'):
            if line[0] != '#':
                #print (line.split('\t'))
                input_protein = line.split('\t')[0]
                domain1 = line.split('\t')[3]
                domain2 = line.split('\t')[7]
                domains = domain1+'+'+domain2
                if domain1+'+'+domain2 in dd:
                    #print (domain1+domain2, dd[domain1+'+'+domain2])
                    if (input_protein+'+'+domain1+'+'+domain2) not in done:
                        main(domain1, domain2)
                        done.append(input_protein+'+'+domain1+'+'+domain2)
                        #sys.exit()
                elif domain2+'+'+domain1 in dd:
                    #print (domain2+domain1, dd[domain2+'+'+domain1])
                    if (input_protein+'+'+domain2+'+'+domain1) not in done:
                        main(domain2, domain1)
                        done.append(input_protein+'+'+domain2+'+'+domain1)
                        #sys.exit()
print ('done')
