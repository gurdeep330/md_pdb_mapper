#!/usr/bin/env python3

## pdb_mapper.py takes a CoV2 protein accession along with its motif's match, start and end positions, and interacting domain and its source as an input
## It then superimpoees the given peptide (motif'match) onto a chain along with its interacting domain onto another chain within same complex

## Import built-in libraries
import os, sys, gzip
sys.path.append('/home/gurdeep/pymol/')
#from pymol import cmd
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB import *
import urllib.request

## Header
__author__ = "Gurdeep Singh"
__license__ = "GPL"
__email__ = "gurdeep330@gmail.com"

import os, sys, gzip

'''
input_protein = 'P0DTC2'
motif = 'DOC_MAPK_MEF2A_6'
domain = 'PF00069'
startm = 1211
endm = 1220
source = 'ELM'
'''

def main(input_protein, motif, domain, startm, endm, source):
    startm = int(startm)
    endm = int(endm)

    ## Make an output folder
    outdir = 'output/'+input_protein+'+'+motif+'+'+domain+'/'
    print (outdir)
    if os.path.isdir(outdir) == False:
        os.system('mkdir '+outdir)

    fasta = ''
    for line in open('../covid19/data/fasta_ncbi/covid_fasta/'+input_protein+'.fasta', 'r'):
        if line[0] != '>':
            fasta += line.replace('\n', '')

    match = fasta[startm-1:endm]
    #print (match)
    #sys.ext()


    def get_FASTA_from_pdb(pdb, given_chain, map):
        if pdb not in map:
            map[pdb] = {}
        map[pdb][given_chain] = {}
        fasta = ''
        num = 1

        try:
            parser = MMCIFParser()
            structure = parser.get_structure(pdb, 'pdbs/'+pdb+'.cif')
        except:
            parser = PDBParser()
            structure = parser.get_structure(pdb, 'pdbs/'+pdb+'.pdb')

        flag = 0
        for model in structure:
            for chain in model:
                if (chain.id == given_chain):
                    fasta += '>'+pdb+'_'+given_chain+'\n'
                    for residue in chain:
                        for atom in residue:
                            if (atom.id == 'CA'):
                                if residue.id[0] == ' ':
                                    #print (pdb, chain.id, residue.get_resname(), residue.id)
                                    AA = protein_letters_3to1[residue.get_resname()]
                                    fasta += AA
                                    map[pdb][given_chain][num] = residue.id[1]
                                    num += 1
                    fasta += '\n'
                    flag = 1
                    break
            if flag == 1:
                break
        return fasta, map

    def run_and_read_elm_motif_chain_blast(pdbs, map):
        ## Run
        for pdb in pdbs:
            mot_chain_fasta, map = get_FASTA_from_pdb(pdb, pdbs[pdb]['motif_chain'], map)
            acci = pdbs[pdb]['acci']
            if os.path.isfile('fasta/'+acci+'.fasta') == False:
                acci_fasta = ''
                with urllib.request.urlopen('https://www.uniprot.org/uniprot/'+acci+'.fasta') as response:
                    acci_fasta = str(response.read().decode('UTF-8'))
                #print (acci_fasta)
                open('fasta/'+acci+'.fasta', 'w').write(acci_fasta)
            #print (html)
            #sys.exit()
            open(outdir+acci+'_motif_chain.fasta', 'w').write(mot_chain_fasta)
            os.system('makeblastdb -in '+outdir+acci+'_motif_chain.fasta -dbtype "prot" -out '+outdir+acci+'_motif_chain_pdb -title '+outdir+acci+'_motif_chain_pdb')
            os.system('blastp -query fasta/'+acci+'.fasta -db '+outdir+acci+'_motif_chain_pdb -out '+outdir+acci+'_motif_chain_blastp.txt')
            #sys.exit()
            ## Read
            dic = {}
            for line in open(outdir+acci+'_motif_chain_blastp.txt', 'r'):
                #print (line.split())
                if len(line.split()) > 0:
                    if line[0] == '>':
                        pdb_chain = line.split('>')[1].replace('\n', '')
                    elif line.split()[0] == 'Query':
                        startq = int(line.split()[1])
                        endq = int(line.split()[3].replace('\n', ''))
                        query = line.split()[2]
                    elif line.split()[0] == 'Sbjct':
                        starts = int(line.split()[1])
                        ends = int(line.split()[3].replace('\n', ''))
                        sbjct = line.split()[2]
                        for q, s in zip(query, sbjct):
                            #print (q, s)
                            if q not in ['.', '-'] and s not in ['.', '-']:
                                dic[startq] = starts
                                startq += 1
                                starts += 1
                            elif q not in ['.', '-']:
                                startq += 1
                            else:
                                starts += 1
            #print (dic)
            pdbs[pdb]['acci_fasta_to_pdb_AA_fasta'] = dic
            mot_chain = pdbs[pdb]['motif_chain']
            #print (map[pdb][mot_chain])
            x = []
            for position in pdbs[pdb]['acci_fasta_to_pdb_AA_fasta']:
                if int(position) >= pdbs[pdb]['starti'] and position <= pdbs[pdb]['endi']:
                    x.append(map[pdb][mot_chain][dic[position]])
            pdbs[pdb]['motif_chains'] = []
            if len(x) > 0:
                pdbs[pdb]['motif_chains'].append(mot_chain+':'+str(x[0])+'-'+str(x[-1]))
            #sys.exit()

    def run_and_read_elm_domain_chain_hmmscan(pdbs, map):
        ## Run
        #print(pdbs)
        #sys.exit()
        dom_chain_fasta = ''
        for pdb in pdbs:
            for dom_chain in pdbs[pdb]['domain_chain']:
                fasta, map = get_FASTA_from_pdb(pdb, dom_chain, map)
                dom_chain_fasta += fasta
        open(outdir+input_protein+'_domain_chain.fasta', 'w').write(dom_chain_fasta)
        os.system('hmmscan --noali --domtblout '+outdir+input_protein+'_domain_chain_hmmscan.txt /home/gurdeep/projects/DB/PFAM/Pfam-A.hmm '+outdir+input_protein+'_domain_chain.fasta')

        ## Read
        #pdbs[pdb]['domain_chains'] = []
        for line in open(outdir+input_protein+'_domain_chain_hmmscan.txt', 'r'):
            if line[0] != '#':
                #print (line)
                pfam = line.split()[1].split('.')[0]
                if pfam == domain:
                    pdb = line.split()[3].split('_')[0]
                    if 'domain_chains' not in pdbs[pdb]:
                        pdbs[pdb]['domain_chains'] = ''
                    chain = line.split()[3].split('_')[1]
                    chain_start = int(line.split()[17])
                    chain_end = int(line.split()[18])
                    #print (map[pdb][chain])
                    #pdbs[pdb]['domain_chains'].append(chain+':'+str(map[pdb][chain][chain_start])+'-'+str(map[pdb][chain][chain_end]))
                    pdbs[pdb]['domain_chains'] += chain+':'+str(map[pdb][chain][chain_start])+'-'+str(map[pdb][chain][chain_end]) + ';'
                    #print (pdbs[pdb]['domain_chains'])
                    #sys.exit()

    def read_SIFTS_chain_to_domain():
        domain_pdb_and_chain = {}
        for line in gzip.open('/home/gurdeep/projects/DB/SIFTS/pdb_chain_pfam.tsv.gz', 'rt'):
            if line[0] != '#' and line.split()[0] != 'PDB':
                pfamid = line.split()[3]
                if pfamid == domain:
                    pdbid = line.split()[0].upper()
                    chainid = line.split()[1]
                    if pdbid not in domain_pdb_and_chain:
                        domain_pdb_and_chain[pdbid] = []
                    domain_pdb_and_chain[pdbid].append(chainid)
        return domain_pdb_and_chain

    def read_SIFTS_acc_to_chain(domain_pdb_and_chain):
        pdb_acc_chain = {}
        for line in gzip.open('/home/gurdeep/projects/DB/SIFTS/pdb_chain_uniprot.tsv.gz', 'rt'):
            if line[0] != '#' and line.split()[0] != 'PDB':
                pdbid = line.split()[0].upper()
                chainid = line.split()[1]
                acc = line.split()[2]
                if pdbid in domain_pdb_and_chain:
                    if pdbid not in pdb_acc_chain:
                        pdb_acc_chain[pdbid] = {}
                    pdb_acc_chain[pdbid][acc] = chainid
        return pdb_acc_chain

    def download_pdb_cif(pdb):
        if os.path.isfile('pdbs/'+pdb+'.cif') == False:
            os.system('wget https://files.rcsb.org/view/'+pdb+'.pdb -O '+'pdbs/'+pdb+'.pdb')
            os.system('wget https://files.rcsb.org/view/'+pdb+'.cif -O '+'pdbs/'+pdb+'.cif')

    pdbs={}
    if source == '3DID':
        flag  = 0
        for line in gzip.open('/home/gurdeep/projects/DB/3did/3did_dmi_flat.gz', 'rt'):
            print (line)
            if line.split()[0] == '#=ID':
                pdbs = {}
                #print (domain, line)
                if line.split()[1] == domain and line.split()[3] == motif:
                    flag = 1
            elif line.split()[0] == '#=3D' and flag == 1:
                pdb = line.split()[1].upper()
                if os.path.isfile(outdir+pdb+'_'+input_protein+'_'+motif+'_'+domain+'_'+str(startm)+'_'+str(endm)+'_'+source+'.pse') == True:
                    print('Found'+pdb+'_'+input_protein+'_'+motif+'_'+domain+'_'+str(startm)+'_'+str(endm)+'_'+source+'.pse')
                    print (line)
                if pdb not in pdbs:
                    pdbs[pdb] = {}
                    pdbs[pdb]['motif_chains'] = []
                    pdbs[pdb]['domain_chains'] = ''
                pdbs[pdb]['motif_chains'].append(line.split()[3])
                pdbs[pdb]['domain_chains'] += line.split()[2]+';'
                download_pdb_cif(pdb)
            elif line[:2] == '//' and flag == 1:
                #print (pdbs)
                #run_blast(input_protein, pdbs, domain, motif)
                flag = 0
                break
    elif source == 'ELM':
        domain_pdb_and_chain = read_SIFTS_chain_to_domain()
        pdb_acc_chain = read_SIFTS_acc_to_chain(domain_pdb_and_chain)
        elm_pdbs = {}
        for line in open('/home/gurdeep/projects/DB/ELM/elm_instances.tsv', 'r'):
            if line[0] != '#':
                if line.split('\t')[2].replace('"', '') == motif:
                    #print (line)
                    #elm_pdbs += line.split('\t')[11].replace('"', '').split()
                    for pdb in line.split('\t')[11].replace('"', '').split():
                        #print (pdb)
                        if pdb in domain_pdb_and_chain:
                            #print (domain_pdb_and_chain[pdb])
                            #sys.exit()
                            starti = int(line.split('\t')[6].replace('"', ''))
                            endi = int(line.split('\t')[7].replace('"', ''))
                            acci = line.split('\t')[4].replace('"', '')
                            print (acci)
                            if acci in pdb_acc_chain[pdb]:
                                motif_chain = pdb_acc_chain[pdb][acci]
                                #elm_pdbs[pdb]['domain_chain'] = ''
                                for chain in domain_pdb_and_chain[pdb]:
                                    if chain != motif_chain:
                                        if pdb not in elm_pdbs:
                                            elm_pdbs[pdb] = {}
                                            elm_pdbs[pdb]['domain_chain'] = []
                                            elm_pdbs[pdb]['motif_chain'] = pdb_acc_chain[pdb][acci]
                                            elm_pdbs[pdb]['starti'] = starti
                                            elm_pdbs[pdb]['endi'] = endi
                                            elm_pdbs[pdb]['acci'] = acci
                                            download_pdb_cif(pdb)
                                        elm_pdbs[pdb]['domain_chain'].append(chain)
                                        #break
        print (elm_pdbs)
        #sys.exit()
        map = {}
        if len(elm_pdbs) > 0:
            run_and_read_elm_motif_chain_blast(elm_pdbs, map)
            run_and_read_elm_domain_chain_hmmscan(elm_pdbs, map)
            pdbs = elm_pdbs
        else:
            return [motif]
        #pdbs = filter_pdbs(elm_pdbs, domain)
        #sys.exit()

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
    return []

'''
input_protein = 'P0DTD1_nsp9'
motif = 'Dynein_light_LIG_0-14'
domain = 'Dynein_light'
startm = 1127
endm = 129
source = '3DID'

main(input_protein, motif, domain, startm, endm, source)
sys.exit()
'''

elm_motifs_wo_pdb = []
num = 0
for line in gzip.open('../covid19/data/interactions/phosphosite/analysis_elm_motif_domain_kinase.txt.gz', 'rt'):
    if line[0] != '#':
        #print (line.split('\t'))
        input_protein = line.split('\t')[0]
        motif = line.split('\t')[2]
        startm = line.split('\t')[6]
        endm = line.split('\t')[7]
        source = line.split('\t')[-1].replace('\n', '')
        if source == '3DID':
            domain = line.split('\t')[9]
        else:
            domain = line.split('\t')[8]
        #if source == '3DID':
        if motif not in elm_motifs_wo_pdb:
            elm_motifs_wo_pdb += main(input_protein, motif, domain, startm, endm, source)
                #break
print ('done')
