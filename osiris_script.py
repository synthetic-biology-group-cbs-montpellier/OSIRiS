#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 15:11:55 2019

@author: sarahguiziou

"""

__author__  = 'Sarah Guiziou <guiziou@cbs.cnrs.fr>'
__license__ = 'MIT'
__version__ = '1.0'

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation
import string
import pandas as pd

# To reverse complement sequences
old_chars = "ACGTacgt"
replace_chars = "TGCAtgca"
tab = string.maketrans(old_chars,replace_chars)
   
BF_site=[]
PF_site=[]
LF_site=[]
RF_site=[]
BR_site=[]
PR_site=[]
LR_site=[]
RR_site=[]
    
" IMPORT A GENBANK FILE AND RETURN THE DNA SEQUENCE AND LIST OF FEATURES "

def import_gb(directory, name):
    
    from Bio import SeqIO
    gb_file = directory+'/'+name+'.gb'
    for gb_record in SeqIO.parse(open(gb_file,"r"), "genbank") :
        # now do something with the record
        #print "Name %s, %i features" % (gb_record.name, len(gb_record.features))
        #print repr(gb_record.seq)
        DNA_seq=str(gb_record.seq)
        feat=gb_record.features
        return DNA_seq, feat
    
    
#[DNA_seq, feat]=import_gb('/Users/sarahguiziou/Desktop/OSIRIS','3-inputV3')


" GENERATE THE INTEGRASE SITE LIST FROM THE EXCEL FILE AND THE LIST OF INTEGRASE OF INTEREST"

def integrase_sites(int_list, directory_seq):
    
    list_seq=pd.read_csv(directory_seq, sep=";",usecols=[0,1], names=['name','seq'])
    dico_seq={}
    
    for X in range(0,len(list_seq)):
        dico_seq[list_seq['name'][X].lower()]=list_seq['seq'][X]
    
    global BF_site
    global PF_site
    global LF_site
    global RF_site
    global BR_site
    global PR_site
    global LR_site
    global RR_site
    
    Error=False
    for integrase in int_list:
        if dico_seq.has_key('attb_'+integrase.lower()):
            BF_site.append(dico_seq.get('attb_'+integrase.lower()))
        else:
            Error=True
        if Error!=True and dico_seq.has_key('attp_'+integrase.lower()):
            PF_site.append(dico_seq.get('attp_'+integrase.lower()))
        else:
            Error=True
        if Error!=True and dico_seq.has_key('attl_'+integrase.lower()):
            LF_site.append(dico_seq.get('attl_'+integrase.lower()))
        else:
            Error=True
        if Error!=True and dico_seq.has_key('attr_'+integrase.lower()):
            RF_site.append(dico_seq.get('attr_'+integrase.lower()))
        else:
            Error=True
    
    if Error!=True:       
        for a in range(len(BF_site)):
            BR_site.append(BF_site[a].translate(tab)[::-1])
        for a in range(len(PF_site)):
            PR_site.append(PF_site[a].translate(tab)[::-1])
        for a in range(len(LF_site)):
            LR_site.append(LF_site[a].translate(tab)[::-1])
        for a in range(len(RF_site)):
            RR_site.append(RF_site[a].translate(tab)[::-1])
            
    return Error

#print integrase_sites(['Tp901', 'bxB1'], '/Users/sarahguiziou/Desktop/osiris/List_integrases.csv')

            
" FROM A LIST OF INTEGRASE, GENERATE THE PERMUTATION OF INTEGRASES"

def generate_list_int_occurence(number_int):
    import itertools
    list_int=range(0, number_int)
    return list(itertools.permutations(list_int))
    
#print generate_list_int_occurence(4)

" RECOMBINE A DNA SEQUENCE BY ONE INTEGRASE"     
   
def recombination_DNA_int(DNA_seq, int_num, feat):
    
    end='n'
    
    if DNA_seq.find(BF_site[int_num])!=-1:
        position_site1=DNA_seq.find(BF_site[int_num])
        if DNA_seq.find(PF_site[int_num])!=-1:
            position_site2=DNA_seq.find(PF_site[int_num])
            type_rec='excision'
            new_site='L'
        elif DNA_seq.find(PR_site[int_num])!=-1: 
            position_site2=DNA_seq.find(PR_site[int_num])
            type_rec='inversion'
            new_site='L'
            
        else:
            end='y'
            
    elif DNA_seq.find(BR_site[int_num])!=-1:
        position_site1=DNA_seq.find(BR_site[int_num])
        if DNA_seq.find(PF_site[int_num])!=-1:
            position_site2=DNA_seq.find(PF_site[int_num])
            type_rec='inversion'
            new_site='R'
        elif DNA_seq.find(PR_site[int_num])!=-1: 
            position_site2=DNA_seq.find(PR_site[int_num])
            type_rec='excision'
            new_site='R'
            
        else:
            end='y'
            
    else:
        end='y'
    
    if end=='n':
        if type_rec=='inversion':
            if position_site1<position_site2:
                #DNA_seq[(position_site1+len(BF_site[int_num])):position_site2].translate(tab)[::-1]
                DNA_seq =string.replace(DNA_seq, DNA_seq[(position_site1+len(BF_site[int_num])):position_site2], DNA_seq[(position_site1+len(BF_site[int_num])):position_site2].translate(tab)[::-1])
                #DNA_seq[(position_site1+len(BF_site[int_num])):position_site2].reverse_complement()
                if new_site=='L':
                    DNA_seq =string.replace(DNA_seq, BF_site[int_num], LF_site[int_num])
                    DNA_seq=string.replace(DNA_seq, PR_site[int_num], RR_site[int_num])
                else:
                    DNA_seq =string.replace(DNA_seq, BR_site[int_num], RR_site[int_num])
                    DNA_seq=string.replace(DNA_seq, PF_site[int_num], LF_site[int_num])                
            else:
                DNA_seq =string.replace(DNA_seq, DNA_seq[(position_site2+len(PF_site(int_num))):position_site1], DNA_seq[(position_site2+len(PF_site(int_num))):position_site1].translate(tab)[::-1])
                #DNA_seq[(position_site2+len(PF_site(int_num))):position_site1].translate(tab)[::-1]
                #DNA_seq[(position_site2+len(PF_site(int_num))):position_site1].reverse_complement()
                if new_site=='L':
                    DNA_seq =string.replace(DNA_seq, BF_site[int_num], RF_site[int_num])
                    DNA_seq=string.replace(DNA_seq, PR_site[int_num], LR_site[int_num]) 
                else:
                    DNA_seq =string.replace(DNA_seq, BR_site[int_num], LR_site[int_num])
                    DNA_seq=string.replace(DNA_seq, PF_site[int_num], RF_site[int_num])
        
        if type_rec=='excision':
            if position_site1<position_site2:
                if new_site=='L':
                    DNA_seq=DNA_seq[:position_site1]+LF_site[int_num]+DNA_seq[position_site2:]
                if new_site=='R':
                    DNA_seq=DNA_seq[:position_site1]+RR_site[int_num]+DNA_seq[position_site2:]
            else:
                if new_site=='L':
                    DNA_seq=DNA_seq[:position_site2]+RF_site[int_num]+DNA_seq[position_site1:]
                if new_site=='R':
                    DNA_seq=DNA_seq[:position_site2]+LR_site[int_num]+DNA_seq[position_site1:]  
     
    return end, DNA_seq, feat

" GENERATE INTEGRASE SITE FEATURE"

def integrase_site_feature(DNA_seq, name_int):
    
    list_site_feat=[]
    # location beginning, location end, Name, strand
    for integrase in range(0,len(name_int)):
        if DNA_seq.find(BF_site[integrase])!=-1:
            list_site_feat.append([DNA_seq.find(BF_site[integrase]), len(BF_site[integrase])+DNA_seq.find(BF_site[integrase]), 'attB_'+name_int[integrase], 1])
        elif DNA_seq.find(BR_site[integrase])!=-1:
            list_site_feat.append([DNA_seq.find(BR_site[integrase]), len(BR_site[integrase])+DNA_seq.find(BR_site[integrase]), 'attB_'+name_int[integrase], -1])
        if DNA_seq.find(PF_site[integrase])!=-1:
            list_site_feat.append([DNA_seq.find(PF_site[integrase]), len(PF_site[integrase])+DNA_seq.find(PF_site[integrase]), 'attP_'+name_int[integrase], 1])
        elif DNA_seq.find(PR_site[integrase])!=-1:
            list_site_feat.append([DNA_seq.find(PR_site[integrase]), len(PR_site[integrase])+DNA_seq.find(PR_site[integrase]), 'attP_'+name_int[integrase], -1])
        if DNA_seq.find(LF_site[integrase])!=-1:
            list_site_feat.append([DNA_seq.find(LF_site[integrase]), len(LF_site[integrase])+DNA_seq.find(LF_site[integrase]), 'attL_'+name_int[integrase], 1])
        elif DNA_seq.find(LR_site[integrase])!=-1:
            list_site_feat.append([DNA_seq.find(LR_site[integrase]), len(LR_site[integrase])+DNA_seq.find(LR_site[integrase]), 'attL_'+name_int[integrase], -1])
        if DNA_seq.find(RF_site[integrase])!=-1:
            list_site_feat.append([DNA_seq.find(RF_site[integrase]), len(RF_site[integrase])+DNA_seq.find(RF_site[integrase]), 'attR_'+name_int[integrase], 1])
        elif DNA_seq.find(RR_site[integrase])!=-1:
            list_site_feat.append([DNA_seq.find(RR_site[integrase]), len(RR_site[integrase])+DNA_seq.find(RR_site[integrase]), 'attR_'+name_int[integrase], -1])
            
    return list_site_feat  
   
" CREATE A GB FILE FROM A DNA SEQ AND A LIST OF FEATURES  "         

def gb_create(DNA_seq, list_feat, name, directory):

    # creation of the formated DNA sequence for the genbank file
    seq_final=Seq(DNA_seq, IUPAC.unambiguous_dna)
    record = SeqRecord(seq_final,
                   id='NA', # random accession number
                   name='',
                   description='Synthetic sequence')

    # loop in sequence list to generate the genbank feature of the sequence
    for feat in list_feat:
        
        feature = SeqFeature(FeatureLocation(start=feat[0], end=feat[1]), id=feat[2], type=feat[2], strand=feat[3])
        record.features.append(feature)

    # Save as GenBank file
    output_file = open(directory+'/'+name+'.gb', 'w')
    SeqIO.write(record, output_file, 'genbank')
    output_file.close()

" GENERATE THE DIFFERENT DNA SEQUENCES FROM AN INTEGRASE COMBINATION " 
  
def search_for_integrase(DNA_seq, feat, combination, directory, input_file, name_int):
    
    list_integrase_used=[]
    list_end=[]
    
    for x in combination:
        [end, DNA_seq, feat]=recombination_DNA_int(DNA_seq, x, feat)
        list_integrase_used.append(name_int[x])
        list_end.append(end)
        
        name=input_file
        for y in list_integrase_used:
            name+='_'+y
            
        feat=integrase_site_feature(DNA_seq, name_int)
        gb_create(DNA_seq, feat, name, directory)

        f=open(directory+'/'+name+'txt', 'w') 
        f.write('>'+name+'\n')
        f.write(str(DNA_seq))
        f.write('\n\n')
        f.close()
        
    return list_end

" FINAL SCRIPT "
def OSIRIS_script(directory, input_file, name_int, directory_seq):
    
    [DNA_seq, feat]=import_gb(directory,input_file)
    
    Error=integrase_sites(name_int, directory_seq)
    
    if Error:
        print 'The sites of the integrase list have not been found in the corresponding csv file.'
    else:
        g=open(directory+'/results_of_recombination.txt', 'w')
        g.write('For the sequence named '+input_file+': \n\n')
        list_order_int=generate_list_int_occurence(len(name_int))
    
        for order in list_order_int:
            list_feat=search_for_integrase(DNA_seq, feat, order, directory, input_file, name_int)
            for x in range(0,len(name_int)):
                for y in range(0,x+1):
                    id_int=order[y]
                    g.write(name_int[id_int])
                    if y!=x and x!=0:
                        g.write(' then ')
                if list_feat[x]=='n':
                    g.write(': recombination.\n')
                else:
                    g.write(': no recombination.\n')
        g.close()


""" TO RUN THE OSIRiS SCRIPT """
"""Input of the script: 
    - directory: the directory where you want the DNA sequence to be generated and where the input DNA sequence is.
    - input_file: name of the input DNA sequence of interest to be recombined.
    - name_int: list of name of integrases of interest.
    - directory_list_sites: path of the csv file containing the list of integrase sites.

Output of the script:
    - Generation of all the DNA sequences in FASTA format corresponding to the different intermediate recombination sites 
    of the different orders of occrences of integrases.
    - Generation of a text file resuming the occurence or not of recombination in specific integrase occurence states.
"""

"EXAMPLE"
#name_int=['Bxb1', 'Tp901', 'Int5']
#directory='/Users/sarahguiziou/Desktop/osiris'
#input_file='3-inputV3'
#directory_list_sites=directory+'/List_Integrases.csv'

#OSIRIS_script(directory, input_file, name_int, directory_list_sites)
