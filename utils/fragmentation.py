# -*- coding: utf-8 -*-
"""
Created on Thu May 20 13:09:29 2021

@author: user
"""

import rdkit
from rdkit import Chem
from itertools import combinations
from rdkit.Chem import GetShortestPath
import argparse
import re


def parse_args():

    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-input_file", "-i", required=True,
                        help="file of the macrocycles represented as SMILES")

    parser.add_argument("-output_file", "-o", required=True,
                        help="Generated acyclic-macrocyclic SMILES pairs")

    return parser.parse_args()


def singlebonds_in_macring(mol):
    
    singlebonds=[]
    
    rings=list(mol.GetRingInfo().AtomRings())
    rings.sort(key=lambda x:len(x))
    maxring=rings[-1]
    
    for b in mol.GetBonds():
        if b.GetBondType() != Chem.rdchem.BondType.SINGLE:
            continue
        else:
            id1=b.GetBeginAtomIdx()
            id2=b.GetEndAtomIdx()
            if id1 in maxring and id2 in maxring:
                singlebonds.append(b.GetIdx())
                
    return(singlebonds)
 
    
def filter_chain_length(chain):
    
    allatoms=chain.GetNumHeavyAtoms()
    
    dummy_ids=[a.GetIdx() for a in chain.GetAtoms() if a.GetAtomicNum()==0]
    length=len(GetShortestPath(chain,dummy_ids[0],dummy_ids[1]))-2
    
    action= (allatoms>2) and (2<length<10) and (length/allatoms>=0.6)
    
    return action,length


def filter_chain_ring(chain):
    
    rings=chain.GetRingInfo().AtomRings()
    
    if len(rings)==0:
        action=True
        ring_numatom=0
    
    elif len(rings)==1:
        ring_numatom=max([len(k) for k in rings])
        if ring_numatom<7:
            action=True
        
    else:
        action=False
        ring_numatom=8
    
    return action,ring_numatom
    
    
def fragmentation(mol):
    
    frags=[]
    bonds=singlebonds_in_macring(mol)
    
    if len(bonds)==0: return frags
    
    for b in combinations(bonds,2):
        
        frags_mol=Chem.FragmentOnBonds(mol,list(b))
        frags_mols=Chem.GetMolFrags(frags_mol,asMols=True)
        
        if len(frags_mols) !=2: continue
        
        frag1_mol=frags_mols[0]
        frag2_mol=frags_mols[1]
        frag1_numatoms=frag1_mol.GetNumHeavyAtoms()
        frag2_numatoms=frag2_mol.GetNumHeavyAtoms()
        
        if frag1_numatoms<frag2_numatoms:
            chain=frag1_mol
            frag=frag2_mol
        else:
            chain=frag2_mol
            frag=frag1_mol
            
        if (chain.GetNumHeavyAtoms()/mol.GetNumHeavyAtoms())<0.25:
            
            if filter_chain_length(chain)[0] and filter_chain_ring(chain)[0]:
    
                chain_smiles=remove_dummy_index(chain)
                frag_smiles=remove_dummy_index(frag)
                s='N_'+str(filter_chain_length(chain)[1])+'.'+str(chain_smiles)+'.'+str(frag_smiles)
                frags.append(s)
                
    return list(set(frags))
    

def remove_dummy_index(frag,pattern=re.compile(r'\[\d+.\]')):
    
    smiles=Chem.MolToSmiles(frag)
    
    dummy_index = pattern.findall(smiles)
    
    for i in dummy_index:
        smiles = smiles.replace(i, '[*]')
    
    cano_smiles=Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
   
    return cano_smiles    


def macs_fragmentation(f_in,f_out_mmps): 
    
    f=open(f_in)
    line=f.readline()
    
    fw=open(f_out_mmps,'w')
#    fw2=open(f_out_mmps_num,'w')
    
    while line:
        
#        chem_id=line.strip().split()[0]
#        smiles=line.strip().split()[1]
        smiles=line.strip().split()[1]
        
        if len(smiles)<201:   
            mol=Chem.MolFromSmiles(smiles)
            
            if mol is not None:
                frags=fragmentation(mol)
#                fw2.write(line.strip()+' '+str(len(mmps))+'\n')
            
                for frag in frags:
                    s=frag+'>'+smiles
                    fw.write(s+'\n')
    
        line=f.readline()
                
    f.close()
    fw.close()
#    fw2.close()
    
def main():
    opt = parse_args()
    macs_fragmentation(opt.input_data_path,opt.output_data_path)

if __name__ == "__main__":
    main()      