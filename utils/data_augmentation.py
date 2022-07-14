# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 21:20:12 2022

@author: yydiao
"""

import argparse
import rdkit
from rdkit import Chem
import random

def parse_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-input_src_file', '-is', required=True,                        
                        help='file of canonical source SMILES')    
    parser.add_argument('-input_tgt_file', '-it', required=True,                        
                        help='file of canonical target SMILES')
    parser.add_argument('-output_path', '-o', required=True,
                        help='path of the augmented src/tgt file')
    parser.add_argument('-augmentation_number', '-n', required=True,
                        help='the number for data augmentation')
    
#   if True, the output file will contain one copy of canonical SMILES and n-1 randomized SMILES
#   else, the output file will only contain n-1 randomized SMILES
    parser.add_argument('-copy_canonical_smiles', '-c', action='store_true', required=True,
                        help='copy the canonical SMILES to the output file')

    return parser.parse_args()


def randomize_src(smiles, random_type="restricted"):   
    
    frag=Chem.MolFromSmiles(smiles)
    
#    if frag is None:
#        return smiles
    
    if random_type == "unrestricted":
        return Chem.MolToSmiles(frag, canonical=False, doRandom=True, isomericSmiles=False)
    
    if random_type == "restricted":
        new_atom_order = list(range(frag.GetNumAtoms()))
        random.shuffle(new_atom_order)
        random_smiles = Chem.RenumberAtoms(frag, newOrder=new_atom_order)
        return Chem.MolToSmiles(random_smiles, canonical=False, isomericSmiles=False)
    
    
def substructure_aligned_tsmiles(randomized_src,mol_smiles):

#    Given an acyclic source randomized SMILES, generate a randomized macrocyclic target SMILES, 
#    the atom numbers of which are reordered according to the acyclic source randomized SMILES,
#    this operation makes the input and output SMILES align better.
 
    randomized_src=randomized_src.replace('n1*','[nH]1').replace('n2*','[nH]2').\
                replace('n3*','[nH]3').replace('n(*)','[nH]').\
                replace('*n','[nH]').replace('(*)','').replace('*','')
    frag=Chem.MolFromSmiles(randomized_src)
    
    if frag is None:
        return(mol_smiles)
    
    mol=Chem.MolFromSmiles(mol_smiles)
    mol_atom_list=list(range(mol.GetNumAtoms()))    
    matches=mol.GetSubstructMatches(frag)
    
    if len(matches)==0:
        return(mol_smiles)
        
    elif len(matches)==1:
        match=list(matches[0])
        link_atom_list=[a for a in mol_atom_list if a not in match]
        new_order=match+link_atom_list
        random_mol = Chem.RenumberAtoms(mol, newOrder=new_order)
        new_mol_smiles=Chem.MolToSmiles(random_mol, canonical=False, isomericSmiles=False)
        return(new_mol_smiles)
        
    elif len(matches)>1:
        for match in matches:
            match=list(match)
            match_sort = sorted(match, reverse=True)
            mol_rw = Chem.RWMol(mol)
            for idx in match_sort:
                mol_rw.RemoveAtom(idx)
                linkers = Chem.Mol(mol_rw)
                if len(Chem.rdmolops.GetMolFrags(linkers)) == 1:
                    link_atom_list=[a for a in mol_atom_list if a not in match]  
                    new_order=match+link_atom_list
        random_mol = Chem.RenumberAtoms(mol, newOrder=new_order)
        new_mol_smiles=Chem.MolToSmiles(random_mol, canonical=False, isomericSmiles=False)
        return(new_mol_smiles)
        
def write_file(input_src_file,input_tgt_file,dir,n,copy_canonical_smiles=True):
#    opt = parse_args()
    
    output_src_file = dir + 'src' + '-aug' + str(n)
    output_tgt_file = dir + 'tgt' + '-aug' + str(n)
    
    src_lines=open(input_src_file).readlines()    
    cano_srcs=[''.join(line.strip().split()[1:]) for line in src_lines]
    lens_srcs=[line.strip().split()[0] for line in src_lines]
    
    tgt_lines=open(input_tgt_file).readlines()  
    cano_tgts=[''.join(line.strip().split()) for line in tgt_lines]
    
    nums=len(cano_srcs)
    
    w1=open(output_src_file,'w')
    w2=open(output_tgt_file,'w')
    
    for i in range(nums):
        cano_smiles_src=cano_srcs[i]
        cano_smiles_tgt=cano_tgts[i]
        
        if copy_canonical_smiles:
            w1.write(src_lines[i])
            w2.write(tgt_lines[i])
            
        for j in range(int(n)-1):            
            randomized_smiles_src=randomize_src(cano_smiles_src, random_type="restricted")
            w1.write(lens_srcs[i]+' '+' '.join(randomized_smiles_src)+'\n')
            randomized_smiles_tgt=substructure_aligned_tsmiles(randomized_smiles_src,cano_smiles_tgt)
            w2.write(' '.join(randomized_smiles_tgt)+'\n')
    w1.close()
    w2.close()
    
    
def main():
    opt = parse_args()
    input_src_file=opt.input_src_file
    input_tgt_file=opt.input_tgt_file
    dir=opt.output_path
    n=opt.augmentation_number
    copy_canonical_smiles=opt.copy_canonical_smiles
    write_file(input_src_file,input_tgt_file,dir,n,copy_canonical_smiles)
        
    
if __name__ == '__main__':
    main()   