# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 10:49:43 2022

@author: yydiao
"""

import argparse
from rdkit import Chem
#import pandas as pd
#from rdkit.Chem import AllChem
#from rdkit import DataStructs
import numpy
from numpy import *


def parse_args():
    parser = argparse.ArgumentParser(description='')
    
    parser.add_argument('-input_tgt_train_file', '-itrain', required=True,                        
                        help='file of target SMILES of the train dataset without augmentation')   
    
    parser.add_argument('-input_tgt_test_file', '-itest', required=True,                        
                        help='file of target SMILES of the test dataset without augmentation')
    
    parser.add_argument('-input_pred_file', '-ipred', required=True,                        
                        help='file of predicted SMILES with or without augmentation')
    
    parser.add_argument('-output_file', '-o', required=True,
                        help='output file with evaluation metrics reported')
    
    parser.add_argument('-augmentation_number', '-n', required=True,
                        help='the number for data augmentation')
    
    parser.add_argument('-best_number', '-bn', required=True,
                        help='the output number for each input SMILES sequence')
    
    return parser.parse_args()

def canonicalize_smiles(s):
    mol=Chem.MolFromSmiles(s)
    if mol is not None:
        return Chem.MolToSmiles(mol, isomericSmiles=True)
    else:
        return None
 
       
def get_train_set(file_train,cano=True):
    f=open(file_train,'r')
    train_smiles=[''.join(line.strip().split()) for line in f.readlines()]
    if cano:
        train_smiles=[canonicalize_smiles(s) for s in train_smiles]
    return set(train_smiles)


def get_test_smiles(file_test,cano=True):
    f=open(file_test,'r')
    test_smiles=[''.join(line.strip().split()) for line in f.readlines()]
    if cano:
        test_smiles=[canonicalize_smiles(s) for s in test_smiles]
    return test_smiles

def get_pred_smiles(file_pred,cano=True):
    f=open(file_pred,'r')
    pred_smiles=[''.join(line.strip().split()) for line in f.readlines()]
    if cano:
        pred_smiles=[canonicalize_smiles(s) for s in pred_smiles]
    return pred_smiles

def get_pred_dict(smiles,best_num):
    pred={i:[] for i in range(1,int(best_num)+1)}
    for i,smile in enumerate(smiles):
        pred[(i%int(best_num))+1].append(smile)        
    return pred


def get_rank(pred_dict,test_smiles,best_num):
    rank=[0]*len(test_smiles)
    for i in range(len(test_smiles)):
        target=test_smiles[i]
        for j in range(1,int(best_num)+1):
            if target==pred_dict[j][i]:
                rank[i]=j
                break
    return rank


def get_rank_sim(pred_dict,test_smiles,best_num,cutoff=1.0):
    rank=[0]*len(test_smiles)
    for i in range(len(test_smiles)):
        target=test_smiles[i]
#        mol_target=Chem.MolFromSmiles(target)
#        fp1 = AllChem.GetMorganFingerprint(mol_target,2)
        for j in range(1,int(best_num)+1):
            if pred_dict[j][i] == target:
                rank[i]=j

    return rank
        
def mac_num(smiles):
    mac=[]
    for s in smiles:
        mol=Chem.MolFromSmiles(s)
        if mol is not None:
            rings=mol.GetRingInfo().AtomRings()
            ring_numatom=[len(k) for k in rings]
            mac_num=[k for k in ring_numatom if k>11]
            if len(mac_num)>0:
                mac.append(s)
    return mac


file_train='F:/subjects/maclink_AI/chembl29/tgt_reorder_from_frags2/tgt-train-ysrc'

file_test='F:/subjects/maclink_AI/chembl29/tgt_reorder_from_frags2/tgt-external-zinc3'
file_pred='F:/subjects/maclink_AI/chembl29/tgt_reorder_from_frags2/predictions_tgt_reorder_from_frags_cano_model_average.pt_on_zinc3a10_best10.txt'


def calculate_metrics(file_train,file_test,file_pred,output_file,best_num,n):
    
#    num_pre_smiles=best_num*n
    train_set=get_train_set(file_train,cano=True)
    test_smiles=get_test_smiles(file_test,cano=True)
    pred_smiles=get_pred_smiles(file_pred,cano=True)
    
    fw=open(output_file,'w')
    total=len(test_smiles)

    pred_smiles2=numpy.array(pred_smiles).reshape(total,int(n),int(best_num))

    pred_smiles3=[[] for i in range(int(n))]
    for i in range(int(n)):
        for j in range(total):
            pred_smiles3[i].extend(list(pred_smiles2[j][i]))
        

    dict_recovery={i:[] for i in range(1,int(best_num)+1)}
    dict_valid={i:[] for i in range(1,int(best_num)+1)}
    dict_unique={i:[] for i in range(1,int(best_num)+1)}
    dict_novelty={i:[] for i in range(1,int(best_num)+1)}
    dict_macror={i:[] for i in range(1,int(best_num)+1)}


    for k in range(1,int(n)+1):
        smiles=pred_smiles3[k-1]
    
        pred_dict=get_pred_dict(smiles,best_num=int(best_num))
        rank=get_rank(pred_dict,test_smiles,best_num=int(best_num))
#    rank_smi=get_rank_sim(pred_dict,test_smiles,best_num,cutoff=1.0) 
    

        correct=0
#    simmol=0
        pred=[]
        mac=[]
        for i in range(1,int(best_num)+1):
            correct+=rank.count(i)
#        simmol+=rank_smi.count(i)
            pred.extend(pred_dict[i])
            valid=len(pred)-pred.count(None)
            unique=len(set(pred)-{None})
#            pred_unique=list(set(pred)-{None})
            mac.extend(mac_num(list(set(pred_dict[i])-{None})))
            macros=len(set(mac))
            novelty=len(set(pred)-{None}-train_set)

            dict_recovery[i].append(correct/total)
            dict_valid[i].append(valid/(total*i))
            dict_unique[i].append(unique/valid)
            dict_novelty[i].append(novelty/unique)
            dict_macror[i].append(macros/unique)

    for i in range(1,int(best_num)+1):
        fw.write('top_{} : recovery {:.2f}±{:.2f}% || valid {:.2f}±{:.2f}% || uniqueness {:.2f}±{:.2f}% || novelty {:.2f}±{:.2f}% || macror {:.2f}±{:.2f}%'
                .format(i, mean(dict_recovery[i])*100, std(dict_recovery[i])*100, 
                   mean(dict_valid[i])*100, std(dict_valid[i])*100, 
                   mean(dict_unique[i])*100, std(dict_unique[i])*100, 
                   mean(dict_novelty[i])*100, std(dict_novelty[i])*100, 
                   mean(dict_macror[i])*100, std(dict_macror[i])*100))
        fw.write('\n')
        

def main():
    opt = parse_args()
    file_train=opt.input_tgt_train_file
    file_test=opt.input_tgt_test_file
    file_pred=opt.input_pred_file
    output_file=opt.output_file
    best_num=int(opt.best_number)
    n=int(opt.augmentation_number)
    calculate_metrics(file_train,file_test,file_pred,output_file,best_num,n)
        

if __name__ == '__main__':
    main()      


