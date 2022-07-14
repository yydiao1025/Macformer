# -*- coding: utf-8 -*-
'''
Created on Fri May 21 10:43:22 2021

@author: user
'''

import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-input_file', '-i', required=True,
                        help='file of acyclic-macrocyclic SMILES pairs')
    parser.add_argument('-output_path', '-o', required=True,
                        help='path of the generated src/tgt file')

    return parser.parse_args()


def read_pairs(file):
    pairs=[line.strip().split()[0] for line in open(file).readlines()]
    return pairs

def spliting(pairs):
    
#   split the data into train, validate, and test datasets with the ratio of 8:1:1
    
    length = [pair.split('.')[0] for pair in pairs]
    groups = list(set(length))
    data = list(zip(length, pairs))
    df1 = pd.DataFrame(data, columns=['length', 'pairs'])
    df2 = df1.groupby('length')

    test = []
    val = []
    train = []
    for g in groups:
        idx = list(df2.groups[g])

        l = len(idx)
        start = int(l * 0.1)
        end = int(l * 0.2)
        idx_test = idx[: start]
        idx_val = idx[start: end]
        idx_train = idx[end:]
        test = test + idx_test
        val = val + idx_val
        train = train + idx_train

    train_pairs = list(df1.iloc[train]['pairs'])
    val_pairs = list(df1.iloc[val]['pairs'])
    test_pairs = list(df1.iloc[test]['pairs'])

    return train_pairs, val_pairs, test_pairs

def write_src_tgt(pairs, name, dir):
    src_file = dir + 'src' + '-' + name
    tgt_file = dir + 'tgt' + '-' + name
    src = [pair.split('.')[0] + ' ' + ' '.join(pair.split('>')[0].split('.')[2]) for pair in pairs]
    tgt = [' '.join(pair.split('>')[1]) for pair in pairs]
    with open(src_file, 'w') as w:
        for j in src:
            w.write(j)
            w.write('\n')
    w.close()

    with open(tgt_file, 'w') as w2:
        for j in tgt:
            w2.write(j)
            w2.write('\n')
    w2.close()


def main():
    opt = parse_args()

    pairs = read_pairs(opt.input_file)
    train_pairs, val_pairs, test_pairs = spliting(pairs)
    write_src_tgt(train_pairs, name='train', dir = opt.output_path)
    write_src_tgt(val_pairs, name='val', dir = opt.output_path)
    write_src_tgt(test_pairs, name='test', dir = opt.output_path)


if __name__ == '__main__':
    main()
