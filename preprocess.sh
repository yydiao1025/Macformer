#!/usr/bin/env bash


# dataset_name=ChEMBL_a10
python preprocess.py -train_src data/ChEMBL/a10/src-train-a10 \
                     -train_tgt data/ChEMBL/a10/tgt-train-a10 \
                     -valid_src data/ChEMBL/a10/src-val \
                     -valid_tgt data/ChEMBL/a10/tgt-val \
                     -save_data data/ChEMBL/a10/${dataset} \
                     -src_seq_length 500 -tgt_seq_length 500 \
                     -src_vocab_size 1000 -tgt_vocab_size 1000 -share_vocab
