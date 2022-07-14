#!/usr/bin/env bash

dataset_name=ChEMBL_a10
model=${dataset_name}_model_average.pt

python translate.py -model checkpoints/${dataset_name}/${model} \
                    -src data/ChEMBL/a10/src-test \
                    -output checkpoints/predictions_${model}_on_ChEMBL_test_best10.txt \
                    -batch_size 32 -replace_unk -max_length 500 -beam_size 10 -verbose -n_best 10 -gpu 1
