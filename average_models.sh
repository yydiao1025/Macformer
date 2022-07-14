#!/usr/bin/env bash

dataset_name=data/ChEMBL/a10

datadir=checkpoints/ChEMBL_a10/

model=$datadir/${dataset_name}_model_average.pt

python average_models.py -p $datadir -o $model
