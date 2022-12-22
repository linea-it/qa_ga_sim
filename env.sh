#!/bin/bash

export CONDAPATH=/home/adriano.pieres/miniconda3/bin
export GA_ROOT=/lustre/t1/cl/lsst/gawa_project/adriano.pieres/qa_ga_sim/qa_ga_sim
export PYTHONPATH=$PYTHONPATH:$GA_ROOT

source $CONDAPATH/activate
conda activate ga_sim
python /lustre/t1/cl/lsst/gawa_project/adriano.pieres/qa_ga_sim/qa_ga_sim.py
