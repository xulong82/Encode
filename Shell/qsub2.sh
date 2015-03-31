#!/bin/bash

qsub -v index=1 run_r_on_cluster.sh
qsub -v index=2 run_r_on_cluster.sh
qsub -v index=3 run_r_on_cluster.sh
qsub -v index=4 run_r_on_cluster.sh
qsub -v index=5 run_r_on_cluster.sh
qsub -v index=6 run_r_on_cluster.sh
qsub -v index=7 run_r_on_cluster.sh
qsub -v index=8 run_r_on_cluster.sh
qsub -v index=9 run_r_on_cluster.sh
qsub -v index=10 run_r_on_cluster.sh

