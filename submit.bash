#!/bin/bash
#SBATCH -t 0-02:00
#SBATCH -c 1

make
time ./wp_DDQ_HHG_v01
