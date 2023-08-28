#!/bin/bash
#SBATCH -t 0-02:00
#SBATCH -c 1

make
time $HOME/projects/def-fdion/DDQ/DDQ_HHG_pull_2003-08-07/wp_DDQ_HHG_v01
