#!/bin/bash
#SBATCH -t 0-02:00
#SBATCH -c 1

make
time ./travail
