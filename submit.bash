#!/bin/bash
#SBATCH -t 0-12:00
#SBATCH -c 1

make
time ./travail
