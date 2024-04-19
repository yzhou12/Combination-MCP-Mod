#!/bin/bash
#BSUB -e log_%J_%I.err
#BSUB -o log_%J_%I
#BSUB -n 1
#BSUB -W 1:00
Rscript connect_results.R $LSB_JOBINDEX