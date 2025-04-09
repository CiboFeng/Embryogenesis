#!/bin/bash
./Obs_Exp_P.x ../PC_100000_iced_chr14_dense.matrix_P.dat 577
python ice_do.py
octave Po_e.m
