clear;clc;close all

PES=fillinviasymmetry(diffusionPES('17_1vasp.dat'),'c3v');
PES.plotMEP2(  {[129,246]},     {[129,246.1]})