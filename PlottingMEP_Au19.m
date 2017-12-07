clear;clc;close all

%PES=fillinviasymmetry(diffusionPES('18_1vasp.dat'),'c2v');
%PES=diffusionPES('18_2vasp.dat');
%PES=diffusionPES('18_3vasp.dat');
%PES=fillinviasymmetry(diffusionPES('19_1vasp.dat'),'c3v');
%PES=diffusionPES('19_2vasp.dat');
%PES=fillinviasymmetry(diffusionPES('20_2vasp.dat'),'c3v');
PES=diffusionPES('20_2vasp.dat');

PES.plotMEP2(  {[60,80]},     {[100,30]})