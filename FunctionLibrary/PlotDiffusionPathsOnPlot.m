function PlotDiffusionPathsOnPlot(PESFull,KMCresults)

PES=fillinviasymmetry(diffusionPES('18_1vasp.dat'),'c2v');
Barriers=GetDiffusionParameters(PESFull.Class,PESFull.Mins);




end