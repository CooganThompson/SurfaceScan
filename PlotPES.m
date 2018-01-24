clear; clc; close all
%% Pick the Surface Scan
% %17_1
% PESFull.Class=fillinviasymmetry(diffusionPES('17_1_vasp.dat'),'c3v');
%PES.plotGoogleMapsView

% %17_1 Co Old
% PES=diffusionPES('17_1_Co_OLDINCAR_vasp.dat');
% PES.plotGoogleMapsView

% %Still has some things that I don't know what are wrong with
% %17_1 Co 
% PES=diffusionPES('17_1_Co_vasp.dat');
% PES.plotGoogleMapsView

% %17_2
% PES=diffusionPES('17_2vasp.dat');
% PES.plotGoogleMapsView

% %17_2 Co Old
% PES=diffusionPES('17_2_Co_OLDINCAR_vasp.dat');
% PES.plotGoogleMapsView

% %17_2 Co 
% PES=diffusionPES('17_2_Co_vasp.dat');
% PES.plotGoogleMapsView

% %17_3  
% PES=fillinviasymmetry(diffusionPES('17_3vasp.dat'),'cs');
% PES.plotGoogleMapsView

% %17_3 Co
% PES=diffusionPES('17_3_Co_vasp.dat');
% PES.plotGoogleMapsView

% %18_1
% PES=fillinviasymmetry(diffusionPES('18_1_vasp.dat'),'c2v');
% PES.plotGoogleMapsView

% %18_1 Co Old
%  PES=diffusionPES('18_1_Co_vasp.dat');
% PES.plotGoogleMapsView

% %18_2
% PESFull.Class=diffusionPES('18_2_vasp.dat');
% PES.plotGoogleMapsView
% PESFull.Symmetry='Sn';

% %18_2 Co Old
% PES=diffusionPES('18_2_Co_vasp.dat');
% PES.plotGoogleMapsView

% %18_3  
% PESFull.Class=diffusionPES('18_3_vasp.dat');
% PESFull.Class.plotGoogleMapsView

% %18_3 Co
% PES=diffusionPES('18_3_Co_vasp.dat');
% PES.plotGoogleMapsView

% %19_1
% PESFull.Class=fillinviasymmetry(diffusionPES('19_1_vasp.dat'),'c3v');
% PES.plotGoogleMapsView
% PESFull.Symmetry='C3v';

% %19_1 Co 
% PES=diffusionPES('19_1_Co_vasp.dat');
% PES.plotGoogleMapsView

% %19_2
% PESFull.Class=diffusionPES('19_2_vasp.dat');
% PESFull.Symmetry='Sn';

% %19_2 Co
% PES=diffusionPES('19_2_Co_vasp.dat');
% PES.plotGoogleMapsView

%20_1
PESFull.Class=fillinviasymmetry(diffusionPES('20_1_vasp.dat'),'c3v');
PESFull.Class.plotGoogleMapsView
PESFull.Symmetry='C3v';

% %20_1 Co
% PES=diffusionPES('20_1_Co_vasp.dat');
% PES.plotGoogleMapsView

% %20_2
% PESFull.Class=fillinviasymmetry(diffusionPES('20_2_vasp.dat'),'cs');
% PESFull.Symmetry='Cs';
% PES.plotGoogleMapsView

% %20_2 Co
% PES=diffusionPES('20_2_Co_vasp.dat');
% PES.plotGoogleMapsView


%% Add other details
% %Diffusion Part
% PES.plotMEP2(  {[40,100]},     {[40,280]})
%PES.plotMEP2(  {[MinsLoc(4,:)]},{[MinsLoc(10,:)]})
%
%PESFull.Class=PES;

PESFull.Mins=GetMins(PESFull.Class);

PESFull.Class.plotContourMapsView
PlotMinsOnPES(PESFull);
PESFull.Barriers=GetDiffusionParameters(PESFull);

%to remove mins
%PESFull.Mins(1,:)=[]

PESFull.Class.plotContourMapsView;
PlotMinsOnPES(PESFull);
SpiderWebPlot(PESFull);

PlotTSEnergiesOnPES(PESFull)

MakeZacrosInputForDiffusions(PESFull.Class,PESFull.Mins,PESFull.Barriers)

% save('PESFull_Au18_1',PESFull)
% load('PESFull_Au18_1',PESFull)


%KMCresults=LoadSurfaceKMC(ZacrosOutputName) File needs to be in data folder
% KMCresults=LoadSurfaceKMC('18_1_Diffusion_KMC_Trial.txt');

%Load all the neccessary files
% load('data/PESFull_Au18_1.mat');

%AddingAtomDots(run from bash window)

%Plot Mins with Lines and Mins
PESFull.Class.plotContourMapsView;
PESFull.Class.plotAtomRegionsOverlay(PESFull.AlignmentData{1},PESFull.AlignmentData{2});
PlotMinsOnPES(PESFull);

%KMC simulation stuff
PESFull.Results=LoadSurfaceKMC('SAMPLERIGHT'); %20_2_A_procstat_output.txt

%PlotBarReactions(KMCresults);
PlotMapReactions(PESFull)

PESFull.Class.plot3DRepresentation(PESFull.AlignmentData{1},PESFull.AlignmentData{2});

%ShowReactionMap
PlotMinsOnPES(PESFull)
PlotMapReactions(PESFull)

PlotBarReactions(PESFull.Results)
PlotBarReactions(PESFull.Results,'log')

