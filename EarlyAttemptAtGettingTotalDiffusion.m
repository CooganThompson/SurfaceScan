clear;clc;close all

fID=fopen('Au17_1_H_Diffusion.dat');
text=textscan(fID,'%s %f32');
fclose(fID);

load Barriers_Au17_3_15.mat

DiffusionResults.Barriers = Barriers;

for i=2:length(text{1})
    temp=char(text{1}(i));
    DiffusionResults.Reaction.Name{i-1}=temp([8:end]);
    DiffusionResults.Reaction.Occurences{i-1}=text{2}(i);
    
    nametimep=char(DiffusionResults.Reaction.Name(i-1))
    ndex=strcmp(Barriers.Name,)
    
    %DiffusionResults.Reaction.Distance
end