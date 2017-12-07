function KMCresults=LoadSurfaceKMC(ZacrosOutput)
%Reads in the results of the KMC run
%fid=fopen('data/18_1_Diffusion_KMC_Trial.txt','r'); %debugging lines
%ZacrosOutput='18_1_Diffusion_KMC_Trial.txt'
fid=fopen(['data/' ZacrosOutput],'r');

NamesAll=textscan(fid,'%s',1,'Delimiter','\n');

NamesAll=char(NamesAll{1});
Names=strsplit(NamesAll,'  ');

KMCresults.Names=Names(2:end);

for i=2:length(Names)
    name=cell2mat(Names(i));
    Locations(:,i-1)=sscanf(name,' diff/H/%d_TO_%d');
end
KMCresults.Locations=Locations;

test=1;
counter=0;

while test
    configuration=textscan(fid,'%s',1,'Delimiter','\n');
    residentstimes=textscan(fid,'%s',1,'Delimiter','\n');
    events=textscan(fid,'%s',1,'Delimiter','\n');
    [test,~]=size(configuration{1});
    counter=counter+1;
    
    if 1-test
        break
    end
    
    configuration=char(configuration{1});
    configuration=strsplit(configuration,'  ') ;
    
    KMCresults.Configs(counter,1)=str2num(cell2mat(configuration(2)));
    KMCresults.Steps(counter,1)=str2num(cell2mat(configuration(3)));
    KMCresults.Time(counter,1)=str2num(cell2mat(configuration(4)));
    
    residentstimes=char(residentstimes{1});
    residentstimes=strsplit(residentstimes,'  ') ;
    
    KMCresults.ResidenceTimesAll(counter,1)=str2num(cell2mat(residentstimes(1)));
    temp=str2num(cell2mat(residentstimes([2:end])));
    KMCresults.ResidenceTimesIndividual(counter,:)=temp;
    
    events=char(events{1});
    events=strsplit(events,'  ') ;
    
    KMCresults.EventsAll(counter,1)=str2num(cell2mat(events(1)));
    temp2=str2double(cellstr(events([2:end])));
    KMCresults.EventsIndividual(counter,:)=temp2;
end

