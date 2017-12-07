function HistoryData = RetrieveHistoryData(ZacrosOutput)
%ZacrosOutput=('17_1_CN5_history_output.txt');
fid=fopen(['data/' ZacrosOutput],'r');

GasSpecies=textscan(fid,'%s',1,'Delimiter','\n');
GasSpecies=char(GasSpecies{1});
GasSpecies=strsplit(GasSpecies,' ');
GasSpecies(1)=[];

SurfaceSpecies=textscan(fid,'%s',1,'Delimiter','\n');
SurfaceSpecies=char(SurfaceSpecies{1});
SurfaceSpecies=strsplit(SurfaceSpecies,' ');
SurfaceSpecies(1)=[];

SiteTypes=textscan(fid,'%s',1,'Delimiter','\n');
SiteTypes=char(SiteTypes{1});
SiteTypes=strsplit(SiteTypes);
SiteTypes(1)=[];

line=textscan(fid,'%s',1,'Delimiter','\n');


test=~isempty(line{1});
i=1;
while test
    
    line=char(line{1});
    line=strsplit(line);
    
    ConfigsWritten(i)   =str2double(cell2mat(line(2)));
    KMCEventsHappened(i)=str2double(cell2mat(line(3)));
    Time(i)             =str2double(cell2mat(line(4)));
    Temperature(i)      =str2double(cell2mat(line(5)));
    Energy(i)           =str2double(cell2mat(line(6)));
    
    for j=1:length(SiteTypes)
        temp=textscan(fid,'%s',1,'Delimiter','\n');
        temp=char(temp{1});
        temp=strsplit(temp);
        temp=str2double(temp);
        data(i,j,:)=temp;   %indices are time,row,column
    end
    
    line=textscan(fid,'%s',1,'Delimiter','\n');
    line=char(line{1});
    line=strsplit(line);
    MoleculesProduced(i)=str2double(cell2mat(line));
    
    
    line=textscan(fid,'%s',1,'Delimiter','\n');
    if ~isempty(line{1})
        i=i+1
    else
        test=0;
    end
end

HistoryData.GasSpecies=GasSpecies';
HistoryData.SurfaceSpecies=SurfaceSpecies';
HistoryData.ConfigsWritten=ConfigsWritten;
HistoryData.KMCEventsHappened=KMCEventsHappened;
HistoryData.Time=Time;
HistoryData.Temperature=Temperature;
HistoryData.Energy=Energy;
HistoryData.PositionData=data;
HistoryData.MoleculesProduced=MoleculesProduced;

end