function Barriers=AddToBarriers(Barriers,MEPObject,TSEnergy,FirstPoint,SecondPoint,Class,Mins)

%check to see if there is already a barrier there
if isnan(Barriers.Grid{FirstPoint,SecondPoint})
    AtoBcount=1;
else
    AtoBcount=length(Barriers.Grid{FirstPoint,SecondPoint})+1;
end

Barriers.Grid{FirstPoint,SecondPoint}(AtoBcount)=TSEnergy;
Barriers.Grid{SecondPoint,FirstPoint}(AtoBcount)=TSEnergy;

if isempty(Barriers.Name{1}) %the first one has to be difficult sometimes
    
    Barriers.Pathway{1}=[MEPObject.XCoords MEPObject.YCoords];
    Barriers.Name(1)={sprintf('%.0f_TO_%.0f',FirstPoint,SecondPoint)};
    Barriers.StartMinIndex{1}=FirstPoint;
    Barriers.EndMinIndex{1}=SecondPoint;
    Barriers.Energies{1}=[getMEPEnergies(Class,Mins(FirstPoint,1),Mins(FirstPoint,2)) TSEnergy getMEPEnergies(Class,Mins(SecondPoint,1),Mins(SecondPoint,2))];
    Barriers.TSLocaton{1}=[MEPObject.TS];
    
else
    count=FindCorrectIndexForPathway(Barriers,FirstPoint,SecondPoint);
    Barriers.Pathway=[ { Barriers.Pathway{1:count-1}} [MEPObject.XCoords MEPObject.YCoords] { Barriers.Pathway{count:end}}];
    if AtoBcount==1
        Barriers.Name=[ { Barriers.Name{1:count-1}} sprintf('%.0f_TO_%.0f',FirstPoint,SecondPoint) { Barriers.Name{count:end}}];
    elseif AtoBcount>1
        Barriers.Name=[ { Barriers.Name{1:count-1}} sprintf('%.0f_TO_%.0f_%.0f',FirstPoint,SecondPoint,AtoBcount) { Barriers.Name{count:end}}];
    else
        disp('AtoBcount error')
    end
    Barriers.StartMinIndex=[ { Barriers.StartMinIndex{1:count-1}} FirstPoint { Barriers.StartMinIndex{count:end}}];
    Barriers.EndMinIndex=[ { Barriers.EndMinIndex{1:count-1}} SecondPoint { Barriers.EndMinIndex{count:end}}];
    EnergiesTemp=[getMEPEnergies(Class,Mins(FirstPoint,1),Mins(FirstPoint,2)) TSEnergy getMEPEnergies(Class,Mins(SecondPoint,1),Mins(SecondPoint,2))];
    Barriers.Energies=[ { Barriers.Energies{1:count-1}} EnergiesTemp { Barriers.Energies{count:end}}];
    Barriers.TSLocaton=[ { Barriers.TSLocaton{1:count-1}} [MEPObject.TS] { Barriers.TSLocaton{count:end}}];
end
end
% Barriers.Pathway{count}=[MEPObject.XCoords MEPObject.YCoords];
% if AtoBcount>1
%     Barriers.Name{count}=sprintf('%.0f_TO_%.0f_%0.f',i,j,AtoBcount);
% elseif AtoBcount==1
%     Barriers.Name{count}=sprintf('%.0f_TO_%.0f',i,j);
% else
%     disp('Barriers.Name error')
% end
% Barriers.StartMinIndex{count}=i;
% Barriers.EndMinIndex{count}=j;
%
% Barriers.Energies{count}=[getMEPEnergies(PES,MinsLoc(i,1),MinsLoc(i,2)) TSEnergy getMEPEnergies(PES,MinsLoc(j,1),MinsLoc(j,2))];