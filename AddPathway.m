function PESFull=AddPathway(PESFull)

PESFull.Class=PESFull.Class.interpolatePESGrid;
MEPObject=minimumEnergyPathway(PESFull.Class.alphaGrid,PESFull.Class.betaGrid,PESFull.Class.energiesGrid);

FirstPoint=input('Please input the number of the first minimum \n');
SecondPoint=input('Please input the number of the second minimum \n');

%Make sure inputs are valid
checkflag=0;
while checkflag==0
    
    if FirstPoint~=floor(FirstPoint)||FirstPoint>length(PESFull.Mins)
        disp('The first entry is invalid \n')
        FirstPoint=input('Please input the number of the first minimum \n');
        
    elseif SecondPoint~=floor(SecondPoint)||SecondPoint>length(PESFull.Mins)
        disp('These are not integers')
        SecondPoint=input('Please input the number of the second minimum \n');
        
    elseif FirstPoint==SecondPoint
        disp('These are the same point \n')
        SecondPoint=input('Please input the number of the second minimum \n');
        
    elseif FirstPoint<SecondPoint
        checkflag=1;
        
    elseif FirstPoint>SecondPoint
        temp=FirstPoint;
        FirstPoint=SecondPoint;
        SecondPoint=temp;
        checkflag=1;
        
    end
    
end

disp('-----------------------------------------')
try
    MEPObject=MEPObject.calculateMEPMultiTS({PESFull.Mins(FirstPoint,:)},{PESFull.Mins(SecondPoint,:)});
    [NumberTSs,~]=size(MEPObject.TS);
    [NumberMins,~]=size(MEPObject.Minima);
    MEPObject.plotMEPLines
    drawnow
catch
    NumberTSs=0;
    NumberMins=0;
    disp('Maybe try drawing it yourself?')
end

X1=sind(PESFull.Mins(FirstPoint,1)).*cosd(PESFull.Mins(FirstPoint,2));
X2=sind(PESFull.Mins(SecondPoint,1)).*cosd(PESFull.Mins(SecondPoint,2));
Y1=sind(PESFull.Mins(FirstPoint,1)).*sind(PESFull.Mins(FirstPoint,2));
Y2=sind(PESFull.Mins(SecondPoint,1)).*sind(PESFull.Mins(SecondPoint,2));
Z1=cosd(PESFull.Mins(FirstPoint,1));
Z2=cosd(PESFull.Mins(SecondPoint,1));
dist=sqrt((X1-X2)^2+(Y1-Y2)^2+(Z1-Z2)^2);

exitflag=0;
while exitflag==0
    disp('    Start     End       Distance  TS        Mins')
    disp([FirstPoint     SecondPoint   dist     NumberTSs NumberMins])
    accepted=input('Is this an acceptable MEP? y/n (or O for options)','s');
    
    if length(accepted)~=1
        disp('Please enter one and only one entry')
        
    elseif accepted=='y'||accepted=='Y'||accepted=='n'||accepted=='N'
        exitflag=1;
        
    elseif accepted=='O'||accepted=='o'
        disp('Options')
        disp('Remove line, r')
        disp('Use a Different Line, d')
        
    elseif accepted=='R'||accepted=='r'
        children = get(gca, 'children');
        delete(children(2));
    elseif accepted=='D'||accepted=='d'
        
        points=[];
        while isempty(points)
            try
                points=input('What are the points [alpha1 beta1 alpha2 beta2]');
            catch
                disp('Not a valid entry')
            end
        end
        
        children = get(gca, 'children');
        delete(children(1));
        try
            MEPObject=MEPObject.calculateMEPMultiTS({points([1,2])},{points([3,4])});
        catch
            disp('Try different points')
        end
        [NumberTSs,~]=size(MEPObject.TS);
        [NumberMins,~]=size(MEPObject.Minima);
        MEPObject.plotMEPLines
        drawnow
        
    else
        disp('Invalied')
    end
    
end

if accepted=='y'||accepted=='Y'
    
    [TScount,~]=size(MEPObject.TS);
    TSEnergies=ones(1,TScount)*nan;
    for k=1:TScount
        TSEnergies(k)=getMEPEnergies(PESFull.Class,MEPObject.TS(k,1),MEPObject.TS(k,2));
    end
    TSEnergy=max(TSEnergies);
    
    if isnan(PESFull.Barriers.Grid{FirstPoint,SecondPoint})
        AtoBcount=1;
    else
        AtoBcount=length(PESFull.Barriers.Grid{FirstPoint,SecondPoint})+1;
    end
    
    PESFull.Barriers.Grid{FirstPoint,SecondPoint}(AtoBcount)=TSEnergy;
    PESFull.Barriers.Grid{SecondPoint,FirstPoint}(AtoBcount)=TSEnergy;
    
    %Finds out where the new path should go
    foundStartSpot=0;
    foundEndSpot=0;
    
    count=0;
    while foundStartSpot==0
        count=count+1;
        
        if count==length(PESFull.Barriers.StartMinIndex)
            count=count+1;
            foundStartSpot=1;
            foundEndSpot=1;
        elseif PESFull.Barriers.StartMinIndex{count}==FirstPoint
            foundStartSpot=1;
        elseif PESFull.Barriers.StartMinIndex{count}>=FirstPoint
            foundStartSpot=1;
            foundEndSpot=1;
        end
    end
    
    while foundEndSpot==0
        count=count+1;
        
        if PESFull.Barriers.EndMinIndex{count}>=SecondPoint
            foundEndSpot=1;
        end
        
    end
    
    PESFull.Barriers.Pathway=[ { PESFull.Barriers.Pathway{1:count-1}} [MEPObject.XCoords MEPObject.YCoords] { PESFull.Barriers.Pathway{count:end}}];
    if AtoBcount==1
        PESFull.Barriers.Name=[ { PESFull.Barriers.Name{1:count-1}} sprintf('%.0f_TO_%.0f',FirstPoint,SecondPoint) { PESFull.Barriers.Name{count:end}}];
    elseif AtoBcount>1
        PESFull.Barriers.Name=[ { PESFull.Barriers.Name{1:count-1}} sprintf('%.0f_TO_%.0f_%.0f',FirstPoint,SecondPoint,AtoBcount) { PESFull.Barriers.Name{count:end}}];
    else
        disp('AtoBcount error')
    end
    
    PESFull.Barriers.StartMinIndex=[ { PESFull.Barriers.StartMinIndex{1:count-1}} FirstPoint { PESFull.Barriers.StartMinIndex{count:end}}];
    PESFull.Barriers.EndMinIndex=[ { PESFull.Barriers.EndMinIndex{1:count-1}} SecondPoint { PESFull.Barriers.EndMinIndex{count:end}}];
    EnergiesTemp=[getMEPEnergies(PESFull.Class,PESFull.Mins(FirstPoint,1),PESFull.Mins(FirstPoint,2)) TSEnergy getMEPEnergies(PESFull.Class,PESFull.Mins(SecondPoint,1),PESFull.Mins(SecondPoint,2))];
    PESFull.Barriers.Energies=[ { PESFull.Barriers.Energies{1:count-1}} EnergiesTemp { PESFull.Barriers.Energies{count:end}}];
    
    children = get(gca, 'children');
    delete(children(1));
end

end