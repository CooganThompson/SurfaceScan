function Barriers=GetDiffusionParameters(PESFull)
%Note this is not well tested as I left right as I added some new stuff
%Use the backup version, it is well tested
%Must put this into matrix form
%Currently only works for Cv symmetries

if nargin==1 %Allows you to just input the the Full PES class
    MinsLoc=PESFull.Mins;
    Symmetry=PESFull.Symmetry;
    PESFull=PESFull.Class;
elseif nargin==2 %if you don't give it a symmetry, you have to pay the price
    Symmetry='C1';
end

if strcmpi(Symmetry,'Cs')||strcmpi(Symmetry,'cs') %easier than making it it's own class
    Symmetry='C1v';
end

children = get(gca, 'children'); %gets the number of things that should be there so I can delete the rest later
NumOfPlotThings=length(children);
if length(children)~=length(MinsLoc)+2
    warning('The plot this is on might produce unusual results')
    Answer=input('Would you like me to try to clean the plot y/n','s');
    if Answer=='y'||Answer=='Y'
        NumOfPlotThings=length(MinsLoc)+2;
        CleanPlot(NumOfPlotThings)
    end
end


PESFull=PESFull.interpolatePESGrid;
MEPObject=minimumEnergyPathway(PESFull.alphaGrid,PESFull.betaGrid,PESFull.energiesGrid);

%initalizations
Barriers.Grid=cell(length(MinsLoc));
Barriers.Pathway=cell(1);
Barriers.Name=cell(1);
Barriers.StartMinIndex=cell(1);
Barriers.EndMinIndex=cell(1);
Barriers.Energies=cell(1);
Barriers.TSLocaton=cell(1);

skipto=NaN;
count=0;
secondentryflag=0;
FirstPoint=0;
while FirstPoint<length(MinsLoc) %loop over all the points
    FirstPoint=FirstPoint+1;
    
    skipflag=0;
    
    SecondPoint=FirstPoint;
    while SecondPoint<length(MinsLoc) %loop over all the points
        if skipto==SecondPoint+1
            skipflag=0;
            skipto=NaN;
        end
        
        SecondPoint=SecondPoint+1;
        
        AtoBcount=0;
        
        X1=sind(MinsLoc(FirstPoint,1)).*cosd(MinsLoc(FirstPoint,2));
        X2=sind(MinsLoc(SecondPoint,1)).*cosd(MinsLoc(SecondPoint,2));
        Y1=sind(MinsLoc(FirstPoint,1)).*sind(MinsLoc(FirstPoint,2));
        Y2=sind(MinsLoc(SecondPoint,1)).*sind(MinsLoc(SecondPoint,2));
        Z1=cosd(MinsLoc(FirstPoint,1));
        Z2=cosd(MinsLoc(SecondPoint,1));
        dist=sqrt((X1-X2)^2+(Y1-Y2)^2+(Z1-Z2)^2);
        
        BoxIsEmptyFlag=isempty(Barriers.Grid{FirstPoint,SecondPoint});
        
        if BoxIsEmptyFlag||secondentryflag%dist<1.2 %~skipflag && BoxIsEmptyFlag%dist<1.2 %mindistnace on perfect sphere/radius
            
            
            if ~skipflag
                try
                    MEPObject=MEPObject.calculateMEPMultiTS({MinsLoc(FirstPoint,:)},{MinsLoc(SecondPoint,:)});
                    [NumberTSs,~]=size(MEPObject.TS);
                    [NumberMins,~]=size(MEPObject.Minima);
                    
                    PlotMEPAllSymmetry(MEPObject,Symmetry)
                    
                    [TScount,~]=size(MEPObject.TS);
                    TSEnergies=0;
                    for k=1:TScount
                        TSEnergies(k)=getMEPEnergies(PESFull,MEPObject.TS(k,1),MEPObject.TS(k,2));
                    end
                    TSEnergy=max(TSEnergies);
                    
                catch
                    NumberTSs=NaN;
                    NumberMins=NaN;
                    TSEnergy=NaN;
                    disp('Maybe try drawing it yourself?')
                end
            end
            exitflag=0;
            if skipflag
                exitflag=1;
                accepted='n';
            end
            
            while exitflag==0
                disp('-----------------------------------------')
                disp('    Start     End       Distance  TS        Mins    HighestTS')
                disp([FirstPoint     SecondPoint   dist     NumberTSs NumberMins TSEnergy])
                accepted=input('Is this an acceptable MEP? y/n (or O for options)','s');
                try
                    if length(accepted)>=2&&(accepted(1)=='S'||accepted(1)=='s')
                        skipflag=1;
                        skipto=str2double(accepted(3:end));
                        exitflag=1;
                        accepted='n';
                    elseif accepted=='EXIT'
                        return
                    elseif length(accepted)~=1
                        disp('Please enter one and only one entry')
                    elseif accepted=='y'||accepted=='Y'||accepted=='n'||accepted=='N'
                        exitflag=1;
                        secondentryflag=0;
                        AtoBcount=AtoBcount+1;
                    elseif accepted=='O'||accepted=='o'
                        disp('Options')
                        disp('"Save" and exit, EXIT')
                        %                     disp('Remove line, r')
                        disp('Skip to Next Starting Point, s')
                        disp('Use a Different Line, d')
                        disp('Ues this one and include another path between the mins, a')
                    elseif accepted=='a'||accepted=='A'
                        exitflag=1;
                        secondentryflag=1;
                        AtoBcount=AtoBcount+1;
                        
                        %                 elseif accepted=='R'||accepted=='r'
                        %                     children = get(gca, 'children');
                        %                     delete(children(2));
                    elseif accepted=='S'||accepted=='s'
                        skipflag=1;
                        exitflag=1;
                        accepted='n';
                    elseif accepted=='D'||accepted=='d'
                        
                        points=[];
                        while isempty(points)
                            try
                                points=input('What are the points [alpha1 beta1 alpha2 beta2]');
                            catch
                                disp('Not a valid entry')
                            end
                        end
                        
                        CleanPlot(NumOfPlotThings)
                        
                        try
                            MEPObject=MEPObject.calculateMEPMultiTS({points([1,2])},{points([3,4])});
                        catch
                            disp('Try different points')
                        end
                        
                        [NumberTSs,~]=size(MEPObject.TS);
                        [NumberMins,~]=size(MEPObject.Minima);
                        
                        PlotMEPAllSymmetry(MEPObject,Symmetry)
                        
                        [TScount,~]=size(MEPObject.TS);
                        TSEnergies=0;
                        for k=1:TScount
                            TSEnergies(k)=getMEPEnergies(PESFull,MEPObject.TS(k,1),MEPObject.TS(k,2));
                        end
                        TSEnergy=max(TSEnergies);
                        
                    else
                        disp('Invalied')
                        
                    end
                catch
                end
                
            end
            
            if accepted=='y'||accepted=='Y'||accepted=='a'||accepted=='A'
                
                count=count+1;
                %Barriers=AddToBarriers(Barriers,MEPObject,TSEnergy,FirstPoint,SecondPoint,PES,MinsLoc);
                
                Points=GetAllSymmetryPoints(MinsLoc,FirstPoint,SecondPoint,Symmetry);
                Points=Cleanpoints(Points);
                
                MEPObjectSSS=GiveSymmetryVersions(MEPObject,Symmetry);
                
                [Barrierslength,~]=size(Points);
                for l=1:Barrierslength%length(MEPObjectSSS)
                    Barriers=AddToBarriers(Barriers,MEPObjectSSS{l},TSEnergy,Points(l,1),Points(l,4),PESFull,MinsLoc);
                end
                
                CleanPlot(NumOfPlotThings)
                
            else
                
                Points=GetAllSymmetryPoints(MinsLoc,FirstPoint,SecondPoint,Symmetry);
                Points=Cleanpoints(Points);
                
                [Barrierslength,~]=size(Points);
                for l=1:Barrierslength
                    Barriers.Grid{Points(l,1),Points(l,4)}=NaN;
                    Barriers.Grid{Points(l,4),Points(l,1)}=NaN;
                end
                CleanPlot(NumOfPlotThings)
                
            end
            %         else
            %             Barriers.Grid{FirstPoint,SecondPoint}=NaN;
            %             Barriers.Grid{SecondPoint,FirstPoint}=NaN;
        end
        if secondentryflag==1
            SecondPoint=SecondPoint-1;
        end
        
    end
end
end