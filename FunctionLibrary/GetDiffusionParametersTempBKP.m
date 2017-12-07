function Barriers=GetDiffusionParameters(PES,MinsLoc)

%Must put this into matrix form
PES=PES.interpolatePESGrid;
MEPObject=minimumEnergyPathway(PES.alphaGrid,PES.betaGrid,PES.energiesGrid);

Barriers.Grid=cell(length(MinsLoc));
count=0;
for i=1:length(MinsLoc)-1
    skipflag=0;
    for j=i+1:length(MinsLoc)
        
        X1=sind(MinsLoc(i,1)).*cosd(MinsLoc(i,2));
        X2=sind(MinsLoc(j,1)).*cosd(MinsLoc(j,2));
        Y1=sind(MinsLoc(i,1)).*sind(MinsLoc(i,2));
        Y2=sind(MinsLoc(j,1)).*sind(MinsLoc(j,2));
        Z1=cosd(MinsLoc(i,1));
        Z2=cosd(MinsLoc(j,1));
        dist=sqrt((X1-X2)^2+(Y1-Y2)^2+(Z1-Z2)^2);

        if dist<1.2 %mindistnace on perfect sphere/radius
            disp('-----------------------------------------')
                        
            try
                MEPObject=MEPObject.calculateMEPMultiTS({MinsLoc(i,:)},{MinsLoc(j,:)});
                [NumberTSs,~]=size(MEPObject.TS);
                [NumberMins,~]=size(MEPObject.Minima);
                MEPObject.plotMEPLines
                drawnow
            catch
                NumberTSs=0;
                NumberMins=0;
                disp('Maybe try drawing it yourself?')
            end
            
            exitflag=0;
            if skipflag
                exitflag=1;
                accepted='n';
            end
            
            while exitflag==0
                disp('    Start     End       Distance  TS        Mins') 
                disp([i     j   dist     NumberTSs NumberMins])
                accepted=input('Is this an acceptable MEP? y/n (or O for options)','s');
                
                if length(accepted)~=1
                    disp('Please enter one and only one entry')
                elseif accepted=='y'||accepted=='Y'||accepted=='n'||accepted=='N'
                    exitflag=1;
                elseif accepted=='O'||accepted=='o'
                    disp('Options')
                    disp('Remove line, r')
                    disp('Skip to Next Starting Point, s')
                    disp('Use a Different Line, d')
                elseif accepted=='R'||accepted=='r'
                    children = get(gca, 'children');
                    delete(children(2));
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

                disp('counted')
                count=count+1;
                   
                [TScount,~]=size(MEPObject.TS);
                for k=1:TScount
                TSEnergies(k)=getMEPEnergies(PES,MEPObject.TS(k,1),MEPObject.TS(k,2));
                end
                TSEnergy=max(TSEnergies);
                
                Barriers.Grid{i,j}=TSEnergy;
                Barriers.Grid{j,i}=TSEnergy;
                Barriers.Pathway{count}=[MEPObject.XCoords MEPObject.YCoords];
                Barriers.Name{count}=sprintf('%.0f_TO_%.0f',i,j);
                Barriers.StartMinIndex{count}=i;
                Barriers.EndMinIndex{count}=j;
                
                Barriers.Energies{count}=[getMEPEnergies(PES,MinsLoc(i,1),MinsLoc(i,2)) TSEnergy getMEPEnergies(PES,MinsLoc(j,1),MinsLoc(j,2))];
                children = get(gca, 'children');
                delete(children(1));
            else
                Barriers.Grid{i,j}=NaN;
                Barriers.Grid{j,i}=NaN;
                
                children = get(gca, 'children');
                delete(children(1));
            end
        else
            Barriers.Grid{i,j}=NaN;
            Barriers.Grid{j,i}=NaN;   
        end
    end
end
end