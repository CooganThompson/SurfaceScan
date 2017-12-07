function Barriers=GetDiffusionParameters(PES,MinsLoc,Symmetry)
%Note this is not well tested as I left right as I added some new stuff
%Use the backup version, it is well tested
%Must put this into matrix form
%Currently only works for Cv symmetries

if nargin==1 %Allows you to just input the the Full PES class
    MinsLoc=PES.Mins;
    Symmetry=PES.Symmetry;
    PES=PES.Class;
elseif nargin==2 %if you don't give it a symmetry, you have to pay the price
    Symmetry='C1';
end

if strcmpi(Symmetry,'Cs')||strcmpi(Symmetry,'cs') %easier than making it it's own class
    symmetry='C1v';
end

PES=PES.interpolatePESGrid;
MEPObject=minimumEnergyPathway(PES.alphaGrid,PES.betaGrid,PES.energiesGrid);

Barriers.Grid=cell(length(MinsLoc));
count=0;
secondentryflag=0;
i=0;
while i<length(MinsLoc)
    i=i+1;
    
    skipflag=0;
    
    j=i;
    while j<length(MinsLoc)
        j=j+1;
        
        AtoBcount=0;
        
        X1=sind(MinsLoc(i,1)).*cosd(MinsLoc(i,2));
        X2=sind(MinsLoc(j,1)).*cosd(MinsLoc(j,2));
        Y1=sind(MinsLoc(i,1)).*sind(MinsLoc(i,2));
        Y2=sind(MinsLoc(j,1)).*sind(MinsLoc(j,2));
        Z1=cosd(MinsLoc(i,1));
        Z2=cosd(MinsLoc(j,1));
        dist=sqrt((X1-X2)^2+(Y1-Y2)^2+(Z1-Z2)^2);
        
        if ~skipflag%dist<1.2 %mindistnace on perfect sphere/radius
            disp('-----------------------------------------')
            
            try
                MEPObject=MEPObject.calculateMEPMultiTS({MinsLoc(i,:)},{MinsLoc(j,:)});
                [NumberTSs,~]=size(MEPObject.TS);
                [NumberMins,~]=size(MEPObject.Minima);
                MEPObjectSSS=GiveSymmetryVersions(MEPObject,Symmetry)
                
                for k=1:length(MEPObjectSSS)
                MEPObjectSSS{i}.plotMEPLines
                end
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
                    secondentryflag=0;
                    AtoBcount=AtoBcount+1;
                    
                    
                elseif accepted=='O'||accepted=='o'
                    disp('Options')
                    disp('Remove line, r')
                    disp('Skip to Next Starting Point, s')
                    disp('Use a Different Line, d')
                    disp('Ues this one and include another path between the mins, a')
                elseif accepted=='a'||accepted=='A'
                    exitflag=1;
                    secondentryflag=1;
                    AtoBcount=AtoBcount+1;
                    
                    
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
            
            if accepted=='y'||accepted=='Y'||accepted=='a'||accepted=='A'
                
                disp('counted')
                count=count+1;
                
                [TScount,~]=size(MEPObject.TS);
                for k=1:TScount
                    TSEnergies(k)=getMEPEnergies(PES,MEPObject.TS(k,1),MEPObject.TS(k,2));
                end
                TSEnergy=max(TSEnergies);
                
                Barriers.Grid{i,j}(AtoBcount)=TSEnergy;
                Barriers.Grid{j,i}(AtoBcount)=TSEnergy;
                
                Barriers.Pathway{count}=[MEPObject.XCoords MEPObject.YCoords];
                if AtoBcount>1
                    Barriers.Name{count}=sprintf('%.0f_TO_%.0f_%0.f',i,j,AtoBcount);
                elseif AtoBcount==1
                    Barriers.Name{count}=sprintf('%.0f_TO_%.0f',i,j);
                else
                    disp('Barriers.Name error')
                end
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
        if secondentryflag==1
            j=j-1;
        end
        
        
    end
end
end