function Barriers=GetDiffusionParameters(PES,MinsLoc)

%Must put this into matrix form
PES=PES.interpolatePESGrid;
MEPObject=minimumEnergyPathway(PES.alphaGrid,PES.betaGrid,PES.energiesGrid);

Barriers.Grid=cell(length(MinsLoc));
count=0;
for i=1:length(MinsLoc)
    for j=i+1:length(MinsLoc)
        X1=sind(MinsLoc(i,1)).*cosd(MinsLoc(i,2));
        X2=sind(MinsLoc(j,1)).*cosd(MinsLoc(j,2));
        Y1=sind(MinsLoc(i,1)).*sind(MinsLoc(i,2));
        Y2=sind(MinsLoc(j,1)).*sind(MinsLoc(j,2));
        Z1=cosd(MinsLoc(i,1));
        Z2=cosd(MinsLoc(j,1));
        dist=(X1-X2)^2+(Y1-Y2)^2+(Z1-Z2).^2;
        if dist<2 %mindistnace on perfect sphere/radius
            %MEPObject=MEPObject.calculateMEP({MinsLoc(i,:)},{MinsLoc(j,:)});
            MEPObject=MEPObject.calculateMEPMultiTS({MinsLoc(i,:)},{MinsLoc(j,:)});
               %doesn't work?         %If the TS is less than the ends, moslty because it just bad
            if isempty(MEPObject.TS)||getMEPEnergies(PES,MEPObject.TS(1),MEPObject.TS(2))<getMEPEnergies(PES,MinsLoc(i,1),MinsLoc(i,2))
                Barriers.Grid{i,j}=NaN;
                Barriers.Grid{j,i}=NaN;
            else
                count=count+1;
                
                Barriers.Grid{i,j}=MEPObject.TS;
                Barriers.Grid{j,i}=MEPObject.TS;
                
                Barriers.Name{count}=sprintf('%.0f_TO_%.0f',i,j);
                Barriers.StartMinIndex{count}=i;
                Barriers.EndMinIndex{count}=j;
                
                Barriers.Energies{count}=[getMEPEnergies(PES,MinsLoc(i,1),MinsLoc(i,2)) getMEPEnergies(PES,MEPObject.TS(1),MEPObject.TS(2)) getMEPEnergies(PES,MinsLoc(j,1),MinsLoc(j,2))];
            end
        else
            Barriers.Grid{i,j}=NaN;
            Barriers.Grid{j,i}=NaN;
            
        end
    end
end
end