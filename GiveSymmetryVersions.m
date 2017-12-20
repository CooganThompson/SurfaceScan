function MEPObjectSSS=GiveSymmetryVersions(MEPObject,Symmetry)
%Currently works for just Cn symmterys

Cnumber=str2num(Symmetry(2));

for i=1:Cnumber
    MEPObjectSSS{i}=MEPObject;
    
    MEPObjectSSS{i}.YCoords=mod(MEPObject.YCoords+(i-1)*360/Cnumber,360);
    MEPObjectSSS{i}.Minima(:,2)=mod(MEPObject.Minima(:,2)+(i-1)*360/Cnumber,360);
    MEPObjectSSS{i}.TS(:,2)=mod(MEPObject.TS(:,2)+(i-1)*360/Cnumber,360);
end

CurrentSyms=length( MEPObjectSSS);
try
    if Symmetry(3)=='v'
        for i=CurrentSyms+1:2*CurrentSyms
            
            MEPObjectSSS{i}=MEPObjectSSS{i-CurrentSyms};
            
            MEPObjectSSS{i}.YCoords=360-MEPObjectSSS{i}.YCoords;
            MEPObjectSSS{i}.Minima(:,2)=360-MEPObjectSSS{i}.Minima(:,2);
            MEPObjectSSS{i}.TS(2)=360-MEPObjectSSS{i}.TS(2);
        end
    end
catch
end

CurrentSyms=length( MEPObjectSSS);
try
    if Symmetry(3)=='h'
        for i=CurrentSyms+1:2*CurrentSyms
            
            MEPObjectSSS{i}=MEPObjectSSS{i-CurrentSyms};
            
            MEPObjectSSS{i}.XCoords=180-MEPObjectSSS{i}.XCoords;
            MEPObjectSSS{i}.Minima(:,1)=180-MEPObjectSSS{i}.Minima(:,1);
            MEPObjectSSS{i}.TS(1)=180-MEPObjectSSS{i}.TS(1);
        end
    end
catch
end

end
