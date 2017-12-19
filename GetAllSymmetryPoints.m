function [Points]=GetAllSymmetryPoints(MinsLoc,FirstPoint,SecondPoint,Symmetry)

Cnumber=str2num(Symmetry(2));     

%Number1 Location1(X) Location1(Y) Number2 Location2(X) Location2(Y)
Points=[FirstPoint MinsLoc(FirstPoint,:) SecondPoint MinsLoc(SecondPoint,:)]; 
    
for i=2:Cnumber
    Points(i,:)=Points(1,:);
    Points(i,3)=mod(Points(i,3)+(i-1)*360/Cnumber,360);
    Points(i,6)=mod(Points(i,6)+(i-1)*360/Cnumber,360);
    
    Points(i,1)=GetClosestMinimum(Points(i,[2 3]),MinsLoc);
    Points(i,4)=GetClosestMinimum(Points(i,[5 6]),MinsLoc);   
end

[CurrentSyms,~]=size(Points);
try
    if Symmetry(3)=='v'
        for i=CurrentSyms+1:2*CurrentSyms
            Points(i,:)=Points(i-CurrentSyms,:);
            Points(i,3)=360-Points(i,3);
            Points(i,6)=360-Points(i,6);
            
            Points(i,1)=GetClosestMinimum(Points(i,[2 3]),MinsLoc);
            Points(i,4)=GetClosestMinimum(Points(i,[5 6]),MinsLoc);
        end
    end
catch
end

[CurrentSyms,~]=size(Points);
try
    if Symmetry(3)=='h'
        for i=CurrentSyms+1:2*CurrentSyms
            Points(i,:)=Points(i-CurrentSyms,:);
            Points(i,2)=180-Points(i,2);
            Points(i,5)=180-Points(i,5);
            
            Points(i,1)=GetClosestMinimum(Points(i,[2 3]),MinsLoc);
            Points(i,4)=GetClosestMinimum(Points(i,[5 6]),MinsLoc);
        end
    end
catch
end

end