function MinsLoc=GetMins(PES)

tol=1e-6;
%Must put this into matrix form
PES=PES.interpolatePESGrid;
%get mins and puts them into better coordinates
BW=imregionalmin(PES.energiesGrid);
scale=(size(PES.energiesGrid)-1)./[360 180];
if length(unique(scale)) ~= 1
    error('Something happened in the scaling, x and y scales are not equal')
end
[r,c]=find(BW);
R=1/scale(1)*(r-1);
C=1/scale(1)*(c-1);
i=1;
while i < length(R)
    X=sind(C).*cosd(R);
    Y=sind(C).*sind(R);
    Z=cosd(C);
    
    dist=(X(i)-X(i+1:end)).^2+(Y(i)-Y(i+1:end)).^2+(Z(i)-Z(i+1:end)).^2;
        
    distmin=.007;
    if sum(dist==0)>0
        for j=find(~dist)
            R(i+j)=[];
            C(i+j)=[];
        end
    elseif sum(dist<=distmin)>0
        
        Energies=nan(length(R),1);
        Energies(i)=PES.energiesGrid(R(i)*scale(1)+1,C(i)*scale(2)+1);
        for j=find(dist<=distmin)'
            Energies(i+j)=PES.energiesGrid(R(i+j)*scale(1)+1,C(i+j)*scale(2)+1);
        end
        
        tempR=R(Energies-min(Energies)<1e-4);
        tempC=C(Energies-min(Energies)<1e-4);
        
        tolerance=10;
        Point=AveragePointsAroundACylinder([tempC tempR],tolerance);

%         if length(tempR)>1
%             if sum(tempR<10)>0 && sum(tempR>350)>0
%                 for k=1:length(tempR)
%                     if tempR(k)>350
%                         tempR(k)=tempR(k)-360;
%                     end
%                 end
%             end
%             R(i)=mean(tempR);
%             C(i)=mean(tempC);
%             if R(i)<0
%                 R(i)=R(i)+360;
%             end
% %         else
%             
%         end
        
        R(i)=Point(2);
        C(i)=Point(1);
        
        C(i)=round(scale(1)*Point(1))/scale(1);
        R(i)=round(scale(1)*Point(2))/scale(1);
        
        
%         [C,R]
%         if C(2)==58.7500
%             3
%         end
        
        
        for j=find(dist<=distmin)
            R(i+j)=[];
            C(i+j)=[];
        end
        
    else
        i=i+1;
    end
end

MinsLoc=[C,R];

end