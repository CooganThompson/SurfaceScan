function MinimumNumber=GetClosestMinimum(Point,MinsLoc)

r=ones(length(MinsLoc),1);
[x1,y1,z1] = sph2cart(MinsLoc(:,2)*pi/180,(90-MinsLoc(:,1))*pi/180,r);
[x2,y2,z2] = sph2cart(Point(:,2)*pi/180,(90-Point(:,1))*pi/180,1);

%  [x1,y1,z1] = sph2cart(MinsLoc(:,1),MinsLoc(:,2),r);
%  [x2,y2,z2] = sph2cart(Point(:,1),Point(:,1),2);

diff=[x1 y1 z1]-[x2 y2 z2];

% 
% diff=MinsLoc-Point;
% diffoffset1=MinsLoc-(Point-[0,360]); %checks just in case the point is on the wrong side of the wrap around (i.e. so 0 and 360 are counted equally
% diffoffset2=MinsLoc-(Point+[0,360]);



diffsum=sum(diff.^2,2);
%diffoffset1sum=sum(diffoffset1.^2,2);
%diffoffset2sum=sum(diffoffset2.^2,2);

%method I did before
%temp=find([diffsum; diffoffset1sum; diffoffset2sum]==min([diffsum; diffoffset1sum; diffoffset2sum]));
temp=find(diffsum==min(diffsum));
MinimumNumber=mod(temp,length(MinsLoc));

if MinimumNumber==0 %Because the modulo function isn't exactly what I want
    MinimumNumber=length(MinsLoc);
end



end