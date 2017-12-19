function Points=Cleanpoints(Points)

[NumPoints,~]=size(Points);

for i=1:NumPoints
    if Points(i,1)>Points(i,4)
        Points(i,:)=circshift(Points(i,:),3);
    end
end

Points(Points(:,3)==360,3)=0;   %Cause 0=360
Points(Points(:,6)==360,6)=0;

Points(Points(:,2)<=1,3)=Points(1,3);   %To make sure all the top pieces are treated as the same point. Might need adjusting if the top isn't top enough
Points(Points(:,5)>=179,6)=Points(1,6); %Similar thing for the bottom

flag=0;
counter=1;
while counter<=NumPoints
    removeindex=find(sum((Points(counter,:)-Points([counter+1:NumPoints],:)).^2,2)==0); %This might also benefit from adding some precision
    %removeindex=find(sum(Points(counter,:)-Points([counter+1:NumPoints],:),2)==0);
    if isempty(removeindex)
        counter=counter+1;
    else
        Points(removeindex+counter,:)=[];
    end
    [NumPoints,~]=size(Points);
end
end