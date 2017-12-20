function Point=AveragePointsAroundACylinder(Points,tolerance)
[NumberOfPoints,~]=size(Points);
if NumberOfPoints==1
    Point=Points;
else
    
    if sum(Points(:,2)<tolerance)>0 && sum(Points(:,2)>360-tolerance)>0
        for k=1:NumberOfPoints
            if Points(k,2)>360-tolerance
                Points(k,2)=Points(k,2)-360;
            end
        end
    end
    Point=mean(Points,1);
    
    if Point(2)<0
        Point(2)=Point(2)+360;
    end    
    
end

end