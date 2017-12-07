function distance=GetDistanceBetweenTwoPoints(Point1,Point2,PositionData)

MatForm=cell2mat(PositionData');
Point1Coords=GetClosestPosition(Point1,MatForm);
Point2Coords=GetClosestPosition(Point2,MatForm);
vectorfrom1to2=mean(Point1Coords,1)-mean(Point2Coords,1);
distance=norm(vectorfrom1to2);

end

function Positions=GetClosestPosition(Point,MatForm)

DifferencesInPositions=Point-MatForm(:,[1,2]);
MagnitudeOfDifferences=(DifferencesInPositions(:,1).^2+DifferencesInPositions(:,2).^2);
IndexOfMin=min(MagnitudeOfDifferences)==MagnitudeOfDifferences;
Positions=MatForm(IndexOfMin,[3 4 5]);

end