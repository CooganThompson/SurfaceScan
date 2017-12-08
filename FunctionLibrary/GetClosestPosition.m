function Positions=GetClosestPosition(Point,MatForm)

DifferencesInPositions=Point-MatForm(:,[1,2]);
MagnitudeOfDifferences=(DifferencesInPositions(:,1).^2+DifferencesInPositions(:,2).^2);
IndexOfMin=min(MagnitudeOfDifferences)==MagnitudeOfDifferences;
Positions=MatForm(IndexOfMin,[3 4 5]);

end