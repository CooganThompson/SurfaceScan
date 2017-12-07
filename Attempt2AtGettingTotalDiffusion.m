
Time=PESFull.Results.Time(end);

DistanceTraveled=0;
for i=2:2:length(PESFull.Results.Names)
    
    timesdiffused=PESFull.Results.EventsIndividual(end,i-1)+PESFull.Results.EventsIndividual(end,i);
    distance=GetDistanceBetweenTwoPoints(PESFull.Mins(PESFull.Results.Locations(1,i),:),PESFull.Mins(PESFull.Results.Locations(2,i),:),PESFull.HpositionData);
    DistanceTraveled=distance*timesdiffused+DistanceTraveled;
    
end

Mobility=DistanceTraveled/Time %ang/s
MobilityInMetersPerSecond=Mobility*10^-10