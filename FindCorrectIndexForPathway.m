function count=FindCorrectIndexForPathway(Barriers,FirstPoint,SecondPoint)

    foundStartSpot=0;
    foundEndSpot=0;
    count=0;
    

    while foundStartSpot==0
        count=count+1;
        
        if count>length(Barriers.StartMinIndex)
            foundStartSpot=1;
            foundEndSpot=1;
            count=count+1;
        elseif FirstPoint==Barriers.StartMinIndex{count}
            foundStartSpot=1;
        elseif FirstPoint<Barriers.StartMinIndex{count}
            foundStartSpot=1;
            foundEndSpot=1;
            count=count+1;

        end
    end
    
    count=count-1;
    
    while foundEndSpot==0
        count=count+1;
        
        if count==length(Barriers.EndMinIndex)
            foundEndSpot=1;
            count=count+1;
        elseif SecondPoint<=Barriers.EndMinIndex{count}||FirstPoint<Barriers.StartMinIndex{count}
            foundEndSpot=1;
        end
        
    end
        
end