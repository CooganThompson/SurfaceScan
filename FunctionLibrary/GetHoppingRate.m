function HoppingRate=GetHoppingRate(HistoryData,Points1,Points2,Points3)

if nargin==4
    HoppingRate1=GetHoppingRate(HistoryData,Points1,Points2);
    HoppingRate2=GetHoppingRate(HistoryData,Points1,Points3);
    HoppingRate3=GetHoppingRate(HistoryData,Points2,Points3);
    
    HoppingRate=(1/HoppingRate1+1/HoppingRate2+1/HoppingRate3)^-1;
    
elseif nargin==3
    
    temp1=sum(HistoryData.PositionData(:,Points1,3),2);
    temp2=sum(HistoryData.PositionData(:,Points2,3),2);
    
    count=0;
    
    STATE=3;
    NEWSTATE=3;
    if temp1(1)==1
        STATE=1;
    elseif temp2(2)==1
        STATE=2;
    else
        count=-0.5;
    end
    
    
    for i=1:length(temp1)
        
        if temp1(i)==1
            NEWSTATE=1;
        elseif temp2(i)==1
            NEWSTATE=2;
        end
        
        if STATE~=NEWSTATE
            count=count+1;
            STATE=NEWSTATE;
        end
    end
    count
    HoppingRate=count/HistoryData.Time(end);
    
end

end