function PlotMapReactions(PESFull)

KMCresults=PESFull.Results;
PESFull.Class.plotContourMapsView
PlotMinsOnPES(PESFull);

MaxSteps=max(KMCresults.EventsIndividual(end,:));
hold on
for i=2:2:length(KMCresults.Names)
%     Mins=KMCresults.Locations(:,i);
%     Min1=PESFull.Mins(Mins(1),:);
%     Min2=PESFull.Mins(Mins(2),:);
    
    Events=KMCresults.EventsIndividual(end,i-1);
    if Events>0
        Coords=cell2mat(PESFull.Barriers.Pathway(i/2));
        plotbrokenlines(Coords(:,1),Coords(:,2),'r','LineWidth',10*Events/MaxSteps)
        
        
        
%         q=quiver( Min1(1),Min1(2),Min2(1)-Min1(1),Min2(2)-Min1(2),0);
%         q.Color='red';
%         q.LineWidth=10*Events/MaxSteps;
    end
%     
%     Min2=PESFull.Mins(Mins(1),:);
%     Min1=PESFull.Mins(Mins(2),:);
    
    Events=KMCresults.EventsIndividual(end,i);
    if Events>0    
        Coords=cell2mat(PESFull.Barriers.Pathway(i/2));
        plotbrokenlines(Coords(:,1),Coords(:,2),'r','LineWidth',10*Events/MaxSteps)
        
        
%         q=quiver( Min1(1),Min1(2),Min2(1)-Min1(1),Min2(2)-Min1(2),0);
%         q.Color='red';
%         q.LineWidth=10*Events/MaxSteps;
    end
end
