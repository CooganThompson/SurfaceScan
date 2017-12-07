function PlotBarReactions(KMCresults,options)

if nargin <= 1
    options=0;
end

%Bar of Site Leaves
count=zeros(1,max(max(KMCresults.Locations)));
for i=2:2:length(KMCresults.Names)
    count(KMCresults.Locations(1,i-1))=KMCresults.EventsIndividual(end,i-1)+count(KMCresults.Locations(1,i-1));
end
figure
bar(count)
title('Number of Site Visits')
xlabel('Site Number')
ylabel('Times Visited')
if options=='log'
    set(gca,'YScale','log');
    ylabel('log of Times Visited')
end
%Time on each site
count2=zeros(1,max(max(KMCresults.Locations)));
for i=2:2:length(KMCresults.Names)
    count2(KMCresults.Locations(1,i-1))=KMCresults.EventsIndividual(end,i-1)*KMCresults.ResidenceTimesIndividual(end,i-1)+count2(KMCresults.Locations(1,i-1));
end
figure
bar(count2)
title('Amount of Time on Each Site')
xlabel('Site Number')
ylabel('Time Spent /s')
if options=='log'
    set(gca,'YScale','log');
    ylabel('log of Time Spent /s')
end
end