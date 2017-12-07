function SpiderWebPlot(PESFull)

hold on
for i=1:length(PESFull.Barriers.Name)
    %PESFull.Barriers.Pathway{i};
    
    plotbrokenlines(PESFull.Barriers.Pathway{i}(:,1),PESFull.Barriers.Pathway{i}(:,2),'w','LineWidth',2)
    %plot(PESFull.Barriers.Pathway{i}(:,1),PESFull.Barriers.Pathway{i}(:,2),'w','LineWidth',4)
end

end