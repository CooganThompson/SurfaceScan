function PlotTSEnergiesOnPES(PESFull)

hold on
for i=1:length(PESFull.Barriers.Energies)
    
    x=PESFull.Barriers.Energies{i}(2);
    x=cellstr(num2str(x'));
    %index=floor(length(PESFull.Barriers.Pathway{i})/2);
    plot(PESFull.Barriers.TSLocaton{i}(1),PESFull.Barriers.TSLocaton{i}(2),'^w','markers',5,'linewidth',5)
    text(PESFull.Barriers.TSLocaton{i}(1),PESFull.Barriers.TSLocaton{i}(2),x,'color','w','FontSize',10,'VerticalAlignment','Bottom','HorizontalAlignment','right')
    %plot(PESFull.Barriers.Pathway{i}(index,1),PESFull.Barriers.Pathway{i}(index,2),'^w','markers',5,'linewidth',5)
    %text(PESFull.Barriers.Pathway{i}(index,1),PESFull.Barriers.Pathway{i}(index,2),x,'color','w','FontSize',10,'VerticalAlignment','Bottom','HorizontalAlignment','right')
    
end
end