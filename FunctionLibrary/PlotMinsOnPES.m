function PlotMinsOnPES(PESFull)

MinsLoc=PESFull.Mins;
hold on
x=1:length(MinsLoc(:,1));
x=cellstr(num2str(x'));
plot(MinsLoc(:,1),MinsLoc(:,2),'ow','markers',5,'linewidth',5)
text(MinsLoc(:,1),MinsLoc(:,2),x,'color','w','FontSize',20,'VerticalAlignment','Bottom','HorizontalAlignment','right')

end