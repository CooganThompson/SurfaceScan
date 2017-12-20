function CleanPlot(NumOfChildrenShouldBeThere)

children = get(gca, 'children');
NumChildThereAre=length(children);
delete(children(1:NumChildThereAre-NumOfChildrenShouldBeThere));

end