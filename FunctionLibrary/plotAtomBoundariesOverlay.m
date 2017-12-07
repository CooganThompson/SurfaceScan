function plotAtomBoundariesOverlay(atomSphericalCoordPositions)
% Plots the boundaries corresponding to all points closest to a certain
% atom.
%{
if isempty(hFig)
    hFig=figure;
else
    figure(hFig)
    hold on
end
%}
Colormap=[...
    143,94,26; ... %brown
    247,104,247; ... %pink
    161,255,74;... % green
    147,10,201]/255; %purple

nAtoms = length(atomSphericalCoordPositions);
for i = 1:nAtoms
    
    % Get all coords which are closest to current atom
    CurrAlphaBeta=atomSphericalCoordPositions{i};
    
    CurrColor=chooseColorBasedOnAtomIdx(i,Colormap);
    
    [xt,yt]=computeAlphaShapeAndBoundaries(CurrAlphaBeta);
    
    for j=1:length(xt)
        plot(xt{j},yt{j},'linewidth',2,'Color',CurrColor);
        hold on
    end
end


axis tight
axis equal
set(gca,'xlim',[0 177])
set(gca,'ylim',[0 360])
view(90,90)
set(gcf, 'Position', get(0,'Screensize')) % maximise figure
end
function CurrColor=chooseColorBasedOnAtomIdx(i,Colormap)

CurrColor=[102,255,255]/255;
%{
switch i
    case {8,17}
        CurrColor=Colormap(1,:);
    case {3,4,5,12,18,13}
        CurrColor=Colormap(2,:);
    case {11,6,2,14}
        CurrColor=Colormap(3,:);
    case {1,10,15,9,16,7}
        CurrColor=Colormap(4,:);
    otherwise
        error('hi')
end
%}
end
function [xt,yt]=computeAlphaShapeAndBoundaries(CurrAlphaBeta)

if isempty(CurrAlphaBeta)
    xt{1}=[0 0];
    yt{1}=[0 0];
    warning('Coogan doesnt like that you have an atom with no closest adsorbate')
    return
end

shp =alphaShape(CurrAlphaBeta(:,1),CurrAlphaBeta(:,2));
shp.Alpha = 2.5;
for z=1:shp.numRegions
    [~, xyz] = boundaryFacets(shp,z);
    xyz(end+1,:)=xyz(1,:);
    x=xyz(:,1);
    y=xyz(:,2);
    
    % Remove all points at border
    y(isnear(x,0))=[];
    x(isnear(x,0))=[];
    
    y(isnear(x,177))=[];
    x(isnear(x,177))=[];
    
    x(isnear(y,0))=[];
    y(isnear(y,0))=[];
    
    x(isnear(y,360))=[];
    y(isnear(y,360))=[];
    
    % Resize the polygons by factor A
    A=0.995;
    xt{z} = A*x+(1-A)*mean(x(1:end-1));
    yt{z} = A*y+(1-A)*mean(y(1:end-1));
    
    
end
end