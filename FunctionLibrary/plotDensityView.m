function plotDensityView
COM = mean(XYZCoordsa0b0,1);
centeredSubstrateCoords = bsxfun(@minus,XYZCoordsa0b0,COM);
atomSphereSize = 2.5;
PESSphereSize = 3;
nSpherePoints = 300;
gridsize = 3;
faceAlpha = 0.7;

% // Self is an object that contains properties
%        alpha
%        beta
%        energies
%        alphaGrid
%        betaGrid
%        energiesGrid

% // Creates alphaGrid/betaGrid/energiesGrid from alpha & beta
self = self.interpolatePESGrid(gridsize);

% // Determines the functions required to transform
% // our spherical coordinates (alpha,beta) to 
% // MATLAB's inputs to sph2cart (azimuth,elevation)
[CorrectionFunctions]=...
    determineMatlabAndOwnCoordsCorrectionFunctions(XYZDirectlyFromZMatConfig1,...
    XYZDirectlyFromZMatConfig2);

% // Converts (alpha,beta) to cartesian coords of a unit sphere
% (first converting to azimuth and elevation)

[sphereCoords, sphereCoordsEnergies] = ...
    getPESUnitSphereRepresentation(CorrectionFunctions);

% // Gets spheres for substrate atoms
atomShapes = ...
    getSubstrateAtomShapes(centeredSubstrateCoords, atomSphereSize, nSpherePoints);

% // Gets "skin" of PES
PESshape = ...
    getPESShape(centeredSubstrateCoords, PESSphereSize, nSpherePoints);

% // Determine the energies on each point of the PES "skin" by finding closest distance
% // from skin to unit sphere.
PESshapeEnergies = ...
    getPESShapeEnergiesAtPoints(PESshape, sphereCoords, sphereCoordsEnergies);

% // Plot PES skin
[figHandle, shapeHandle] = ...
    plotPESshape(PESshape, PESshapeEnergies, faceAlpha);

% // Plot substrate atoms
plotSubstrateAtomShapes(atomShapes, figHandle)

xlim([-10,10])
ylim([-10,10])
zlim([-10,10])
end
%------------------------------------------------------------
%
function [sphereCoords, sphereCoordsEnergies] = ...
    getPESUnitSphereRepresentation(self, CorrectionFunctions)

correctedAlpha = ...
    CorrectionFunctions.CorrectionFunctionAlpha(self.alphaGrid(:));
correctedBeta = ...
    CorrectionFunctions.CorrectionFunctionBeta(self.betaGrid(:));


radius=1;
AlphaBetaCoords=sphericalPESCoords([correctedAlpha,correctedBeta]);
AlphaBetaCoords=AlphaBetaCoords.setRadius(radius);
CartesianSphereCoords=AlphaBetaCoords.createXYZCoords;

sphereCoords = CartesianSphereCoords.coords;
sphereCoordsEnergies= self.energiesGrid(:);
end
%------------------------------------------------------------
%
function atomShapes = ...
    getSubstrateAtomShapes(CenteredSubstrateCoords, atomSphereSize, nSpherePoints)

alphaShapeAlpha = 4;

[smallsphereX,smallsphereY,smallsphereZ] = sphere(nSpherePoints);
smallsphereX = smallsphereX(:);
smallsphereY = smallsphereY(:);
smallsphereZ = smallsphereZ(:);

nAtoms = size(CenteredSubstrateCoords,1);
x = cell(1,nAtoms);
y = cell(1,nAtoms);
z = cell(1,nAtoms);

atomShapes = cell(1,nAtoms);
for i = 1:nAtoms
    AtomCoords = CenteredSubstrateCoords(i,:);
    
    x{i} = smallsphereX*atomSphereSize + AtomCoords(1);
    y{i} = smallsphereY*atomSphereSize + AtomCoords(2);
    z{i} = smallsphereZ*atomSphereSize + AtomCoords(3);
    
    atomShapes{i} = alphaShape(x{i},y{i},z{i},alphaShapeAlpha,...
        'HoleThreshold',4);
    
end
end

%------------------------------------------------------------
%
function PESshape = ...
    getPESShape(CenteredSubstrateCoords, PESSphereSize, nSpherePoints)

[smallsphereX,smallsphereY,smallsphereZ] = sphere(nSpherePoints);
smallsphereX = smallsphereX(:);
smallsphereY = smallsphereY(:);
smallsphereZ = smallsphereZ(:);

nAtoms = size(CenteredSubstrateCoords,1);
x = cell(1,nAtoms);
y = cell(1,nAtoms);
z = cell(1,nAtoms);

for i = 1:nAtoms
    AtomCoords = CenteredSubstrateCoords(i,:);
    
    x{i} = smallsphereX*PESSphereSize + AtomCoords(1);
    y{i} = smallsphereY*PESSphereSize + AtomCoords(2);
    z{i} = smallsphereZ*PESSphereSize + AtomCoords(3);
end

PESshape = alphaShape(vertcat(x{:}),vertcat(y{:}),vertcat(z{:}));
PESshape.Alpha = 2;
PESshape.HoleThreshold = 50;

end
%------------------------------------------------------------
%
function idxNearestAtoms = ...
    getNearestPointsFromPESShapeToPESSphereRep(PESshape, sphereCoords)

nPESshapePoints = size(PESshape.Points,1);
windowSize = 1000;
nWindows = ceil(nPESshapePoints/windowSize);
idxNearestAtoms = zeros(nPESshapePoints,1);

for i = 1:nWindows
    
    windowStart = (i-1)*windowSize+1;
    windowEnd = min([i*windowSize, nPESshapePoints]);
    currWindow = windowStart:windowEnd;
    
    % Single to save memory
    distsTmp = euclideanDist(single(sphereCoords),...
        single(PESshape.Points(currWindow,:)));
    
    [~,idxNearestAtoms(currWindow)] = min(distsTmp);
    
end
end

%------------------------------------------------------------
%
function [figHandle, shapeHandle] = ...
    plotPESshape(PESshape, PESshapeEnergies, faceAlpha)

[triangulation, vertexCoords] = PESshape.boundaryFacets();

% As tri refers to rows of vertexCoords, and not PESShape,
% we need to get the indices of PESShape.Points
% then we can get verticesEnergy
[~,PESshapeIndices] = ismember(vertexCoords,PESshape.Points,'rows');
verticesEnergy = PESshapeEnergies(PESshapeIndices);

figHandle = figure;
light('Position',[1 1 1])
shapeHandle = patch('Faces',triangulation, 'Vertices', vertexCoords, ...
    'FaceVertexCData',...
    verticesEnergy,'FaceColor','interp','linestyle','none',...
    'FaceAlpha',faceAlpha);

shapeHandle.DiffuseStrength = 1;
shapeHandle.SpecularStrength = 0;
shapeHandle.AmbientStrength = 1;

end

%------------------------------------------------------------
%
function plotSubstrateAtomShapes(atomShapes, figHandle)

figure(figHandle)
hold on

for i = 1:length(atomShapes)
    a=atomShapes{i}.plot('linestyle','none','FaceAlpha',0.9,...
        'FaceColor','r') ;
    a.SpecularStrength = 1;
    a.AmbientStrength = 1;
    hold on
end

end

%------------------------------------------------------------
%
function shpEnergies = ...
    getPESShapeEnergiesAtPoints(PESShape, sphereCoords, ...
    sphereCoordsEnergies)

IdxNearestAtoms = ...
    getNearestPointsFromPESShapeToPESSphereRep(PESShape, sphereCoords);

shpEnergies = sphereCoordsEnergies(IdxNearestAtoms);

end
