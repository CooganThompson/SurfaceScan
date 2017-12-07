function [atomSphericalCoordPositions, CenteredSubstrateCoords, GRID_SIZE] = ...
    getSphericalCoordsRegionsForAtoms(CorrectionFunctions, XYZCoordsa0b0, GRID_SIZE)

CenteredSubstrateCoords = centerSubstrateCoordinates(XYZCoordsa0b0);

if ~exist('GRID_SIZE','var')
GRID_SIZE = 0.2;
end
beta = 0:GRID_SIZE:360;
alpha = 0:GRID_SIZE:180;
alphaAndBetaCombinations=makeMatrixCombinations(alpha(:),beta(:));

SphereCoords=determineSphereCoordinatesFromAlphaBeta(alphaAndBetaCombinations,...
    CorrectionFunctions);

% // If something is wrong, run these commands to debug.
% This should give the cluster, in the a0_b0 position, in the
% center of a sphere of dots
%{
figure
scatter3(SphereCoords.coords(1:1:end,1),SphereCoords.coords(1:1:end,2),...
    SphereCoords.coords(1:1:end,3))
hold on
scatter3(CenteredSubstrateCoords(:,1),CenteredSubstrateCoords(:,2),...
    CenteredSubstrateCoords(:,3),'filled')
%}


IdxNearestAtoms=...
    calculateNearestAtomsToPointsOnSphere(CenteredSubstrateCoords,SphereCoords);


nAtoms=size(XYZCoordsa0b0,1);
atomSphericalCoordPositions = cell(1, nAtoms);
for i = 1:nAtoms
    idxPositions = IdxNearestAtoms == i;
    atomSphericalCoordPositions{i} =  alphaAndBetaCombinations(idxPositions,:);
end

end



function CenteredSubstrateCoords = centerSubstrateCoordinates(XYZCoordsa0b0)
COM=mean(XYZCoordsa0b0,1);
CenteredSubstrateCoords=bsxfun(@minus,XYZCoordsa0b0,COM);
end

function CartesianCoords=...
    determineSphereCoordinatesFromAlphaBeta(alphaAndBetaCombinations,...
    CorrectionFunctions)



radius=6;
alphas = alphaAndBetaCombinations(:,1);
betas = alphaAndBetaCombinations(:,2);
if ~isempty(CorrectionFunctions)
    alphas = CorrectionFunctions.CorrectionFunctionAlpha(alphas);
    betas = CorrectionFunctions.CorrectionFunctionBeta(betas);
end

AlphaBetaCoords=sphericalPESCoords([alphas,betas]);
AlphaBetaCoords=AlphaBetaCoords.setRadius(radius);
CartesianCoords=AlphaBetaCoords.createXYZCoords;

end


function IdxNearestAtoms=...
    calculateNearestAtomsToPointsOnSphere(CenteredAuCoords,SphereCoords)
% // Calculates distances of cluster atoms to points on a sphere around it



% Compute the shortest distance, and determine the atom closest to the points
% of the sphere
dists=euclideanDist(CenteredAuCoords,SphereCoords.coords);
[~,IdxNearestAtoms]=min(dists);

end


