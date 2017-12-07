function [AlphaFinestGrid,BetaFinestGrid,InterpolEnergies]=...
    interpolateGrid(alpha,beta,NormalizedEnergies,gridsize)
% Makes a finer grid, and interpolates the energies. The inputs are all
% vectors, while the outputs will be in matrices




% // Make into a form where alpha,beta, and energies are matrices
% This is necessary to use interp2.
% Grid data does some interpolation too.


alphaVec=unique(alpha);
betaVec=unique(beta);

[alphaGrid,betaGrid]=meshgrid(alphaVec,betaVec);

% Only interpolation, gives NaNs elsewhere
GridDataEnergies = ...
    griddata(alpha,beta,NormalizedEnergies,alphaGrid,betaGrid,'linear');

% For extrapolation
GridDataEnergiesNear = ...
    griddata(alpha,beta,NormalizedEnergies,alphaGrid,betaGrid,'nearest');
posNans = isnan(GridDataEnergies);

GridDataEnergies(posNans) = GridDataEnergiesNear(posNans);


% // Interpolate energies with finer mesh
% Do interpolation
minAlpha=min(alpha);
minBeta=min(beta);
maxAlpha=max(alpha);
maxBeta=max(beta);

AlphaFinestVec=minAlpha:gridsize/2:maxAlpha;
BetaFinestVec=minBeta:gridsize/2:maxBeta;
[AlphaFinestGrid,BetaFinestGrid] = meshgrid(AlphaFinestVec,BetaFinestVec);

InterpolEnergies = ...
    interp2(alphaGrid,betaGrid,GridDataEnergies,...
    AlphaFinestGrid,BetaFinestGrid,'spline');




end