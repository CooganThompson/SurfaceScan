function [CorrectionFunctions, XYZCoordsa0b0]=...
    determineMatlabAndOwnCoordsCorrectionFunctions(Config1,Config2)
% // Determines the functions required to transform
% // our spherical coordinates (alpha,beta) to
% // MATLAB's inputs to sph2cart (azimuth,elevation)


if isempty(Config1) && isempty(Config2)
    CorrectionFunctions.CorrectionFunctionAlpha = @(x) 90-x;
    CorrectionFunctions.CorrectionFunctionBeta = @(x) x;
    XYZCoordsa0b0 = [];
    return
end
    
[Config1, Config1RotatedToZAxis] = determineMatlabAlphaBetaFor(Config1);

[Config2, Config2RotatedToZAxis] = determineMatlabAlphaBetaFor(Config2);

XYZCoordsa0b0 = ...
    determineXYZCoordsa0b0(Config1RotatedToZAxis, Config2RotatedToZAxis, ...
    Config1, Config2);

[CorrectionFunctionAlpha,CorrectionFunctionBeta]=...
    determineAlphaAndBetaCorrectionFunctions(Config1,Config2);

CorrectionFunctions.CorrectionFunctionAlpha=CorrectionFunctionAlpha;
CorrectionFunctions.CorrectionFunctionBeta=CorrectionFunctionBeta;

end

function [CorrectionFunctionAlpha,CorrectionFunctionBeta]=...
    determineAlphaAndBetaCorrectionFunctions(Config1,Config2)

Config1=determineMatlabAlphaBetaFor(Config1);
Config2=determineMatlabAlphaBetaFor(Config2);

determineCorrectionFunctionAlpha = @(Config1,Config2) ...
    determineCorrectionFunction(Config1,Config2,'alpha');
determineCorrectionFunctionBeta = @(Config1,Config2) ...
    determineCorrectionFunction(Config1,Config2,'beta');


%CorrectionFunctionAlpha=determineCorrectionFunctionAlpha(Config1,Config2);
%should always be the case as we start from z axis and go down, but matlab
%starts from xy plane and goes up to z axis (matlab calls it elevation)
CorrectionFunctionAlpha = @(x) -mod(x,180)+90;
CorrectionFunctionBeta=determineCorrectionFunctionBeta(Config1,Config2);

end

function CorrectionFunction=determineCorrectionFunction(Config1,Config2,coord)

if strcmpi(coord,'alpha')
    idx_coord=1;
    ModAmount=180;
else
    idx_coord=2;
    ModAmount=360;
end

testIfReflectionAndTranslationAlpha= @(Config1,Config2) ...
    testIfReflectionAndTranslation(Config1,Config2,'alpha');

testIfReflectionAndTranslationBeta= @(Config1,Config2) ...
    testIfReflectionAndTranslation(Config1,Config2,'beta');




if strcmpi(coord,'alpha')
    
    isReflectionAndTranslation=...
        testIfReflectionAndTranslationAlpha(Config1,Config2);
else
    isReflectionAndTranslation=...
        testIfReflectionAndTranslationBeta(Config1,Config2);
end


isCoordsEqual=Config1.MatlabAlphaBeta.isCoordsEqual(Config2.MatlabAlphaBeta);
if ~isReflectionAndTranslation || isCoordsEqual(idx_coord)
    % // Assume just translation
    TranslatedCoords=Config1.MatlabAlphaBeta-Config1.OwnAlphaBeta;
    TranslationAmount=TranslatedCoords.coords(idx_coord);
    CorrectionFunction=@(x) mod(x+TranslationAmount,ModAmount);
else
    % // translation + reflection
    ReflectedCoords=Config1.OwnAlphaBeta.reflectAboutAxis(idx_coord);
    ReflectedAndTranslatedCoords=Config1.MatlabAlphaBeta-ReflectedCoords;
    TranslationAmount=ReflectedAndTranslatedCoords.coords(idx_coord);
    
    CorrectionFunction=@(x) mod(-x+TranslationAmount,ModAmount);
end
end

function isCoordReflectedAndTranslated=...
    testIfReflectionAndTranslation(Config1,Config2,coord)

if strcmpi(coord,'alpha')
    idx_coord=1;
elseif strcmpi(coord,'beta')
    idx_coord=2;
end

% //
getProposedCoordTranslationFrom = @(Config) ...
    Config.MatlabAlphaBeta-Config.OwnAlphaBeta.reflectAboutAxis(idx_coord);
proposedCoordTranslation=...
    getProposedCoordTranslationFrom(Config1);


% //
getTrialCoordsAfterReflectionAndTranslation =@(Config,ProposedCoordTranslation) ...
    Config.OwnAlphaBeta.reflectAboutAxis(idx_coord)+ProposedCoordTranslation;


trial_coords=...
    getTrialCoordsAfterReflectionAndTranslation(Config2,proposedCoordTranslation);


% //
trial_coords=trial_coords.coords;
verification_coords=Config2.MatlabAlphaBeta.coords;

if isnear(trial_coords(idx_coord),verification_coords(idx_coord))
    isCoordReflectedAndTranslated=true;
else
    isCoordReflectedAndTranslated=false;
end

end

function [Config, CoordsRotatedToZAxis] = determineMatlabAlphaBetaFor(Config)
% //
CoordsRotatedToZAxis = rotateCOMScanningAtomAxisToZAxis(Config);
% //
alphaRange=[-90 90];
betaRange=[0 360];

ScanningAtomTestAlphaBetaCoords=...
    makeAlphaBetaTestCombinations(alphaRange,betaRange);

% //
ScanningAtomTestXYZCoords=...
    convertAlphaBetaCombinationsToXYZ(Config,CoordsRotatedToZAxis,...
    ScanningAtomTestAlphaBetaCoords);

% //
CoordsRotatedToZAxis_ScanningAtom = ...
    cartesianCoords(CoordsRotatedToZAxis(Config.ScanningAtomIdx,:));

MatlabAlphaBetaCombination=...
    determineAlphaBetaCombinationWithSmallestDifference(ScanningAtomTestXYZCoords,...
    ScanningAtomTestAlphaBetaCoords,CoordsRotatedToZAxis_ScanningAtom);

% //
Config.MatlabAlphaBeta=sphericalPESCoords(MatlabAlphaBetaCombination);
end

function CoordsRotatedToZAxis=rotateCOMScanningAtomAxisToZAxis(Config)

AtomLinkedToCOMIdx=Config.AtomLinkedToCOMforZAxisAlignmentIdx;
XYZDirectlyFromZMAT=Config.XYZDirectlyFromZMAT;

MetalPartCOM=Config.COM;

CoordsRotatedToZAxis = ...
    alignWithZAxis(XYZDirectlyFromZMAT,AtomLinkedToCOMIdx,MetalPartCOM);

end

function alphaAndBetaTestCoords=...
    makeAlphaBetaTestCombinations(alphaRange,betaRange)

alphaTest=alphaRange(1):0.5:alphaRange(2);
betaTest=betaRange(1):0.5:betaRange(2);

alphaAndBetaCombinations=makeMatrixCombinations(alphaTest(:),betaTest(:));
alphaAndBetaTestCoords=sphericalPESCoords(alphaAndBetaCombinations);
end

function ScanningAtomTestXYZCoords=...
    convertAlphaBetaCombinationsToXYZ(Config,CoordsRotatedToZAxis,...
    ScanningAtomTestCoords)

ScanningAtomIdx=Config.ScanningAtomIdx;
CoordsRotatedToZAxis_ScanningAtom=CoordsRotatedToZAxis(ScanningAtomIdx,:);
ScanningAtomToCOMDist=norm(CoordsRotatedToZAxis_ScanningAtom);


ScanningAtomTestCoords=ScanningAtomTestCoords.setRadius(ScanningAtomToCOMDist);
ScanningAtomTestXYZCoords=ScanningAtomTestCoords.createXYZCoords;
end

function bestAlphaBetaCombination=...
    determineAlphaBetaCombinationWithSmallestDifference(ScanningAtomTestXYZCoords,...
    ScanningAtomTestAlphaBetaCoords,CoordsRotatedToZAxis_ScanningAtom)

% // Using class method
iBestComb=...
    ScanningAtomTestXYZCoords.getClosestElementsWith...
    (CoordsRotatedToZAxis_ScanningAtom);

bestAlphaBetaCombination=ScanningAtomTestAlphaBetaCoords.coords(iBestComb,:);

end
function XYZCoordsa0b0 = ...
    determineXYZCoordsa0b0(Config1RotatedToZAxis, Config2RotatedToZAxis, ...
    Config1, Config2)


if norm(Config1RotatedToZAxis(Config1.MetalAtomsIdx,:) - ...
        Config2RotatedToZAxis(Config2.MetalAtomsIdx,:)) > 0.2
    error('XYZCoordsa0b0 predicted by the two configs is different! Check COM/Atom linked Idx')
else
    XYZCoordsa0b0 = Config1RotatedToZAxis(Config1.MetalAtomsIdx,:);
end


end







