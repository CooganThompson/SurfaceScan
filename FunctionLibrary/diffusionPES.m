classdef diffusionPES
    
    properties (GetAccess=public, SetAccess=public)
        alpha
        beta
        energies
        filename
        stepsizeAlpha
        stepsizeBeta
        lowestenergy
        alphaGrid
        betaGrid
        energiesGrid
        
    end
    
    properties (GetAccess=private,SetAccess=private)
        
    end
    
    methods
        
        %------------------------------------------------------------
        %
        function self = diffusionPES(varargin)
            
            i = 1;
            AutoFillData = true;
            NormalizeEnergies = true;
            while i <= length(varargin)
                if i == 1 && ischar(varargin{i})
                    summaryFileName = varargin{i};
                    inputAlphaBetaPairs = [];
                    inputEnergies = [];
                    i = i+1;
                elseif i == 1 && ~ischar(varargin{i})
                    inputAlphaBetaPairs = varargin{i};
                    inputEnergies = varargin{i+1};
                    summaryFileName = [];
                    i = i+2;
                elseif strcmpi(varargin{i},'noautofilldata')
                    AutoFillData = false;
                    i = i+1;
                    
                elseif strcmpi(varargin{i},'nonormalizeenergies')
                    NormalizeEnergies = false;
                    i = i+1;
                else
                    error(['Usage: diffusionPES(alphaBetaPair, energies, options)' ...
                        '\n or diffusionPES(summaryFileName, options)'])
                end
            end

            DataFolder='data';
            if ~isempty(summaryFileName)
                [alpha,beta,energies] = ...
                    self.extractAlphaBetaAndEnergiesFromFile(summaryFileName,DataFolder);
            else
                alpha = inputAlphaBetaPairs(:,1);
                beta = inputAlphaBetaPairs(:,2);
                energies = inputEnergies;
                
                inputDataLengths = [length(beta) length(alpha) length(energies)];
                if ~all(inputDataLengths == inputDataLengths(1))
                    error('Inputs have different lengths!')
                end

            end
            %{
            [alpha,beta,energies] = ...
                self.removeZeroEnergyEntries(alpha,beta,energies);
            %}
            if NormalizeEnergies
                energies = ...
                    self.normalizeEnergiesWithRespectToMinimum(energies);
                
                energies = ...
                    self.convertEnergiesToEVs(summaryFileName,energies);
                
            end
            [stepSizeAlpha,stepSizeBeta] = self.calcStepSize(alpha,beta);
            if AutoFillData
                [alpha,beta,energies] = ...
                    self.fillInAlpha0Entries(alpha,beta,energies,stepSizeBeta);
                
                [alpha,beta,energies] = ...
                    self.fillInAlpha180Entries(alpha,beta,energies,stepSizeBeta);
                %}
                
                [alpha,beta,energies] = ...
                    self.fillInBeta360Entries(alpha,beta,energies);
            end
            
            self.alpha = alpha;
            self.beta = beta;
            self.lowestenergy = min(energies);
            self.energies = energies;
            self.filename = summaryFileName;
            self.stepsizeAlpha = stepSizeAlpha;
            self.stepsizeBeta = stepSizeBeta;
        end
        
        %------------------------------------------------------------
        %                %
        function self=fillinviasymmetry(self,symmetry) %Uses the point group symmetry operations (kinda) so you can graph the whole surface with only each unique part specified
        disp('Be sure you defined alpha and beta well')
        disp('Currently only good for C_n and C_nv symmetries, C_nh is untested')
            
            if strcmpi(symmetry,'Cs') 
                symmetry='C1v';
            end

            if symmetry(1)=='C'||symmetry(1)=='c'
                n=str2double(symmetry(2));
                for i=2:n
                selfprime{i}=translateByBeta(self,360/n);
                self=self+selfprime{i};
                end
            end
            
            if symmetry(3)=='v'
                selfprime=reflectAboutBetaZero(self);
                self=self+selfprime;
            end
                
            if symmetry(3)=='h'
                selfprime=reflectAboutAlphaZero(self);
                self=self+selfprime;
            end       
            
%             if (symmetry(1)=='T'||symmetry(1)=='t')&&symmetry(2)=='d'
%                 selfprime1=translateByBeta(self,120);
%                 selfprime2=translateByBeta(self,240);
%                 self=self+selfprime1+selfprime2;
%                 
%                 selfprime=reflectAboutBetaZero(self);
%                 self=self+selfprime;
%                 
%                 
%                 
%             end
        end
 
        function self=combineAndGetMinimum(self,differentPES)
            AllAlpha=[self.alpha; differentPES.alpha];
            AllBeta=[self.beta; differentPES.beta];
            AllEnergies=[self.energies+self.lowestenergy; differentPES.energies+differentPES.lowestenergy];
            
            Coords=[AllAlpha AllBeta];
            [UniqueCoords,iUniqueCoords,iUniqueCoordsOri]=unique(Coords,'rows');
            
            % Get lowest energy out of all unique coords
            UniqueEnergy=zeros(size(UniqueCoords,1),1);
            for i=1:size(UniqueCoords,1)
                IdxOldCoords=iUniqueCoordsOri==i;
                UniqueEnergy(i)=min(AllEnergies(IdxOldCoords));
            end
            
            self.alpha=UniqueCoords(:,1);
            self.beta=UniqueCoords(:,2);
            self.energies=UniqueEnergy-min([self.lowestenergy,differentPES.lowestenergy]);
            self.lowestenergy=min([self.lowestenergy,differentPES.lowestenergy]);
        end

        %------------------------------------------------------------
        %
        function Coords=getMinimumEnergyCoords(self)
            [~,iMinEnergy]=min(self.energies);
            Coords(1)=self.alpha(iMinEnergy);
            Coords(2)=self.beta(iMinEnergy);
        end

        %------------------------------------------------------------
        %
        function plotGoogleMapsView(self)
            figure
            triangulations = delaunay(self.alpha,self.beta);
            trisurf(triangulations,self.alpha,self.beta,...
                self.energies)
            xlabel('Alpha');
            ylabel('Beta');
            zlabel('Energy relative to most stable state / eV')
        end
        %------------------------------------------------------------
        %
        function [hContour, hColorbar] = plotDifferenceContourMap(self,differentPES)
            differencePES = self.takeDifference(differentPES);
            [hContour, hColorbar] = differencePES.plotContourMapsView('redbluecolormap');
            
        end
        
        %------------------------------------------------------------
        %
        function [hContour, hColorbar] = plotDifference3DRepresentation(self, ...
                differentPES, XYZDirectlyFromZMatConfig1, XYZDirectlyFromZMatConfig2)
            differencePES = self.takeDifference(differentPES);
            [hContour, hColorbar] = differencePES.plot3DRepresentation(XYZDirectlyFromZMatConfig1, ...
                XYZDirectlyFromZMatConfig2);
            
        end
        
        %------------------------------------------------------------
        %
        function differencePES = takeDifference(self,differentPES)
            % // DifferentPES is taken as final, self is taken as initial
            
            % Update Grid
            self = self.interpolatePESGrid;
            differentPES = differentPES.interpolatePESGrid;
            
            % Get differences
            deltaEnergies = differentPES.energiesGrid - self.energiesGrid;
            
            %
            alphaBetaPairs = [self.alphaGrid(:),self.betaGrid(:)];
            deltaEnergies = deltaEnergies(:);
            differencePES = diffusionPES(alphaBetaPairs, deltaEnergies, ...
                'noautofilldata', 'nonormalizeenergies');
            
        end
        %------------------------------------------------------------
        %
        function [hFig, hColorbar] = plotContourMapsView(self, varargin)
            
            % Options
            i = 1;
            noLabels = false;
            nContours = 60;
            redBlueColorMap = false;
            while i <= length(varargin)
                if strcmpi(varargin{i}, 'ncontours')
                    nContours = varargin{i+1};
                    i = i + 2;
                elseif strcmpi(varargin{i}, 'nolabels')
                    noLabels = true;
                    i = i + 1;
                elseif strcmpi(varargin{i}, 'redbluecolormap')
                    redBlueColorMap = true;
                    i = i + 1;
                elseif isempty(varargin{i})
                    i = i + 1;
                else
                    error('Flag %s unknown!',varargin{i+1});
                end
            end
            % Update Grid
            hFig = figure;
            self=self.interpolatePESGrid;
            
            contourf(self.alphaGrid,self.betaGrid,self.energiesGrid,nContours);
            
            hold on
            view(90,90) %rotate to correct viewing angle
            set(gcf, 'Position', get(0,'Screensize')) % maximise figure
            set(gca,'fontsize',14)
            
            xlabel('Alpha');
            ylabel('Beta');
            
            if redBlueColorMap
                climits = caxis;
                climits = max(abs(climits))*[-1 1];
                colormap(b2r(climits(1),climits(2)));
            end
            
            % // colorbar
            colorbarLabel = 'Energy Relative to Most Stable State / eV';
            hColorbar = self.plotColorbar(hFig, colorbarLabel);

            
            if noLabels
                xlabel('')
                ylabel('')
                set(gca,'xticklabels','')
                set(gca,'yticklabels','')
            end
            
        end
        
        
        %----------------------------------------------------------------
        %
        function [hFig, hColorbar] = plot3DRepresentation(self, XYZDirectlyFromZMatConfig1,...
                XYZDirectlyFromZMatConfig2, CorrectionFunctions, ...
                centeredXYZCoordsa0b0)
            
            ATOM_RADIUS = 1.7;
            N_SPHERE_POINTS = 150;
            hFig = figure;

            self = self.interpolatePESGrid();
            
            if nargin == 3
                [CorrectionFunctions, centeredXYZCoordsa0b0]=...
                    determineMatlabAndOwnCoordsCorrectionFunctions(XYZDirectlyFromZMatConfig1,...
                    XYZDirectlyFromZMatConfig2);
            end
            
            
            alphaGridCorrected = CorrectionFunctions.CorrectionFunctionAlpha(self.alphaGrid);
            betaGridCorrected = CorrectionFunctions.CorrectionFunctionBeta(self.betaGrid);
            
            [X, Y, Z] = sphere(N_SPHERE_POINTS);

            
            for i=1:size(centeredXYZCoordsa0b0,1)
                idxAtom = i;
                
                [x, y, z]= ...
                    self.centerSphereOnAtomLocation(X, Y, Z, centeredXYZCoordsa0b0, ...
                    idxAtom, ATOM_RADIUS);

                
                [currAlpha, currBeta] = ...
                    self.convertSphereCartesianCoordsToAlphaAndBeta(x,y,z);
                
                
                currEnergies =  self.getEnergiesOnSphere(alphaGridCorrected, ...
                    betaGridCorrected, self.energiesGrid, currAlpha, currBeta);
                s = surf(x,y,z,currEnergies,'LineStyle','none');
                s.SpecularStrength = 0.1;
                s.AmbientStrength = 0.4;
                s.DiffuseStrength = 0.8;
                hold on
            end
            
            light('position',[0 0 5]);
            lighting flat;
            %lighting gouraud;
            hold off;
            %   colormap parula
            axis off;
            view([45,30]);
            axis vis3d;
            rotate3d;
            camlight(-180,-90);
            axis tight
            axis equal
            
            %
            colorbarLabel = 'Energy Relative to Most Stable State / eV';
            hColorbar = self.plotColorbar(hFig, colorbarLabel);
        end
        %------------------------------------------------------------
        %
        function hFig = plotExoplanetView(self, varargin)
            % First does contour in 2D, and then transforms levels to
            % spherical form
            
            hFig = figure;
            
            % // update Grid
            self=interpolatePESGrid(self);
            
            % // Contour options
            n2DContourLevels=180;
            nSphereContourLevels=60; % less than n2DContourLevels to increase resolution
            
            
            RatioOf2DToSphereContourLevels=n2DContourLevels/nSphereContourLevels;
            if ~isint(RatioOf2DToSphereContourLevels)
                error('The ratio is not an integer!')
            end
            
            contourMatrixOfLines=...
                self.calculateExoplanetContourLines(n2DContourLevels);
            
            ColorMapSphere=...
                self.createExoplanetColorMap(n2DContourLevels,...
                RatioOf2DToSphereContourLevels);
            
            self.transformAndPlotExoplanetContourLines(contourMatrixOfLines,...
                ColorMapSphere,RatioOf2DToSphereContourLevels);
            
        end
        %------------------------------------------------------------
        %
        function plotDensityView(self,XYZDirectlyFromZMatConfig1,...
                XYZDirectlyFromZMatConfig2,CorrectionFunctions, centeredXYZCoordsa0b0)
            
            
            atomSphereSize = 2.5;
            PESSphereSize = 2.5;
            nSpherePoints = 100;
            gridsize = 3;
            faceAlpha = 0.7;
            
            % // update Grid
            self = self.interpolatePESGrid(gridsize);
            
            if nargin == 3
                [CorrectionFunctions, centeredXYZCoordsa0b0]=...
                    determineMatlabAndOwnCoordsCorrectionFunctions(XYZDirectlyFromZMatConfig1,...
                    XYZDirectlyFromZMatConfig2);
            end
            
            [sphereCoords, sphereCoordsEnergies] = ...
                self.getPESUnitSphereRepresentation(CorrectionFunctions);
            
            
            atomShapes = ...
                self.getSubstrateAtomShapes(centeredXYZCoordsa0b0, atomSphereSize, nSpherePoints);
            
            PESshape = ...
                self.getPESShape(centeredXYZCoordsa0b0, PESSphereSize, nSpherePoints);
            
            PESshapeEnergies = ...
                self.getPESShapeEnergiesAtPoints(PESshape, sphereCoords, sphereCoordsEnergies);
            
            [figHandle, shapeHandle] = ...
                self.plotPESshape(PESshape, PESshapeEnergies, faceAlpha);
            
            self.plotSubstrateAtomShapes(atomShapes, figHandle)
            
            xlim([-10,10])
            ylim([-10,10])
            zlim([-10,10])
            axis vis3d;
            view([45,30]);
            
        end
        
        %------------------------------------------------------------
        %
        function plotAtomRegionsOverlay(self,XYZDirectlyFromZMatConfig1,...
                XYZDirectlyFromZMatConfig2)
            
            % //  Each boundary is an area where the points are closest to that atom.
            [CorrectionFunctions, XYZCoordsa0b0]=...
                determineMatlabAndOwnCoordsCorrectionFunctions(XYZDirectlyFromZMatConfig1,...
                XYZDirectlyFromZMatConfig2);
            
            atomSphericalCoordPositions = ...
                getSphericalCoordsRegionsForAtoms(CorrectionFunctions, ...
                XYZCoordsa0b0);
            
            plotAtomBoundariesOverlay(atomSphericalCoordPositions);
        end
        
        %------------------------------------------------------------
        %
        function self=translateByAlpha(self,deltaAlpha)
            
            self.alpha=self.alpha+deltaAlpha;
            self=self.centerPES;
        end
        
        %------------------------------------------------------------
        %
        function self=translateByBeta(self,deltaBeta)
            self.beta=self.beta+deltaBeta;
            self=self.centerPES;
        end
        
        %------------------------------------------------------------
        %
        function self=reflectAboutAlphaZero(self)
            self.alpha = -self.alpha+180;
            
            self=self.centerPES;
        end
        %------------------------------------------------------------
        %
        function self=reflectAboutBetaZero(self)
            self.beta = -self.beta + 360;
            self=self.centerPES;
        end
        
        %------------------------------------------------------------
        %
        function self=centerPES(self)
            self.alpha(self.alpha>180)=mod(self.alpha(self.alpha>180),180);
            self.alpha(self.alpha<0)=mod(self.alpha(self.alpha<0),180);
            
            self.beta(self.beta>=360)=mod(self.beta(self.beta>=360),360);
            self.beta(self.beta<0)=mod(self.beta(self.beta<0),360);
            
            [self.alpha,self.beta,self.energies] = ...
                self.fillInBeta360Entries(self.alpha,self.beta,self.energies);
        end
        
        %------------------------------------------------------------
        %
        function self=keepAlphaValues(self,alphaValues)
            
            edgePoints = self.alpha == 0 | self.alpha == 180 | ...
                self.beta == 0 | self.beta == 360;
            
            idxEntriesToKeep=...
                ismember(self.alpha,alphaValues) |  edgePoints;
            
            self=self.getEntries(idxEntriesToKeep);
        end
        
        %------------------------------------------------------------
        %
        function self=keepBetaValues(self, betaValues)
            
            edgePoints = self.alpha == 0 | self.alpha == 180 | ...
                self.beta == 0 | self.beta == 360;
            
            idxEntriesToKeep=...
                ismember(self.beta,betaValues) |  edgePoints;
            
            self=self.getEntries(idxEntriesToKeep);
        end
        %------------------------------------------------------------
        %
        function self=keepAlphaValuesThatAreMultiplesOf(self,BaseValue)
            
            edgePoints = self.alpha == 0 | self.alpha == 180 | ...
                self.beta == 0 | self.beta == 360;
            
            idxEntriesToKeep=...
                mod(self.alpha,BaseValue) == 0 |  edgePoints;
            
            self=self.getEntries(idxEntriesToKeep);
        end
        
        %------------------------------------------------------------
        %
        function self=keepBetaValuesThatAreMultiplesOf(self,BaseValue)
            
            edgePoints = self.alpha == 0 | self.alpha == 180 | ...
                self.beta == 0 | self.beta == 360;
            
            idxEntriesToKeep=...
                mod(self.beta,BaseValue) == 0 | edgePoints;
            
            self=self.getEntries(idxEntriesToKeep);
        end
        
        %------------------------------------------------------------
        %

        function [hContour,hNotation,MEPObject]=...
                plotMEP2(self,StartPoints,EndPoints,Color)
            
            if nargin==3
                Color='w';
            end
            % // update grid
            self=interpolatePESGrid(self);
            
            MEPObject=...
                minimumEnergyPathway(self.alphaGrid,self.betaGrid,self.energiesGrid);
            MEPObject=MEPObject.calculateMEP(StartPoints,EndPoints);
            
            % // plot contour
            hContour=figure;
            
            nContours=60;
            self.plotContourMapsView('nContours',nContours);

            hold on
            axis tight
            
            
            MEPObject.plotMEPLines(Color);
            
            % Only for one side, the others are repeats
            MinimaCutoff.Y=183;
            MinimaCutoff.X=177;
            MEPObject.plotMEPMinimaPoints(MinimaCutoff,Color);
            
            TSCutoff.Y=183;
            TSCutoff.X=177;
            MEPObject.plotMEPTSPoints(TSCutoff,Color);
            xlimits=get(gca,'xlim');
            ylimits=get(gca,'ylim');
            
            
            % // Notation only figure
            hNotation=figure;
            MEPObject.plotMEPLines();
            MEPObject.plotMEPMinimaPoints(MinimaCutoff,Color);
            MEPObject.plotMEPTSPoints(TSCutoff,Color);
            
            xlabel('')
            ylabel('')
            set(gca,'xticklabels','')
            set(gca,'yticklabels','')
            
            view(90,90)
            set(gca,'xlim',xlimits)
            set(gca,'ylim',ylimits)
            set(gcf, 'Position', get(0,'Screensize')) % maximise figure
            
            
        end
        
        %------------------------------------------------------------
        %
        function missingPoints = reportMissingPoints(self)
            allAlpha = [0:self.stepsizeAlpha:180].';
            allBeta = [0:self.stepsizeBeta:360].';

            
            nAlphaPoints=length(allAlpha);
            nBetaPoints=length(allBeta);
            
            allAlpha=repmat(allAlpha,[nBetaPoints,1]);
            allBeta=sort(repmat(allBeta,[nAlphaPoints,1]));
            
            allPoints=[allAlpha allBeta];
            
            foundPoints=[self.alpha self.beta];
            
            missingPoints=allPoints(~ismember(allPoints,foundPoints,'rows'),:);
            
            for i=1:size(missingPoints,1)
                fprintf(['submit_a' num2str(missingPoints(i,1)) '_b' ...
                    num2str(missingPoints(i,2)) '\n'])
            end
            
        end
        
        %------------------------------------------------------------
        %
        function self=interpolatePESGrid(self,gridsize)
            if nargin==1
                gridsize=0.5;
            end
            
            [self.alphaGrid,self.betaGrid,self.energiesGrid]=...
                interpolateGrid(self.alpha,self.beta,...
                self.energies,gridsize);
        end
    end
    
    methods (Hidden)
        
        
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
            
            CartesianCoords=AlphaBetaCoords.createXYZCoords;
            sphereCoords = CartesianCoords.coords;
            sphereCoordsEnergies= self.energiesGrid(:);
        end
        %------------------------------------------------------------
        %
        function self=removeEntries(self,idxEntriesToRemove)
            if isa(idxEntriesToRemove,'logical')
                % Indices in logical mask form
                self.alpha = self.alpha(~idxEntriesToRemove);
                self.beta = self.beta(~idxEntriesToRemove);
                self.energies = self.energies(~idxEntriesToRemove);
            else
                % Indices in counting form
                idxEntriesToKeep = 1:length(self.alpha);
                idxEntriesToKeep = idxEntriesToKeep(~ismembc(idxEntriesToKeep,...
                    sort(idxEntriesToRemove)));
                
                self.alpha = self.alpha(idxEntriesToKeep);
                self.beta = self.beta(idxEntriesToKeep);
                self.energies = self.energies(idxEntriesToKeep);
            end
        end
        
        
        %------------------------------------------------------------
        %
        function self=getEntries(self,indices)
            self.alpha = self.alpha(indices);
            self.beta = self.beta(indices);
            self.energies = self.energies(indices);
            
        end
        %------------------------------------------------------------
        %
        function self=plus(self,differentPES)
            self.alpha=[self.alpha; differentPES.alpha];
            self.beta=[self.beta; differentPES.beta];
            self.energies=[self.energies; differentPES.energies];
            
            if ~strcmpi(self.filename,differentPES.filename)
                warning('Combining PESs of potentially different jobs!')
            end
        end
        
        
        
        %------------------------------------------------------------
        %
        function pescoords=createSphericalPESCoordsObject(self)
            pescoords=sphericalPESCoords([self.alpha self.beta]);
            
        end
        
        %------------------------------------------------------------
        %
        function contourMatrixOfLines=calculateExoplanetContourLines(self,...
                n2DContourLevels)
            
            contourMatrixOfLines=plotContourMapsViewInRadians(self,...
                n2DContourLevels);
            
        end
        
        %------------------------------------------------------------
        %
        function transformAndPlotExoplanetContourLines(~,contourMatrixOfLines,...
                ColorMapSphere,RatioOf2DToSphereContourLevels)
            
            [~,contourMatrixNumCols] = size(contourMatrixOfLines);
            
            
            % Transform contour lines
            hFig=figure;
            k = 1;
            contourNumber = 1;
            contourLevel = contourMatrixOfLines(1,k);
            contourMin = contourLevel;
            while k < contourMatrixNumCols % Draw each contour line.
                kl = contourMatrixOfLines(2,k);
                v = k+1:k+kl;
                xCoord = cos(contourMatrixOfLines(2,v)).*cos(contourMatrixOfLines(1,v));
                yCoord = cos(contourMatrixOfLines(2,v)).*sin(contourMatrixOfLines(1,v));
                zCoord = sin(contourMatrixOfLines(2,v));
                if rounddec(contourMatrixOfLines(1,k),4) ~= rounddec(contourLevel,4)
                    contourNumber = contourNumber+1;
                    contourLevel = contourMatrixOfLines(1,k);
                end
                
                if mod(contourNumber,RatioOf2DToSphereContourLevels)==0
                    fill3(xCoord,yCoord,zCoord,ColorMapSphere(contourNumber,:))
                else
                    fill3(xCoord,yCoord,zCoord,ColorMapSphere(contourNumber,:),...
                        'EdgeColor','none')
                end
                
                hold on;
                k = k+kl+1;
            end
            
            
            caxis([contourMin contourLevel]);
            set(gca,'FontSize',14)
            xlabel('x / a.u.')
            ylabel('y / a.u.')
            zlabel('z / a.u.')
            
            
            axis tight;
            daspect([1 1 1]);
            
            % // colorbar
            colormap(ColorMapSphere);
            hcolorbar=colorbar;
            ylabel(hcolorbar,'Energy relative to most stable state / eV');
            hcolorbar.FontSize=14;
            
            for i=1:length(hFig.Children)
                if any(strfind(class(hFig.Children(i)),'ColorBar'))
                    hFig.Children(i).TickLabelInterpreter='latex';
                    hFig.Children(i).Label.Interpreter='latex';
                    break
                end
            end
        end
        
        %------------------------------------------------------------
        %
        function contourHandle=plotContourMapsViewInRadians(self,...
                nContourLevels)
            DegreesToRadians=2*pi/360;
            % Add 90 to alpha because of matlab definition that alpha is
            % the elevation, but in spherical coordinates, alpha goes down
            % from the vertical axis
            contourHandle = contourf(self.betaGrid*DegreesToRadians,...
                (self.alphaGrid+90)*DegreesToRadians,...
                self.energiesGrid,nContourLevels);
        end
        
        %------------------------------------------------------------
        %
        function shpEnergies = ...
                getPESShapeEnergiesAtPoints(self, PESShape, sphereCoords, ...
                sphereCoordsEnergies)
            
            IdxNearestAtoms = ...
                self.getNearestPointsFromPESShapeToPESSphereRep(PESShape, sphereCoords);
            
            shpEnergies = sphereCoordsEnergies(IdxNearestAtoms);
            
        end
    end
    
    methods (Hidden,Static)
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
        function [figHandle, shapeHandle] = ...
                plotPESshape(PESshape, PESshapeEnergies, faceAlpha)
            
            [triangulation, vertexCoords] = PESshape.boundaryFacets();
            
            % As tri refers to rows of vertexCoords, and not PESShape,
            % we need to get the indices of PESShape.Points
            % then we can get verticesEnergy
            [~,PESshapeIndices] = myismemberrows(vertexCoords,PESshape.Points);
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
        function atomShapes = ...
                getSubstrateAtomShapes(CenteredAuCoords, atomSphereSize, nSpherePoints)
            
            alphaShapeAlpha = 4;
            
            [smallsphereX,smallsphereY,smallsphereZ] = sphere(nSpherePoints);
            smallsphereX = smallsphereX(:);
            smallsphereY = smallsphereY(:);
            smallsphereZ = smallsphereZ(:);
            
            nAtoms = size(CenteredAuCoords,1);
            x = cell(1,nAtoms);
            y = cell(1,nAtoms);
            z = cell(1,nAtoms);
            atomShapes = cell(1,nAtoms);
            for i = 1:nAtoms
                AtomCoords = CenteredAuCoords(i,:);
                
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
        function hcolorbar = plotColorbar(hFig,colorbarLabel)
            hcolorbar = colorbar;
            ylabel(hcolorbar,colorbarLabel);
            hcolorbar.FontSize=14;

        end
        
        %------------------------------------------------------------
        %
        function [stepsizeAlpha,stepsizeBeta]=calcStepSize(alpha,beta)
            stepsizeBeta=min(diff(unique(sort(beta))));
            stepsizeAlpha=min(diff(unique(sort(alpha))));
        end
        
        %------------------------------------------------------------
        %
        function [alpha,beta,energies]=...
                extractAlphaBetaAndEnergiesFromFile(summaryFileName,DataFolder)
            
            % Get energies
            [EntryNames,energies]=textread([DataFolder '/' summaryFileName],'%s%n');
            
            % Get alpha and beta
            nEntries=length(EntryNames);
            alpha=zeros(nEntries,1);
            beta=zeros(nEntries,1);
            for i=1:nEntries
                CurrEntry=EntryNames{i};
                Parameters=strsplit(CurrEntry,{'_','.log'});
                Parameters(cellfun(@isempty,Parameters))=[];
                alpha(i)=str2double(Parameters{end-1}(2:end));
                beta(i)=str2double(Parameters{end}(2:end));
            end
        end
        
        %------------------------------------------------------------
        %
        function [alpha,beta,energies]=removeZeroEnergyEntries(alpha,beta,energies)
            % Remove zero entries
            iZeroEnergies=energies==0;
            alpha(iZeroEnergies)=[];
            beta(iZeroEnergies)=[];
            energies(iZeroEnergies)=[];
        end
        
        %------------------------------------------------------------
        %
        function NormalizedEnergies=normalizeEnergiesWithRespectToMinimum(energies)
            [minEnergy,~]=min(energies);
            
            NormalizedEnergies=energies-minEnergy;
        end
        
        %------------------------------------------------------------
        %
        function NormalizedEnergies=convertEnergiesToEVs(summaryFileName,NormalizedEnergies)
            if ~any(strfind(summaryFileName,'vasp'))
                NormalizedEnergies=NormalizedEnergies*27.114;
            end
            NormalizedEnergies=rounddec(NormalizedEnergies,4);
        end
        
        %------------------------------------------------------------
        %
        function  [alpha,beta,NormalizedEnergies]=...
                fillInBeta360Entries(alpha,beta,NormalizedEnergies)

            alpha360=alpha(beta==0);
            beta360=ones(size(alpha360))*360;
            NormalizedEnergies360=NormalizedEnergies(beta==0);
            
            alpha=[alpha(:) ;alpha360(:)];
            beta=[beta(:) ;beta360(:)];
            NormalizedEnergies=[NormalizedEnergies(:); NormalizedEnergies360(:)];
        end
        
        %------------------------------------------------------------
        %
        function [alpha,beta,NormalizedEnergies]=...
                fillInAlpha0Entries(alpha,beta,NormalizedEnergies,stepSizeBeta)
            
            % Fixed bug
            Energy_Alpha0=unique(rounddec(NormalizedEnergies(alpha==0),1));
               EnergyAlpha0=NormalizedEnergies(alpha==0);
            EnergyAlpha0=EnergyAlpha0(1);

            if length(Energy_Alpha0)>1
                error('All alpha zero entries should have same energies, but here they do not!')
            end
            
            Betas_Alpha0 = [0:stepSizeBeta:max(beta)].';
            Alphas_Alpha0 = zeros(size(Betas_Alpha0));
            Energies_Alpha0 = ones(size(Betas_Alpha0))*EnergyAlpha0;
            
            alpha=[alpha(:) ;Alphas_Alpha0(:)];
            beta=[beta(:) ;Betas_Alpha0(:)];
            NormalizedEnergies = [NormalizedEnergies(:); Energies_Alpha0(:)];
            
        end
        
        %------------------------------------------------------------
        %
        function [alpha,beta,NormalizedEnergies]=...
                fillInAlpha180Entries(alpha,beta,NormalizedEnergies,stepSizeBeta)
            
            Energy_Alpha180=unique(rounddec(NormalizedEnergies(alpha==180),2));
                  EnergyAlpha180=NormalizedEnergies(alpha==180);
            EnergyAlpha180=EnergyAlpha180(1);
            
            if isempty(Energy_Alpha180)
                return
            end
            
            if length(Energy_Alpha180)>1
                error('All alpha zero entries should have same energies, but here they do not!')
            end
            
            Betas_Alpha180=[0:stepSizeBeta:360].';
            Alphas_Alpha180=zeros(size(Betas_Alpha180))+180;
            Energies_Alpha180=ones(size(Betas_Alpha180))*EnergyAlpha180;
            
            alpha=[alpha(:) ;Alphas_Alpha180(:)];
            beta=[beta(:) ;Betas_Alpha180(:)];
            NormalizedEnergies=[NormalizedEnergies(:); Energies_Alpha180(:)];
            
        end
        
        %------------------------------------------------------------
        %
        function ColorMapSphere=...
                createExoplanetColorMap(n2DContourLevels,RatioOf2DToSphereContourLevels)
            
            
            ColorMapSphereTmp = ...
                parula(n2DContourLevels/RatioOf2DToSphereContourLevels);
            ColorMapSphere=zeros(n2DContourLevels,3);
            
            % Use the same colors for every RatioOf2DToSphereContourLevels
            for i=1:RatioOf2DToSphereContourLevels
                IdxOfSameColors=i:RatioOf2DToSphereContourLevels:n2DContourLevels-...
                    RatioOf2DToSphereContourLevels+i;
                ColorMapSphere(IdxOfSameColors,:)=ColorMapSphereTmp;
            end
        end
        
        %------------------------------------------------------------
        %
        function [x, y, z]= ...
                centerSphereOnAtomLocation(X, Y, Z, XYZCoordsa0b0, idxAtom, atomRadius)
            x = atomRadius*X+XYZCoordsa0b0(idxAtom, 1);
            y = atomRadius*Y+XYZCoordsa0b0(idxAtom, 2);
            z = atomRadius*Z+XYZCoordsa0b0(idxAtom, 3);
        end
        %------------------------------------------------------------
        %
        function [alpha, beta] = ...
                convertSphereCartesianCoordsToAlphaAndBeta(x,y,z)
            
            [beta, alpha]= cart2sph(x,y,z);
            alpha = alpha*180/pi;
            beta = beta*180/pi;
            
            beta(beta>360)=mod(beta(beta>360),360);
            beta(beta<0)=mod(beta(beta<0),360);
        end
        
        %------------------------------------------------------------
        %
        function currEnergies = getEnergiesOnSphere(alphaGridCorrected, ...
                betaGridCorrected, energiesGrid, currAlpha, currBeta)
            
            idxClosestElement = zeros(size(currBeta));
            
            for zz =1:numel(currBeta)
                idxClosestElement(zz) = getClosestElement(alphaGridCorrected, ...
                    betaGridCorrected,currAlpha(zz),currBeta(zz));
            end
            
            currEnergies = energiesGrid(idxClosestElement);
        end
    end
end
