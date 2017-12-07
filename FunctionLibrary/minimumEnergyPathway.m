classdef minimumEnergyPathway
    
    
    properties (GetAccess=public,SetAccess=private)
        alphaGrid
        betaGrid
        energiesGrid
        XCoords
        YCoords
        Minima
        TS
    end
    
    
    methods
        %------------------------------------------------------------
        %
        function  self=minimumEnergyPathway(alphaGrid,betaGrid,energiesGrid)
            
            sizes{1}=size(alphaGrid);
            sizes{2}=size(betaGrid);
            sizes{3}=size(energiesGrid);
            
            if ~isequaln(sizes{1},sizes{:})
                error('Inputs must have the same sizes!')
            end
            
            self.alphaGrid=alphaGrid;
            self.betaGrid=betaGrid;
            self.energiesGrid=energiesGrid;
            
        end
        
        %------------------------------------------------------------
        %
        function self=calculateMEP(self,StartPoints,EndPoints)
            
            [MEP_X,MEP_Y,MEP_Minima,MEP_TS]=...
                self.getMEPForAllInputPoints(StartPoints,EndPoints);
            
            MEP_Minima=self.removeDuplicatedMinimaPoints(MEP_Minima);
            
            MEP_TS=self.removeDuplicatedTSPoints(MEP_TS);
            
            self.XCoords=MEP_X;
            self.YCoords=MEP_Y;
            self.Minima=MEP_Minima;
            self.TS=MEP_TS;
        end
        
        %------------------------------------------------------------
        %
        function self=calculateMEPMultiTS(self,StartPoint,EndPoint)
            %Only can do one MEP at a time, but gives you all the TS's
            
            [MEP_X,MEP_Y,MEP_Minima,MEP_TS]=...
                self.getMEPAllTS(StartPoint,EndPoint);
            
            MEP_TS=self.removeDuplicatedTSPoints(MEP_TS);
            
            self.XCoords=MEP_X;
            self.YCoords=MEP_Y;
            self.Minima=MEP_Minima;
            self.TS=MEP_TS;
        end
        
        
        %------------------------------------------------------------
        %
        function plotMEPLines(self,Color)
            
            if nargin==1
                Color='w';
            end
            % // plot lines
            if ~iscell(self.XCoords)
                self.XCoords=num2cell(self.XCoords,1);
                self.YCoords=num2cell(self.YCoords,1);
            end
            
            for i=1:length(self.XCoords)
                
                plotbrokenlines(self.XCoords{i},self.YCoords{i},'linestyle','-','color',Color,'linewidth',2)


%                 Ydis=abs(diff(self.YCoords{i}));
%                 skipped=[0 (Ydis>100)'];
%                 blocks = cumsum(skipped);
%                 if sum(blocks)>0
%                     
%                     indexofchange=diff(blocks);
%                     point=find(indexofchange==1);
%                     indexofchange(point+1:end+1)=indexofchange(point:end);
%                     
%                     IndexofFirst=find(indexofchange==1,1);
%                     FirstBeta=self.YCoords{i}(IndexofFirst);
%                     SecondBeta=self.YCoords{i}(IndexofFirst+1);
%                     
%                     if FirstBeta<180
%                         FirstBeta=FirstBeta+360;
%                         flag=1;
%                     elseif SecondBeta<180
%                         SecondBeta=SecondBeta+360;
%                         flag=2;
%                     else
%                         error('Something is wrong with the points')
%                     end
%                     
%                     Yvalues=self.XCoords{i}(indexofchange==1);
%                     slope=diff(Yvalues)/(SecondBeta-FirstBeta);
%                     
%                     if flag==1
%                         NewYvalue=slope*(360-FirstBeta)+Yvalues(1);
%                     elseif flag==2
%                         NewYvalue=slope*(360-SecondBeta)+Yvalues(2);
%                     end
%                     
%                     self.XCoords{i}(point+3:end+2)=self.XCoords{i}(point+1:end);
%                     self.YCoords{i}(point+3:end+2)=self.YCoords{i}(point+1:end);
%                     
%                     self.XCoords{i}(point+1:point+2)=NewYvalue;
%                     
%                     if flag==1
%                         self.YCoords{i}(point+1:point+2)=[0 360];
%                     elseif flag==2
%                         self.YCoords{i}(point+1:point+2)=[360 0];
%                     end
%                     
%                     Ydis=abs(diff(self.YCoords{i}));
%                     skipped=[0 (Ydis>100)'];
%                     blocks = cumsum(skipped);
%                     
%                 end
%                 
%                 for j=0:blocks(end)
%                     plot(self.XCoords{i}(blocks==j),self.YCoords{i}(blocks==j),...
%                         'linestyle','-','color',Color,'linewidth',2);
%                     hold on
%                 end

            end
            
        end
        
        %------------------------------------------------------------
        %
        function plotMEPTSPoints(self,TSCutoff,Color)
            
            if nargin==2
                Color='w';
            end
            
            line(self.TS(:,1),self.TS(:,2),'linestyle','none','linewidth',1.5,'marker','v',...
                'color',Color,'MarkerSize',6,'MarkerFaceColor',Color);
            
            % Labels only for one side, the others are repeats
            iMEP_TSOneSide=self.TS(:,2)<TSCutoff.Y;
            MEP_TSOneSide=self.TS(iMEP_TSOneSide,:);
            
            
            labelsTS=strsplit(num2str(1:size(MEP_TSOneSide,1)),' ');
            labelsTS=cellfun(@(x) [x '*'],labelsTS,'UniformOutput',0);
            
            for i=1:length(labelsTS)
                text(MEP_TSOneSide(i,1), MEP_TSOneSide(i,2), labelsTS{i}, ...
                    'VerticalAlignment','top', ...
                    'HorizontalAlignment','right',...
                    'FontSize',14,'Color',Color,...
                    'FontName','Helvetica')
                hold on
                drawnow
                
            end
            
            % -- Label energies
            TSEnergies=self.getMEPEnergies(MEP_TSOneSide(:,1),...
                MEP_TSOneSide(:,2));
            
            labelsts_energy=strsplit(num2str(TSEnergies,'%3.2f '),' ');
            
            % Make the labels for energy lower than labels for identifier
            YOffset=3;
            text(MEP_TSOneSide(:,1)+YOffset, MEP_TSOneSide(:,2), ...
                labelsts_energy, ...
                'VerticalAlignment','top', 'HorizontalAlignment','left', ...
                'FontSize',14,'Color',Color,'FontName','Helvetica','FontAngle',...
                'italic');
            
        end
        
        %------------------------------------------------------------
        %
        function plotMEPMinimaPoints(self,MinimaCutoff,Color)
            
            if nargin==2
                Color='w';
            end
            
            line(self.Minima(:,1),self.Minima(:,2),'linestyle','none','marker','.',...
                'color',Color,'MarkerSize',24);
            % -- Label with identifying numbers
            
            
            iMEP_MinimaOneSide=self.Minima(:,2)<MinimaCutoff.Y;
            MEP_MinimaOneSide=self.Minima(iMEP_MinimaOneSide,:);
            labelsminima=strsplit(num2str(1:size(MEP_MinimaOneSide,1)),' ');
            
            for i=1:length(labelsminima)
                text(MEP_MinimaOneSide(i,1), MEP_MinimaOneSide(i,2), labelsminima{i}, ...
                    'VerticalAlignment','top', 'HorizontalAlignment','left', ...
                    'FontSize',14,'Color',Color,'FontName','Helvetica')
                hold on
                drawnow
            end
            % -- Label energies
            MinimaEnergies=self.getMEPEnergies(MEP_MinimaOneSide(:,1),...
                MEP_MinimaOneSide(:,2));
            
            labelsminima_energy=strsplit(num2str(MinimaEnergies,'%3.2f '),' ');
            
            % Make the labels for energy lower than labels for identifier
            YOffset=3;
            text(MEP_MinimaOneSide(:,1)+YOffset, MEP_MinimaOneSide(:,2), ...
                labelsminima_energy, ...
                'VerticalAlignment','top', 'HorizontalAlignment','left', ...
                'FontSize',14,'Color',Color,'FontName','Helvetica','FontAngle',...
                'italic')
        end
        
    end
    
    methods (Hidden)
        %------------------------------------------------------------
        %
        function [MEP_X,MEP_Y,MEP_Minima,MEP_TS]=...
                getMEPForAllInputPoints(self,StartPoints,EndPoints)
            % Loop over all start and end points to get MEPs for each pair of points
            
            PlaceHolderCoords=-1;
            nPathways=length(StartPoints);
            MEP_Minima=zeros(nPathways*2,2)+PlaceHolderCoords;
            MEP_TS=zeros(nPathways,2)+PlaceHolderCoords;
            MEP_X=cell(1,nPathways);
            MEP_Y=cell(1,nPathways);
            for i=1:nPathways
                
                [MEP_X,MEP_Y,MEP_Minima,MEP_TS]=...
                    self.calcFullTrajectory(StartPoints,EndPoints,i,...
                    MEP_X,MEP_Y,MEP_Minima,MEP_TS);
            end
            
            MEP_TS(ismember(MEP_TS,[PlaceHolderCoords PlaceHolderCoords],'rows'),:)=[];
            MEP_Minima(ismember(MEP_Minima,[PlaceHolderCoords PlaceHolderCoords],'rows'),:)=[];
            
        end
        
        function [MEP_X,MEP_Y,MEP_Minima,MEP_TS]=...
                getMEPAllTS(self,StartPoints,EndPoints)
            % Loop over all start and end points to get MEPs for each pair of points
            
            [MEP_X,MEP_Y,MEP_Minima,MEP_TS]=...
                self.calcFullTrajectoryMultTS(StartPoints,EndPoints);
            
            %             MEP_TS(ismember(MEP_TS,[PlaceHolderCoords PlaceHolderCoords],'rows'),:)=[];
            %             MEP_Minima(ismember(MEP_Minima,[PlaceHolderCoords PlaceHolderCoords],'rows'),:)=[];
            %
        end
        
        
        %------------------------------------------------------------
        %
        function  [MEP_X,MEP_Y,MEP_Minima,MEP_TS]=  calcFullTrajectoryMultTS(self,...
                StartPoints,EndPoints)
            
            CurrStartPoint=StartPoints;
            CurrEndPoint=EndPoints;
            
            [CurrMEP_Coords,CurrMinima_Coords,CurrTS_Coords]=...
                self.getMEPengineMultiTS(CurrStartPoint,CurrEndPoint);
            
            
            % Save MEP pathway
            MEP_X=CurrMEP_Coords(:,1);
            MEP_Y=CurrMEP_Coords(:,2);
            
            
            MEP_Minima=CurrMinima_Coords;
            MEP_TS=CurrTS_Coords;
            
        end
        
        %------------------------------------------------------------
        %
        function  [MEP_X,MEP_Y,MEP_Minima,MEP_TS]=  calcFullTrajectory(self,...
                StartPoints,EndPoints,i,MEP_X,MEP_Y,MEP_Minima,MEP_TS)
            
            CurrStartPoint=StartPoints{i};
            CurrEndPoint=EndPoints{i};
            
            [CurrMEP_Coords,CurrMinima_Coords,CurrTS_Coords]=...
                self.getMEPengine(CurrStartPoint,CurrEndPoint);
            
            
            % Save MEP pathway
            MEP_X{i}=CurrMEP_Coords(:,1);
            MEP_Y{i}=CurrMEP_Coords(:,2);
            
            
            % Check if both minima are the same
            Acc=0.5;
            UniqueMinima=unique(round(CurrMinima_Coords*Acc)/Acc,'rows');
            nUniqueMinima=size(UniqueMinima,1);
            
            if nUniqueMinima==1
                
                % For "single point" calcs, don't save TS, as both ends of the string
                % are the same minima
                
                % If is single point, use the same minima point for both so as to
                % reduce trouble when rounding etc.
                
                
                MEP_Minima(i*2-1:i*2,:)=[CurrMinima_Coords(1,:);CurrMinima_Coords(1,:)];
            else
                % -- Save minima and TS coords
                MEP_Minima(i*2-1:i*2,:)=CurrMinima_Coords;
                MEP_TS(i,:)=CurrTS_Coords;
            end
        end
        
        
        %------------------------------------------------------------
        %
        function [MEP_Coords,Minima_Coords,TS_Coords]=...
                getMEPengineMultiTS(self,CurrStartPoint,CurrEndPoint,varargin)
            
            X1=sind(CurrStartPoint{1}(1)).*cosd(CurrStartPoint{1}(2));
            Y1=sind(CurrStartPoint{1}(1)).*sind(CurrStartPoint{1}(2));
            Z1=cosd(CurrStartPoint{1}(1));
            X2=sind(CurrEndPoint{1}(1)).*cosd(CurrEndPoint{1}(2));
            Y2=sind(CurrEndPoint{1}(1)).*sind(CurrEndPoint{1}(2));
            Z2=cosd(CurrEndPoint{1}(1));
            
            dist=sqrt((X1-X2)^2+(Y1-Y2)^2+(Z1-Z2)^2);
            nimages=ceil(dist*55);
            if nimages < 3
                nimages=3;
            end
            
            
            if abs(CurrStartPoint{1}(2)-CurrEndPoint{1}(2))>180
                
                PESshifted=TranslateByBeta(self,180);
                %CurrStartPoint{1}(2)=mod(CurrStartPoint{1}(2)+180,360);
                ShiftedBetaStart=mod(CurrStartPoint{1}(2)+180,360);
                ShiftedBetaEnd=mod(CurrEndPoint{1}(2)+180,360);
                %CurrEndPoint{1}(2)=mod(CurrEndPoint{1}(2)+180,360);
                
                [MEP_XCoords,MEP_YCoords]=...
                    ztsMueller(PESshifted.alphaGrid,PESshifted.betaGrid,PESshifted.energiesGrid,...
                    CurrStartPoint{1}(1),...
                    ShiftedBetaStart,CurrEndPoint{1}(1),ShiftedBetaEnd,'nodraw','nimages',nimages);
                
%                 [MEP_XCoords,MEP_YCoords]=...
%                     ztsMueller(PESshifted.alphaGrid,PESshifted.betaGrid,PESshifted.energiesGrid,...
%                     CurrStartPoint{1}(1),...
%                     ShiftedBetaStart,CurrEndPoint{1}(1),ShiftedBetaEnd,'nimages',nimages);
%                 
                %self=TranslateByBeta(self,180);
                MEP_YCoords=mod(MEP_YCoords+180,360);
                %MEP_XCoords=mod(MEP_XCoords+180,360);
                %CurrStartPoint{1}(2)=mod(CurrStartPoint{1}(2)+180,360);
                %CurrEndPoint{1}(2)=mod(CurrEndPoint{1}(2)+180,360);
                
            else
                [MEP_XCoords,MEP_YCoords]=...
                    ztsMueller(self.alphaGrid,self.betaGrid,self.energiesGrid,...
                    CurrStartPoint{1}(1),...
                    CurrStartPoint{1}(2),CurrEndPoint{1}(1),CurrEndPoint{1}(2),'nodraw','nimages',nimages);
            end
            
            % // Get PES from zero temperature string method
            
            MEP_Coords=[MEP_XCoords(:) MEP_YCoords(:)];
            
            % // Get energies of the found MEP Path
            CurrMEP_Energies=self.getMEPEnergies(MEP_XCoords,MEP_YCoords);
            
            % // Find coordinates of the TS
            [~,Idx]=findpeaks(CurrMEP_Energies);
            [~,IdxMin]=findpeaks(-CurrMEP_Energies);
            
            %Compiles the Minima
            Minima_Coords=[MEP_Coords(1,:);MEP_Coords(IdxMin,:);MEP_Coords(end,:)];
            
            TS_Coords=[MEP_Coords(Idx,:)];
        end
        
        function [MEP_Coords,Minima_Coords,TS_Coords]=...
                getMEPengine(self,CurrStartPoint,CurrEndPoint,varargin)
            
            % // Get PES from zero temperature string method
            [MEP_XCoords,MEP_YCoords]=...
                ztsMueller(self.alphaGrid,self.betaGrid,self.energiesGrid,...
                CurrStartPoint(1),...
                CurrStartPoint(2),CurrEndPoint(1),CurrEndPoint(2),'nodraw');
            
            MEP_Coords=[MEP_XCoords(:) MEP_YCoords(:)];
            Minima_Coords=MEP_Coords([1,end],:);
            
            % // Get energies of the found MEP Path
            CurrMEP_Energies=self.getMEPEnergies(MEP_XCoords,MEP_YCoords);
            
            % // Find coordinates of the TS
            [~,Idx]=max(CurrMEP_Energies);
            TS_Coords=[MEP_Coords(Idx,:)];
            
        end
        %------------------------------------------------------------
        %
        function Energies=getMEPEnergies(self,XCoords,YCoords)
            
            iElement=zeros(1,length(XCoords));
            for z=1:length(XCoords)
                CurrX=XCoords(z);
                CurrY=YCoords(z);
                iElement(z)=getClosestElement(self.alphaGrid,self.betaGrid,CurrX,CurrY);
            end
            Energies=self.energiesGrid(iElement);
        end
        
        %------------------------------------------------------------
        %
        function self=TranslateByBeta(self,delb)
            [a,~]=size(self.alphaGrid);
            tempE=self.energiesGrid;
            %tempB=self.betaGrid;
            for i=1:a
                self.energiesGrid(i,:)=tempE(mod(i-1-delb*(a-1)/360,a-1)+1,:);
                %self.betaGrid(i,:)=tempB(mod(i-1-delb*(a-1)/360,a-1)+1,:);
            end
            self.YCoords=mod(self.YCoords+delb,360);
            
            if ~isempty(self.Minima)
            self.Minima(:,2)=mod(self.Minima(:,2)+delb,360);
            end
            
            if ~isempty(self.TS)
            self.TS(:,2)=mod(self.TS(:,2)+delb,360);
            end
            
            
        end
    end
    
    methods (Hidden,Static)
        
        %------------------------------------------------------------
        %
        function MEP_Minima=removeDuplicatedMinimaPoints(MEP_Minima)
            % Add t4 sites minima
            %     MEP_Minima=[MEP_Minima; 7 55.75; 7 305.3; 112,353.8; 112,6.25];
            MEP_MinimaOri=MEP_Minima;
            
            % Round minima coords so that we can use unique
            Acc=2;
            MEP_MinimaRounded=round(MEP_Minima/Acc)*Acc;
            [~,iNewRowsMinima,~]=unique(MEP_MinimaRounded,'rows');
            MEP_Minima=MEP_Minima(iNewRowsMinima,:);
        end
        
        %------------------------------------------------------------
        %
        function MEP_TS=removeDuplicatedTSPoints(MEP_TS)
            MEP_TSOri=MEP_TS;
            
            MEP_TSRounded=round(MEP_TS);
            [~,iNewRowsTS,~]=unique(MEP_TSRounded,'rows');
            MEP_TS=MEP_TS(iNewRowsTS,:);
        end
        
    end
end