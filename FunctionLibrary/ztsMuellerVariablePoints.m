function [MEP_XCoords,MEP_YCoords]=ztsMuellerVariablePoints(alphaGrid,betaGrid,energiesGrid,...
    CurrStartPoint1,CurrStartPoint2,CurrEndPoint1,CurrEndPoint2,varargin)
%Not setup yet, but I neat idea



options.NoDraw=0;
options.nImages=30;

i=1;
if ~isempty(varargin)
    while i<=length(varargin)
        if strcmpi(varargin{i},'nodraw')
            options.NoDraw=1;
            i=i+1;
        elseif strcmpi(varargin{i},'nimages')
            options.nImages = varargin{i+1};
            i=i+2;
        elseif isempty(varargin{i})
            i=i+1;
        else
            error('No such flag "%s"!',varargin{i})
        end
    end
end



set(0,'DefaultTextFontName','TimesRoman')
set(0,'DefaultAxesFontSize',16)



% // Options
options.nStepmax = 5000; % max number of iterations
options.nStepplot = 50; % frequency of plotting
options.PlottingFlag1 = 1;% plot string every nstepplot if flag1 = 1
options.WriteMovieFlag = 0;
options.StoppingTol1 = 1e-7; % parameter used as stopping criterion
options.TimeStep = 10e-0;% time-step; independent of nImages


[MEP_XCoords,MEP_YCoords]=initializeString(CurrStartPoint1,CurrEndPoint1,...
    CurrStartPoint2,CurrEndPoint2,options);

[gradientX,gradientY]=getGradient(alphaGrid,betaGrid,energiesGrid);

%plotInitialString(xx,yy,V1,x,y,NoDraw);

whitedots=initializePlotOfStringEvolution(alphaGrid,betaGrid,energiesGrid,...
    MEP_XCoords,MEP_YCoords,options);

mov=openMovie(options);

[MEP_XCoords,MEP_YCoords,tol,nstep,dVx,dVy]=runZTSMethod(MEP_XCoords,...
    MEP_YCoords,alphaGrid,betaGrid,gradientX,gradientY,mov,whitedots,...
    options);

printOutput(options,tol,nstep)

%plotFinalString(xx,yy,V1,x,y)


estimateErrorBetweenDiscreteAndActualMEP(MEP_XCoords,MEP_YCoords,dVx,dVy)

%reinterpolateStringToGetAccurateEnergy

end
%------------------------------------------------------------
%
function [x,y]=initializeString(xa,xb,ya,yb,options)

%[MEP_XCoords,MEP_YCoords]=initializeString(CurrStartPoint1,CurrEndPoint1
%,CurrStartPoint2,CurrEndPoint2,options);

% // String initialization
% Start and end points of the initial string
% notice that they do NOT have to be at minima of V -- the method finds
% those automatically


% initialization of string. each image is represented by with coords stored
% in "x" and "y"
g1 = linspace(0,1,options.nImages);
%x = (xb-xa)*g1+xa;
%y = (x-xa)*(yb-ya)/(xb-xa)+ya;

x=linspace(xa,xb,options.nImages);%linear interpolation
y=linspace(ya,yb,options.nImages);

dx = x-circshift(x,[0 1]);
dy = y-circshift(y,[0 1]);
dx(1) = 0;
dy(1) = 0;
lxy = cumsum(sqrt(dx.^2+dy.^2));
lxy = lxy/lxy(options.nImages);
x = interp1(lxy,x,g1);
y = interp1(lxy,y,g1);
end
%------------------------------------------------------------
%
function [gradientX,gradientY]=getGradient(alphaGrid,BetaGrid,energiesGrid)
UniqueXX=unique(alphaGrid);
AccXX=UniqueXX(2)-UniqueXX(1);

UniqueYY=unique(BetaGrid);
AccYY=UniqueYY(2)-UniqueYY(1);


[gradientX,gradientY] = gradient(energiesGrid,AccXX,AccYY);
end
%------------------------------------------------------------
%
function plotInitialString(xx,yy,V1,xi,yi,NoDraw)

if ~NoDraw
    
    
    % Initial string plot
    
    figure(1);clf;
    contourf(xx,yy,min(V1,200),60);
    hold on
    axis tight
    plot(xi,yi,'.-w','MarkerSize',14)
    %set(gca,'XTick',-1.5:.5:1,'YTick',0:.5:2);
    xlabel('x','FontAngle','italic');
    ylabel('y','FontAngle','italic');
    title('Initial string');
    drawnow
    
end
end
%------------------------------------------------------------
%
function whitedots=initializePlotOfStringEvolution(xx,yy,V1,x,y,options)
% // String evolution plot
whitedots=[];
if ~options.NoDraw
    figure(1);
    clf;
    axis tight;
    
    contourf(xx,yy,min(V1,200),60);
    whitedots = line(x,y,'linestyle','none','marker','.','color','w','MarkerSize',14);
    %set(gca,'XTick',-1.5:.5:1,'YTick',0:.5:2);
    xlabel('x','FontAngle','italic');
    ylabel('y','FontAngle','italic');
    title('String evolution')
    drawnow
end

end
%------------------------------------------------------------
%
function mov=openMovie(options)

mov=[];
if options.WriteMovieFlag
    
    view(90,90) %rotate to correct viewing angle
    set(gcf, 'Position', get(0,'Screensize')) % maximise figure
    
    
    mov = VideoWriter('Mov.avi','MPEG-4');
    mov.FrameRate = 15;  % Default 30
    mov.Quality = 120;    % Default 75
    open(mov);
    F = getframe(gcf);
    
    % Write frame to avi file
    for i=1:12
        writeVideo(mov,F);
    end
end
end
%------------------------------------------------------------
%
function [MEP_XCoords,MEP_YCoords,tol,nstep,GradientXNear,GradientYNear]=...
    runZTSMethod(MEP_XCoords,MEP_YCoords,alphaGrid,betaGrid,GradientX,...
    GradientY,mov,whitedots,options)
for nstep=1:options.nStepmax
    
    xBeforeStep = MEP_XCoords;
    yBeforeStep = MEP_YCoords;
    
    [GradientXNear,GradientYNear]=GetNearestGradient(alphaGrid,betaGrid,...
        MEP_XCoords,MEP_YCoords,GradientX,GradientY);
    
    [MEP_XCoords,MEP_YCoords]=takeOptimizationStep(MEP_XCoords,MEP_YCoords,...
        GradientXNear,GradientYNear,options);
    
    [MEP_XCoords,MEP_YCoords]=reparameterizeString(MEP_XCoords,MEP_YCoords);
    
    plotIntermediateStrings(MEP_XCoords,MEP_YCoords,whitedots,mov,nstep,options);
    
    [tol,Stop]=checkingForStoppingCriteria(MEP_XCoords,xBeforeStep,MEP_YCoords,yBeforeStep,options);
    
    if Stop
        break
    end
end
end
%------------------------------------------------------------
%

function [GradientXNear,GradientYNear]=GetNearestGradient(alphaGrid,...
    betaGrid,MEP_XCoords,MEP_YCoords,GradientX,GradientY)
% // Calculation of the x and y-components of the force
% dVx and dVy respectively. These are vectors, which store the forces
% on each image.

nImages=length(MEP_XCoords);
iElement=zeros(1,nImages);
for z=1:nImages
    iElement(z)=getClosestElement(alphaGrid,betaGrid,MEP_XCoords(z),MEP_YCoords(z));
end
GradientXNear=GradientX(iElement);
GradientYNear=GradientY(iElement);
end

%------------------------------------------------------------
%
function [MEP_XCoords,MEP_YCoords]=takeOptimizationStep(MEP_XCoords,...
    MEP_YCoords,GradientXNear,GradientYNear,options)

MEP_XCoords = MEP_XCoords - options.TimeStep*GradientXNear;
MEP_YCoords = MEP_YCoords - options.TimeStep*GradientYNear;

% Keep bounds
MEP_YCoords(MEP_YCoords<0)=0;
MEP_XCoords(MEP_XCoords<0)=0;
MEP_YCoords(MEP_YCoords>357)=357;
MEP_XCoords(MEP_XCoords>177)=177;
end
%------------------------------------------------------------
%
function [MEP_XCoords,MEP_YCoords]=reparameterizeString(MEP_XCoords,MEP_YCoords)
[ArcLength]=calcArcLength(MEP_XCoords,MEP_YCoords);
[MEP_XCoords,MEP_YCoords]=interpolateNewPoints(ArcLength,MEP_XCoords,MEP_YCoords);
end
%------------------------------------------------------------
%
function [ArcLength]=calcArcLength(MEP_XCoords,MEP_YCoords)
nImages=length(MEP_XCoords);
% Calculate arc length, lxy
dx = MEP_XCoords-circshift(MEP_XCoords,[0 1]);
dy = MEP_YCoords-circshift(MEP_YCoords,[0 1]);
dx(1) = 0;
dy(1) = 0;
ArcLength = cumsum(sqrt(dx.^2+dy.^2));
ArcLength = ArcLength/ArcLength(nImages);
end
%------------------------------------------------------------
%
function [x,y]=interpolateNewPoints(lxy,x,y)
% Next we use interpolation to obtain the new points i
% n+1 at the uniform grid points i=i/N. This can be done,
% for example, by using cubic spline interpolation for the
% data
nImages=length(x);
g1 = linspace(0,1,nImages);

x = interp1(lxy,x,g1);
y = interp1(lxy,y,g1);
end

%------------------------------------------------------------
%

function plotIntermediateStrings(MEP_XCoords,MEP_YCoords,whitedots,mov,nstep,options)
if ~options.NoDraw && and(options.PlottingFlag1 == 1,mod(nstep,options.nStepplot) == 0)
    set(whitedots,'xdata',MEP_XCoords,'ydata',MEP_YCoords)
    drawnow
    pause(0.025)
    
    
    if options.WriteMovieFlag
        
        F = getframe(gcf);
        
        % Write frame to avi file
        for i=1:12
            writeVideo(mov,F);
        end
    end
end
end

%------------------------------------------------------------
%
function [tol,Stop]=checkingForStoppingCriteria(MEP_XCoords,xBeforeStep,...
    MEP_YCoords,yBeforeStep,options)

nImages=length(MEP_XCoords);
tol = (norm(MEP_XCoords-xBeforeStep)+norm(MEP_YCoords-yBeforeStep))/nImages;

if tol <= options.StoppingTol1
    Stop=true;
else
    Stop=false;
end
end

%------------------------------------------------------------
%
function printOutput(options,tol,nstep)
fprintf('\n')
fprintf('\n')
fprintf('ZTS calculation with %d images\n',options.nImages)
if tol > options.StoppingTol1
    fprintf('The calculation failed to converge after %d iterations\n',nstep)
else
    fprintf('The calculation terminated after %d iterations\n',nstep)
end
end

%------------------------------------------------------------
%
function plotFinalString(xx,yy,V1,x,y)
figure(2);clf;
contourf(xx,yy,min(V1,200),60);
hold on
plot(x,y,'.-w','MarkerSize',14)
set(gca,'XTick',-1.5:.5:1,'YTick',0:.5:2);
xlabel('x','FontAngle','italic');
ylabel('y','FontAngle','italic');
title('Final string');
drawnow
end

%------------------------------------------------------------
%
function estimateErrorBetweenDiscreteAndActualMEP(x,y,dVx,dVy)
% Energy along MEP
tx = circshift(x,[0 -1])-circshift(x,[0 1]);
ty = circshift(y,[0 -1])-circshift(y,[0 1]);

% potential computed as integral of projection of gradV on string tangent
Vz=cumtrapz(tx.*dVx + ty.*dVy);
Vz=0.5*Vz;
Vz=Vz-min(Vz);


ntxy = sqrt(tx.*tx+ty.*ty);
tx = tx./ntxy;
ty = ty./ntxy;


% err is an estimate of the error between disrectized MEP and actual one
% err scale as 1/n1
nImages=length(x);
err = trapz(1-(tx.*dVx+ty.*dVy).^2./(dVx.*dVx+dVy.*dVy))/(nImages-1);
fprintf('Estimate of difference between discretized MEP and actual one: %f\n',err)
end

%------------------------------------------------------------
%
function reinterpolateStringToGetAccurateEnergy
% // Reinterpolate string with lots of points to get accurate energy along it
%{
g0 = linspace(0,1,1e3);
x0 = interp1(lxy,x,g0);
y0 = interp1(lxy,y,g0);
ee = AA(1)*exp(aa(1)*(x0-XX(1)).^2+bb(1)*(x0-XX(1)).*(y0-YY(1))+cc(1)*(y0-YY(1)).^2);
V0=ee;
dVx = (2*aa(1)*(x0-XX(1))+bb(1)*(y0-YY(1))).*ee;
dVy = (bb(1)*(x0-XX(1))+2*cc(1)*(y0-YY(1))).*ee;
for j=2:4
    ee = AA(j)*exp(aa(j)*(x0-XX(j)).^2+bb(j)*(x0-XX(j)).*(y0-YY(j))+cc(j)*(y0-YY(j)).^2);
    V0 = V0 + ee;
    dVx = dVx + (2*aa(j)*(x0-XX(j))+bb(j)*(y0-YY(j))).*ee;
    dVy = dVy + (bb(j)*(x0-XX(j))+2*cc(j)*(y0-YY(j))).*ee;
end;

% Normalize potential energy
V0=V0-min(V0);

figure(3);clf;
hold on
title('Energy along MEP');
plot(g1,Vz)
plot(g0,V0,'r')
box on
legend('Thermodynamic integration along string','Exact','Location','South')
legend('boxoff')
%}
end
