classdef sphericalPESCoords < periodicCoords
    % // Coords in alpha and beta
    % // Beta is longitude (matlab: azimuth)
    % // Alpha is angle coming down from zaxis
    % // (matlab: elevation = 90-alpha)
    % //         |   /
    % //         |^1/   1.alpha
    % //         | /      2.elevation (used by matlab sph2cart)
    % //         |/_)2 _ _ _ _
    
    properties (GetAccess=public,SetAccess=private)
        radius
    end
    
    
    methods
        function sphericalPESCoords=sphericalPESCoords(coords)
            alphaBetaPeriodicity=[180,360];
            sphericalPESCoords@periodicCoords(alphaBetaPeriodicity,coords)
            
            sphericalPESCoords.radius=1;
        end
        
        function self=setRadius(self,radius)
            
            NumPoints=size(self.coords,1);
            
            if length(radius) == 1 || length(radius) ==  NumPoints
                self.radius=radius;
            else
                error('Radius must be either scalar, or vector containing as many points as coords!')
            end
        end
        
    
        
        function cartesianCoordsObj=createXYZCoords(object)
            
            if isempty(object.radius)
                error('Radius property of sphericalPESCoords must be first set!')
            end
            DegreesToRadians=2*pi/360;
            
            Coords=object.coords;
          %  Elevation=(90-Coords(:,1))*DegreesToRadians;
            Elevation=(Coords(:,1))*DegreesToRadians;
            Azimuth=(Coords(:,2))*DegreesToRadians;
            Radius=object.radius;
            
            [xCoords,yCoords,zCoords]=...
                sph2cart(Azimuth,Elevation,Radius);
            
           cartesianCoordsObj=...
                cartesianCoords([xCoords yCoords zCoords]);
        end
    end
end