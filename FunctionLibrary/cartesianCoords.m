classdef cartesianCoords < coordinates
    
    
    methods
        function cartesianCoords=cartesianCoords(coords)
            % //
            if size(coords,2) ~= 3
                error('Cartesian coords can only have 3 dimensions!')
            end
            cartesianCoords@coordinates(coords);
        end
        
        function [iClosestElementOfFirstSetOfCoords,iClosestElementOfSecondSetOfCoords]=...
                getClosestElementsWith(cartesianCoordsOne,cartesianCoordsTwo)
            
            distances=cartesianCoordsOne.getDistancesWith(cartesianCoordsTwo);
            [~,minDistanceIdx]=min(distances(:));
            
            % Get only the first one if there are multiple
            minDistanceIdx=minDistanceIdx(1);
            
            sizeDistances=size(distances);
            
            % //
            [iClosestElementOfFirstSetOfCoords,iClosestElementOfSecondSetOfCoords]=...
                ind2sub(sizeDistances,minDistanceIdx);
        end
        
        
        
        function distances=getDistancesWith(cartesianCoordsOne,cartesianCoordsTwo)
            % // If both sets of coords are matrices, then distances is a
            % // matrix with elements i,j corresponding to row i of the
            % // first matrix and row j of the second matrix
            distances=...
                euclideanDist(cartesianCoordsOne.coords,cartesianCoordsTwo.coords);
            
        end
        
        function normalizedCartesianCoords=normalizeRows(cartesianCoords)
            % // Normalizing each row individually to 1
            coords=cartesianCoords.coords;
            normInputRows=sqrt(sum(coords.^2,2));
            normalizedRows=bsxfun(@times,coords,1./normInputRows);
            
            %//
            normalizedCartesianCoords=...
                cartesianCoords.createObjectOfSimilarClass(normalizedRows);
        end
    end
end