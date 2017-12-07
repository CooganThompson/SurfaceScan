classdef periodicCoords < coordinates
    
    properties (GetAccess=public,SetAccess=public)
        periodicity
    end
    
    methods
        function periodicCoords=periodicCoords(periodicity,coords)
            
            % // Inputs
            periodicCoords@coordinates(coords);
            
            periodicCoords.periodicity=periodicity;
            
            if length(periodicity) ~= size(coords,2)
                error('Length of inputted periodicity and coords must match!')
            end
            
            % //
            periodicCoords=periodicCoords.reduceToUnitCell;
        end
        
        
        function periodicCoords=reduceToUnitCell(periodicCoords)
            
            for i=1:periodicCoords.dimensionality
                PeriodicityOfCurrDimension=periodicCoords.periodicity(i);
                CoordsOfCurrDimension=periodicCoords.coords(i);
                
                while CoordsOfCurrDimension<0
                    CoordsOfCurrDimension=...
                        CoordsOfCurrDimension+PeriodicityOfCurrDimension;
                end
                
                periodicCoords.coords(i)=...
                    mod(CoordsOfCurrDimension,PeriodicityOfCurrDimension);
            end
        end
        
        function periodicCoords=translate(periodicCoords,translationVector)
            if length(translationVector) ~= periodicCoords.dimensionality
                error('Translation vector is not of same dimensions as coords!')
            end
            
            periodicCoords.coords=periodicCoords.coords+translationVector;
            
            
            periodicCoords=periodicCoords.reduceToUnitCell;
        end
        
        
        function periodicCoords=reflectAboutAxis(periodicCoords,reflectionAxis)
            IdxPossibleAxes=1:periodicCoords.dimensionality;
            if ~ismember(reflectionAxis,IdxPossibleAxes)
                error('Specified axis is not possible!')
            elseif length(reflectionAxis) > 1
                error('reflectionAxis must be a scalar')
            end
            
            periodicCoords.coords(reflectionAxis)=...
                -periodicCoords.coords(reflectionAxis);
            
            periodicCoords=periodicCoords.reduceToUnitCell;
        end
        
        function isCoordsEqual=isCoordsEqual(periodicCoords1,periodicCoords2)
            
            periodicCoords1=periodicCoords1.reduceToUnitCell;
            periodicCoords2=periodicCoords2.reduceToUnitCell;
            
            isCoordsEqual=periodicCoords1.coords == periodicCoords2.coords;
        end
    end
    
end

