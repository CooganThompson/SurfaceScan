classdef (Abstract) coordinates < objectsHandler
    properties
        coords
        
    end
    
    properties (Access=protected)
        dimensionality
    end
    
    methods
        function coordinates=coordinates(coords)
            coordinates.coords=coords;
            coordinates.dimensionality=size(coords,2);
        end
        
        function coordinates=getRows(coordinates,rowNum)
            if rowNum > size(coordinates.coords,1)
                error('Row number exceeds number of rows of coordinates!')
            end
            coordinates.coords=coordinates.coords(rowNum,:);
        end
        
        function COM=getCOM(cartesianCoords)
            COMTmp=mean(cartesianCoords.coords,1);
            COM=cartesianCoords.createObjectOfSimilarClass(COMTmp);
        end
        
        function centeredCoordinates=centerAbout(coordinates,center)
            
            % //
            centerVector=...
                centerAbout_DetermineCenterCoordsFromInput(coordinates,center);
            
            % // Centering operation
            centeredCoordinatesTmp=...
                bsxfun(@minus,coordinates.coords,centerVector);
            
            %
            centeredCoordinates=coordinates.createObjectOfSimilarClass(centeredCoordinatesTmp);
            
        end
        
        
        function object_out=plus(object1,object2)
            
            n_points_object1=size(object1.coords,1);
            n_points_object2=size(object2.coords,1);
            
            if n_points_object1 >= n_points_object2
                CoordsDifferenceTmp=bsxfun(@plus,object1.coords,object2.coords);
            else
                error('Second sphericalPESCoords object has more points than first one!')
            end
            
            object_out=object1.createObjectOfSimilarClass(CoordsDifferenceTmp);
        end
        
        function object_out=minus(object1,object2)
            object2.coords=-object2.coords;
            object_out=object1+object2;
        end
    end
    
    methods (Hidden)
        function center=...
                centerAbout_DetermineCenterCoordsFromInput(coordinates,center)
            
            if strcmpi(center,'COM')
                centerTmp=coordinates.getCOM;
                center=centerTmp.coords;
            elseif strcmpi( class(center) , class(coordinates) )
                if length(center) ~= coordinates.dimensionality
                    error('Dimensions of center and coordinates to be centered do not match!')
                end
                center=center.coords;
            elseif strcmpi( class(center) , 'double' )
                center=center;
            else
                error('Center is of type that cannot be handled!')
            end
        end
    end
    
    
end