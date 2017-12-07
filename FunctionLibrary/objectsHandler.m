classdef (Abstract) objectsHandler
    
    methods
        function newObject=...
                createObjectOfSimilarClass(originalObject,varargin)
         
            argumentsToConstructor=varargin;
            classConstructor=str2func(class(originalObject));
            newObject=classConstructor(argumentsToConstructor{:});
        end
    end
    
end