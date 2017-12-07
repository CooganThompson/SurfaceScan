function newObject=createObjectOfSimilarClassAs(originalObject,varargin)

argumentsToConstructor=varargin;
classConstructor=str2func(class(originalObject));
newObject=classConstructor(argumentsToConstructor{:});
end