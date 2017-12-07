function iElement=getClosestElement(XX,YY,testx,testy)
% // Gets indices of element closest in value to testx and testy
% Assumes XX and YY are meshgrids

% Check orientation of XX and YY
if XX(1,1)==XX(2,1) && YY(1,1)==YY(1,2)
    CaseNo=1; %XX changes in col, YY changes in row
    
elseif XX(1,1)==XX(1,2) && YY(1,1)==YY(2,1)
    CaseNo=2;%XX changes in row, YY changes in col
else
    error('hi')
end


% Just need to study one row and one col from XX and YY
switch CaseNo
    case 1
        XXNew=XX(1,:); %row
        YYNew=YY(:,1); %col
    case 2
        YYNew=YY(1,:); %row
        XXNew=XX(:,1); %col
end

XXDiff=abs(XXNew-testx);
[minXXDiff,XXIdx]=min(XXDiff);
YYDiff=abs(YYNew-testy);
[minYYDiff,YYIdx]=min(YYDiff);

% Get linear index
switch CaseNo
    case 1
        iElement=sub2ind(size(XX),YYIdx,XXIdx);
    case 2
        iElement=sub2ind(size(XX),XXIdx,YYIdx);

end
