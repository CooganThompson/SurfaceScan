function TransformationMatrix=diffusion_determineRotation(cluster1,cluster2)
% Determine transformation that has to be applied to cluster 1 to obtain
% cluster 2.
% Multiply coords of cluster 1 by TransformationMatrix
% i.e. cluster1*TransformationMatrix

%cluster1='b1';
%cluster2='pw91';

coords1=getclustercoords(cluster1);
coords2=getclustercoords(cluster2);

% Center
coords1=centerPoints(coords1,'com');
coords2=centerPoints(coords2,'com');

% 
matchingpoints=[12,17,8];
coords2formatching=coords2(matchingpoints,:);
[a,b,TransformationMatrix,d]=determineBestLinearTransform(coords1,coords2formatching,matchingpoints,'rigidtransform');
TransformationMatrix=TransformationMatrix(1:3,1:3);

close
newcoords1=coords1*TransformationMatrix;
plot3(newcoords1(:,1),newcoords1(:,2),newcoords1(:,3),'linestyle','none','marker','o')
hold on
plot3(coords2(:,1),coords2(:,2),coords2(:,3),'linestyle','none','marker','o','markersize',18)


end
function coords=getclustercoords(name)
switch lower(name)
    case 'b1'
        coords=[
   -2.8878   -1.2495   -0.6421
   -1.8838    0.1619    1.5456
    0.1071    1.3838    3.2134
    1.6190    2.6338    1.3090
   -1.1620    2.8730    1.3040
    0.1564    1.8133   -0.9599
    2.8894    1.2509   -0.6343
    4.0307   -0.3434   -2.5149
    1.3285   -0.1150   -2.6283
    2.6375   -1.7210   -0.6311
    1.8779   -0.1589    1.5495
   -0.1207   -1.3875    3.2182
    1.1533   -2.8718    1.3110
   -0.1535   -1.8089   -0.9543
   -1.3125    0.1113   -2.6364
   -2.6300    1.7191   -0.6485
   -4.0157    0.3414   -2.5368
   -1.6291   -2.6337    1.3093];
    case 'pw91'
        coords=[         
  1.9204	-1.8134	-2.0318
-0.3618	0.0382	-2.0559
-2.7806	1.7008	-1.8065
-1.6943	2.9171	0.5378
-0.1384	2.7725	-1.9976
0.9424	1.6745	0.4781
-0.7117	1.6964	2.7795
0	0	4.824
1.6237	-0.1568	2.5818
-1.0182	-1.5359	2.7821
-1.9846	0.1896	0.6024
-3.0586	-1.1335	-1.7985
-2.2212	-2.5438	0.5425
0.6055	-1.8254	0.4789
3.0255	-0.2919	0.2904
2.2365	1.417	-2.029
4.2837	-0.414	-2.1775
-0.6707	-2.6873	-1.9954
];
end
end