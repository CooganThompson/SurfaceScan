clear;clf;clc;close all 
%SAMPLE PES=fillinviasymmetry(diffusionPES('17_1vasp.dat'),'c3v'); 
%SAMPLE NameofPESdata='PESFull_Au18_1.mat'; 
load(NameofPESdata); 
XYZDirectlyFromZMatConfig1.XYZDirectlyFromZMAT=[ 
      0.000000000000      0.000000000000      0.000000000000
      0.000000000000      0.000000000000      2.713320000000
      2.473664473316      0.000000000000      1.117073680081
      0.720816974449      2.363121523825      1.121241635497
      2.537912005852     -0.132091984157      3.932184101695
      0.124899372913      0.104240242692      5.372215876841
      0.388834167318      0.290896711377      8.063711983582
      0.999153031781      2.575704113178      6.727178799572
      0.580227633699      2.478501613866      3.933019887276
      1.651778337121      4.709841936792      5.275204727817
      2.382689643184      6.786431083123      3.687345669333
      4.133165261064      4.711480970123      3.691566727467
      1.507124759793      4.603788504270      2.334248489228
      3.501164325448      2.588566270100      5.431507720775
      5.714649568348      2.565145842947      3.686197252103
      3.326475909181      2.462107453761      2.194711317134
      4.844124406903      0.090785962617      2.331431184176
      7.184385868713      0.285068284326      3.691255607736
      4.989765599184      0.200784913522      5.283107531222
      2.746926312939      0.197558810711      6.722697590349
      0.416716243929      5.928469203201      4.992746803144
      2.344532603659      8.534884270074      3.624453084271
]; 
XYZDirectlyFromZMatConfig1.MetalAtomsIdx=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
XYZDirectlyFromZMatConfig1.AtomLinkedToCOMforZAxisAlignmentIdx=11;
XYZDirectlyFromZMatConfig1.ScanningAtomIdx=21;
XYZDirectlyFromZMatConfig1.OwnAlphaBeta=sphericalPESCoords([30 30]);
XYZDirectlyFromZMatConfig1.COM=mean(XYZDirectlyFromZMatConfig1.XYZDirectlyFromZMAT(XYZDirectlyFromZMatConfig1.MetalAtomsIdx,:),1); 

XYZDirectlyFromZMatConfig2.XYZDirectlyFromZMAT=[ 
      0.000000000000      0.000000000000      0.000000000000
      0.000000000000      0.000000000000      2.713320000000
      2.473664473316      0.000000000000      1.117073680081
      0.720816974449      2.363121523825      1.121241635497
      2.537912005852     -0.132091984157      3.932184101695
      0.124899372913      0.104240242692      5.372215876841
      0.388834167318      0.290896711377      8.063711983582
      0.999153031781      2.575704113178      6.727178799572
      0.580227633699      2.478501613866      3.933019887276
      1.651778337121      4.709841936792      5.275204727817
      2.382689643184      6.786431083123      3.687345669333
      4.133165261064      4.711480970123      3.691566727467
      1.507124759793      4.603788504270      2.334248489228
      3.501164325448      2.588566270100      5.431507720775
      5.714649568348      2.565145842947      3.686197252103
      3.326475909181      2.462107453761      2.194711317134
      4.844124406903      0.090785962617      2.331431184176
      7.184385868713      0.285068284326      3.691255607736
      4.989765599184      0.200784913522      5.283107531222
      2.746926312939      0.197558810711      6.722697590349
      0.092900528983      4.278584394411      4.525655474238
      2.344532603659      8.534884270074      3.624453084271
];
XYZDirectlyFromZMatConfig2.MetalAtomsIdx=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
XYZDirectlyFromZMatConfig2.AtomLinkedToCOMforZAxisAlignmentIdx=11; 
XYZDirectlyFromZMatConfig2.ScanningAtomIdx=21;
XYZDirectlyFromZMatConfig2.OwnAlphaBeta=sphericalPESCoords([45 45]);
XYZDirectlyFromZMatConfig2.COM=mean(XYZDirectlyFromZMatConfig2.XYZDirectlyFromZMAT(XYZDirectlyFromZMatConfig2.MetalAtomsIdx,:),1); 


PES.plotContourWithAtomRegions(XYZDirectlyFromZMatConfig1,XYZDirectlyFromZMatConfig2,1)
PES.plot3DRepresentation(XYZDirectlyFromZMatConfig1,XYZDirectlyFromZMatConfig2)
PESFull.AlignmentData={XYZDirectlyFromZMatConfig1,XYZDirectlyFromZMatConfig2};
save(NameofPESdata,'PESFull')