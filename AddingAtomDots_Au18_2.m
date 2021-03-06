clear;clf;clc;close all 
%SAMPLE PES=fillinviasymmetry(diffusionPES('17_1vasp.dat'),'c3v'); 
%SAMPLE NameofPESdata='PESFull_Au18_1.mat'; 
load(NameofPESdata); 
XYZDirectlyFromZMatConfig1.XYZDirectlyFromZMAT=[ 
      0.000000000000      0.000000000000      0.000000000000
      0.000000000000      0.000000000000      2.773160000000
      2.474514283724      0.000000000000      1.202300782964
      1.811327294578      1.683778413573     -1.051228457569
      2.095286751597      1.396326101452      3.756225627457
      0.456643558419      2.395662737201      1.534281849860
      2.287066759678      3.981673769448      2.777748654971
      4.215121320020      2.919456715019      4.482140194766
      2.943854647029      6.362198661615      1.444035567047
      4.795652894270      5.167446091663      3.052327196540
      6.581814408538      1.543026767565      4.271448694903
      6.136357030295      2.907703775992      1.984167647603
      4.613602287144      0.581752153852      2.676669762534
      4.297837274330      1.813133704926      0.211789698476
      3.629219267442      3.459630039792     -1.970256918130
      4.899666194980      4.540727513382      0.270648677598
      1.615825073363      4.294877644682      0.076436171167
      3.201569475697      6.101684806231     -1.280498002076
      4.574163941657      0.409207476305      4.431898244302
      7.894782519429      1.093285496157      5.337444718405
]; 
XYZDirectlyFromZMatConfig1.MetalAtomsIdx=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 ];
XYZDirectlyFromZMatConfig1.AtomLinkedToCOMforZAxisAlignmentIdx=11;
XYZDirectlyFromZMatConfig1.ScanningAtomIdx=19;
XYZDirectlyFromZMatConfig1.OwnAlphaBeta=sphericalPESCoords([30 30]);
XYZDirectlyFromZMatConfig1.COM=mean(XYZDirectlyFromZMatConfig1.XYZDirectlyFromZMAT(XYZDirectlyFromZMatConfig1.MetalAtomsIdx,:),1); 

XYZDirectlyFromZMatConfig2.XYZDirectlyFromZMAT=[ 
      0.000000000000      0.000000000000      0.000000000000
      0.000000000000      0.000000000000      2.773160000000
      2.474514283724      0.000000000000      1.202300782964
      1.811327294578      1.683778413573     -1.051228457569
      2.095286751597      1.396326101452      3.756225627457
      0.456643558419      2.395662737201      1.534281849860
      2.287066759678      3.981673769448      2.777748654971
      4.215121320020      2.919456715019      4.482140194766
      2.943854647029      6.362198661615      1.444035567047
      4.795652894270      5.167446091663      3.052327196540
      6.581814408538      1.543026767565      4.271448694903
      6.136357030295      2.907703775992      1.984167647603
      4.613602287144      0.581752153852      2.676669762534
      4.297837274330      1.813133704926      0.211789698476
      3.629219267442      3.459630039792     -1.970256918130
      4.899666194980      4.540727513382      0.270648677598
      1.615825073363      4.294877644682      0.076436171167
      3.201569475697      6.101684806231     -1.280498002076
      3.413322982088      0.634475014851      4.710883814961
      7.894782519429      1.093285496157      5.337444718405
];
XYZDirectlyFromZMatConfig2.MetalAtomsIdx=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 ];
XYZDirectlyFromZMatConfig2.AtomLinkedToCOMforZAxisAlignmentIdx=11; 
XYZDirectlyFromZMatConfig2.ScanningAtomIdx=19;
XYZDirectlyFromZMatConfig2.OwnAlphaBeta=sphericalPESCoords([45 45]);
XYZDirectlyFromZMatConfig2.COM=mean(XYZDirectlyFromZMatConfig2.XYZDirectlyFromZMAT(XYZDirectlyFromZMatConfig2.MetalAtomsIdx,:),1); 


PES.plotContourWithAtomRegions(XYZDirectlyFromZMatConfig1,XYZDirectlyFromZMatConfig2,1)
PES.plot3DRepresentation(XYZDirectlyFromZMatConfig1,XYZDirectlyFromZMatConfig2)
PESFull.AlignmentData={XYZDirectlyFromZMatConfig1,XYZDirectlyFromZMatConfig2};
save(NameofPESdata,'PESFull')