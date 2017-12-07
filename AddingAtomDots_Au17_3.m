clear;clf;clc;close all 
PES=fillinviasymmetry(diffusionPES('17_3vasp.dat'),'cs'); 
NameofPESdata='PESFull_Au17_3.mat'; 
load(['data/' NameofPESdata]); 
XYZDirectlyFromZMatConfig1.XYZDirectlyFromZMAT=[ 
      0.000000000000      0.000000000000      0.000000000000
      0.000000000000      0.000000000000      2.742230000000
      2.575215534212      0.000000000000      1.522895984092
      1.093251006913      2.350637539536      0.884618894460
     -1.690207969995      1.706906608711      1.179369224561
      3.811948755357      2.361111654204      1.428529247385
      2.290103034386     -0.909722279516      4.013262448830
      4.721449137281      4.842212398413      2.075910300914
      3.808962591032      3.273069056189      4.120595638726
      2.026285682649      4.610992817735      2.168620513728
      1.105673948250      2.719878380572      4.102301844285
      2.922128176817      1.330036532867      5.643534362072
     -0.592444752976      4.131959651226      2.324491873979
     -3.234224232027      3.471587248936      2.510471503129
     -1.618530835196      1.922991494272      4.071429357838
      0.093294356112      0.432060794038      5.471730286175
      2.013920485450     -1.014694366367      6.717396268991
      0.170683322552      4.199444795963      4.220972640327
      1.067059972005      3.812604976584      5.468668393776
]; 
XYZDirectlyFromZMatConfig1.MetalAtomsIdx=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 ];
XYZDirectlyFromZMatConfig1.AtomLinkedToCOMforZAxisAlignmentIdx=11;
XYZDirectlyFromZMatConfig1.ScanningAtomIdx=18;
XYZDirectlyFromZMatConfig1.OwnAlphaBeta=sphericalPESCoords([30 30]);
XYZDirectlyFromZMatConfig1.COM=mean(XYZDirectlyFromZMatConfig1.XYZDirectlyFromZMAT(XYZDirectlyFromZMatConfig1.MetalAtomsIdx,:),1); 

XYZDirectlyFromZMatConfig2.XYZDirectlyFromZMAT=[ 
      0.000000000000      0.000000000000      0.000000000000
      0.000000000000      0.000000000000      2.742230000000
      2.575215534212      0.000000000000      1.522895984092
      1.093251006913      2.350637539536      0.884618894460
     -1.690207969995      1.706906608711      1.179369224561
      3.811948755357      2.361111654204      1.428529247385
      2.290103034386     -0.909722279516      4.013262448830
      4.721449137281      4.842212398413      2.075910300914
      3.808962591032      3.273069056189      4.120595638726
      2.026285682649      4.610992817735      2.168620513728
      1.105673948250      2.719878380572      4.102301844285
      2.922128176817      1.330036532867      5.643534362072
     -0.592444752976      4.131959651226      2.324491873979
     -3.234224232027      3.471587248936      2.510471503129
     -1.618530835196      1.922991494272      4.071429357838
      0.093294356112      0.432060794038      5.471730286175
      2.013920485450     -1.014694366367      6.717396268991
     -0.031475624399      5.301717233549      3.542697176541
      1.067059972005      3.812604976584      5.468668393776
];
XYZDirectlyFromZMatConfig2.MetalAtomsIdx=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 ];
XYZDirectlyFromZMatConfig2.AtomLinkedToCOMforZAxisAlignmentIdx=11; 
XYZDirectlyFromZMatConfig2.ScanningAtomIdx=18;
XYZDirectlyFromZMatConfig2.OwnAlphaBeta=sphericalPESCoords([45 45]);
XYZDirectlyFromZMatConfig2.COM=mean(XYZDirectlyFromZMatConfig2.XYZDirectlyFromZMAT(XYZDirectlyFromZMatConfig2.MetalAtomsIdx,:),1); 


PES.plotContourWithAtomRegions(XYZDirectlyFromZMatConfig1,XYZDirectlyFromZMatConfig2,1)
PES.plot3DRepresentation(XYZDirectlyFromZMatConfig1,XYZDirectlyFromZMatConfig2)
PESFull.AlignmentData={XYZDirectlyFromZMatConfig1,XYZDirectlyFromZMatConfig2};
%save(['data/' NameofPESdata],'PESFull')