clear;clf;clc 
PES=fillinviasymmetry(diffusionPES('HCOO+H_2.dat'),'c1v'); 
XYZDirectlyFromZMatConfig1.XYZDirectlyFromZMAT=[ 
      0.000000000000      0.000000000000      0.000000000000
      0.000000000000      0.000000000000      2.724320000000
      0.168067400701      0.000000000000      5.398654141599
      0.526958280249     -0.072914126531      8.088841099178
      2.227435553026      1.126883205818      1.101196883547
      2.234595249414      1.243026206969      3.920376416396
      2.585906642241      1.058267320845      6.711532703880
      4.393704606120      2.126778637292      2.320444479453
      4.584559812598      2.060087279268      5.249405305040
      6.591868896771      2.978547861578      3.663278513267
      1.707160458997     -1.739602255756      1.128870474130
      1.679914141308     -1.902549545170      3.971603827779
      2.090406483139     -1.827406069654      6.714249200178
      4.119242318166     -0.669348943172      2.195330416472
      4.245942211862     -0.787865502628      5.273685473913
      6.265843032678      0.277923068375      3.672918263363
      3.427203401264     -3.350657049085      2.372547809128
      3.588155634446     -3.545501374732      5.286217955095
      5.890766455691     -2.389005456755      3.613960688452
      4.704237836248     -5.380860260908      4.840795441794
      6.436178651247     -4.514051481805      3.591235153012
      5.775869584763     -5.432472330672      4.165862222236
      6.201117915614     -6.449306399489      4.046951421631
      7.010099980824     -8.000613853366      4.085593372858
]; 
XYZDirectlyFromZMatConfig1.MetalAtomsIdx=1:19;
XYZDirectlyFromZMatConfig1.AtomLinkedToCOMIdx=24;
XYZDirectlyFromZMatConfig1.ScanningAtomIdx=24;
XYZDirectlyFromZMatConfig1.OwnAlphaBeta=sphericalPESCoords([0 0]);
XYZDirectlyFromZMatConfig2.XYZDirectlyFromZMAT=[
      0.000000000000      0.000000000000      0.000000000000
      0.000000000000      0.000000000000      2.724320000000
      0.168067400701      0.000000000000      5.398654141599
      0.526958280249     -0.072914126531      8.088841099178
      2.227435553026      1.126883205818      1.101196883547
      2.234595249414      1.243026206969      3.920376416396
      2.585906642241      1.058267320845      6.711532703880
      4.393704606120      2.126778637292      2.320444479453
      4.584559812598      2.060087279268      5.249405305040
      6.591868896771      2.978547861578      3.663278513267
      1.707160458997     -1.739602255756      1.128870474130
      1.679914141308     -1.902549545170      3.971603827779
      2.090406483139     -1.827406069654      6.714249200178
      4.119242318166     -0.669348943172      2.195330416472
      4.245942211862     -0.787865502628      5.273685473913
      6.265843032678      0.277923068375      3.672918263363
      3.427203401264     -3.350657049085      2.372547809128
      3.588155634446     -3.545501374732      5.286217955095
      5.890766455691     -2.389005456755      3.613960688452
      4.704237836248     -5.380860260908      4.840795441794
      6.436178651247     -4.514051481805      3.591235153012
      5.775869584763     -5.432472330672      4.165862222236
      6.201117915614     -6.449306399489      4.046951421631
      2.488075894915     -4.399804861375      3.475546297604
    ];
XYZDirectlyFromZMatConfig2.MetalAtomsIdx=1:19;
XYZDirectlyFromZMatConfig2.AtomLinkedToCOMIdx=24;
XYZDirectlyFromZMatConfig2.ScanningAtomIdx=24;
XYZDirectlyFromZMatConfig2.OwnAlphaBeta=sphericalPESCoords([45 45]);
XYZCoordsa0b0=17*[ 
0.56368289882 0.49592502824 0.58053420588
0.56336635412 0.61091353235 0.31514951294
0.56039972941 0.40356272118 0.44084209588
0.54450208588 0.32632575471 0.29430860000
0.47978584118 0.56545628000 0.45447611882
0.47090902353 0.47507743353 0.29211826059
0.46014783882 0.36001175941 0.56206734412
0.45194079882 0.27650914176 0.42966596412
0.44424270059 0.20104530353 0.29106144471
0.40699281000 0.62639156882 0.31494055706
0.38380585294 0.51308761647 0.58039840941
0.37015704000 0.34310202412 0.29460902882
0.36952441235 0.42188143941 0.44112351353
0.31581030765 0.57559689882 0.45193969000
0.30524985588 0.48651600471 0.30266197765
0.24706980706 0.63426525059 0.31972688706
0.72154871059 0.58803998529 0.31977907353
0.64242674765 0.54412021706 0.45174677824
0.63571058118 0.45459201176 0.30266774706
0.47155149059 0.48059692882 0.73682664941
0.47058826471 0.47059476824 0.80125779941
0.47058823529 0.47058823529 0.90419897529
0.53955829706 0.47869860941 0.70614488000
0.40455796118 0.49121913294 0.70587213412
    ];
PES.plotContourWithAtomRegions(XYZDirectlyFromZMatConfig1,...
    XYZDirectlyFromZMatConfig2,XYZCoordsa0b0)