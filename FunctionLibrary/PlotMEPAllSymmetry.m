function PlotMEPAllSymmetry(MEPObject,Symmetry)

MEPObjectSSS=GiveSymmetryVersions(MEPObject,Symmetry);
for k=1:length(MEPObjectSSS)
    MEPObjectSSS{k}.plotMEPLines
end

end