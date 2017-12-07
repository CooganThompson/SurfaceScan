function Energies=getMEPEnergies(self,XCoords,YCoords)
%%Copied out of BC's class defintion

iElement=zeros(1,length(XCoords));
for z=1:length(XCoords)
    CurrX=XCoords(z);
    CurrY=YCoords(z);
    iElement(z)=getClosestElement(self.alphaGrid,self.betaGrid,CurrX,CurrY);
end
Energies=self.energiesGrid(iElement);
end