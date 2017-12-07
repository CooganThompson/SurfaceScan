function Speed=CalculateMobility(PES,MinsLoc,FileName)
%Calculates the net diffusion speed

fid=fopen('data/18_1_Diffusion_KMC_Trial.txt','r');
a=textscan(fid,'%s %f');

Names=a{1};
Times=a{2};
Names(1)=[];
Times(1)=[];
for i=2:2:length(Names)
    Name=char(Names(i));
    n=sscanf(Name,'diff/H/%d_TO_%d');
    Data(i/2).Name=Name;
    Data(i/2).Points=n;
    Data(i/2).Times=Times(i)+Times(i-1);
        
        X1=sind(MinsLoc(i,1)).*cosd(MinsLoc(i,2));
        X2=sind(MinsLoc(j,1)).*cosd(MinsLoc(j,2));
        Y1=sind(MinsLoc(i,1)).*sind(MinsLoc(i,2));
        Y2=sind(MinsLoc(j,1)).*sind(MinsLoc(j,2));
        Z1=cosd(MinsLoc(i,1));
        Z2=cosd(MinsLoc(j,1));
        dist=(X1-X2)^2+(Y1-Y2)^2+(Z1-Z2).^2;
    
    
    
    
    Data(i/2).LengthBetween=
    
    
    
    
end



CumDistance=0;
for i=1:length(Data.Name)

        X1=sind(MinsLoc(i,1)).*cosd(MinsLoc(i,2));
        X2=sind(MinsLoc(j,1)).*cosd(MinsLoc(j,2));
        Y1=sind(MinsLoc(i,1)).*sind(MinsLoc(i,2));
        Y2=sind(MinsLoc(j,1)).*sind(MinsLoc(j,2));
        Z1=cosd(MinsLoc(i,1));
        Z2=cosd(MinsLoc(j,1));
        dist=(X1-X2)^2+(Y1-Y2)^2+(Z1-Z2).^2;





end