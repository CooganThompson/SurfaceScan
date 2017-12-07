% %Load PES
% %19_1
% clear
% clc
% close all
% PES=fillinviasymmetry(diffusionPES('19_1vasp.dat'),'c3v');
% PES.plotGoogleMapsView
% 
% %Convert to matrix
% for i=1:length(PES.alpha)
% Mat(PES.alpha(i)/3+1,PES.beta(i)/3+1)=PES.energies(i);
% end
% 
% %Get minimums
% [n,m]=size(Mat);
% Minimums=[];
% for i=2:n-1
%     for j=1:m
%         if Mat(1+mod(i-2,n),j)>=Mat(i,j)&&Mat(1+mod(i,n),j)>=Mat(i,j)&& ...
%             Mat(i,1+mod(j-2,m))>=Mat(i,j)&&Mat(i,1+mod(j,m))>=Mat(i,j)&& ...
%             Mat(1+mod(i-2,n),1+mod(j-2,m))>=Mat(i,j)&&Mat(1+mod(i,n),1+mod(j,m))>=Mat(i,j)&& ...
%             Mat(1+mod(i,n),1+mod(j-2,m))>=Mat(i,j)&&Mat(1+mod(i-2,n),1+mod(j,m))>=Mat(i,j)
%         %disp('Minimum found!')
%         Minimums=[Minimums;(i-1)*3,(j-1)*3,Mat(i,j)];
%         end 
%     end
% end
% load MinimumsForAu19.mat
% Minimums=[Minimums;[Minimums(:,1) Minimums(:,2)+120 Minimums(:,3) Minimums(:,4) Minimums(:,5)];[Minimums(:,1) Minimums(:,2)+240 Minimums(:,3) Minimums(:,4) Minimums(:,5)]];
% Minimums=[Minimums;[Minimums(:,1) -Minimums(:,2)+360 Minimums(:,3) Minimums(:,4) Minimums(:,5)]];
% hold on
% scatter3(Minimums(:,1),Minimums(:,2),Minimums(:,3)+.01,'ro')
% 
% AddingAtomDots_Au19_1
% hold on
% scatter3(Minimums(:,1),Minimums(:,2),Minimums(:,3)+.01,'ro')
% 
% %manually remove points 
% % MinimumsReduced=Minimums
% % scatter(MinimumsReduced(:,4),MinimumsReduced(:,3),100,MinimumsReduced(:,5),'filled')
% % title('Energy Relative to Best vs CN number')
% 
% 
% 
% 
% 




%Load PES
%19_1+HCOO
clear
clc
close all
PES=fillinviasymmetry(diffusionPES('HCOO+Hvasp.dat'),'c1v');
PES.plotGoogleMapsView


%Convert to matrix
for i=1:length(PES.alpha)
Mat(PES.alpha(i)/3+1,PES.beta(i)/3+1)=PES.energies(i);
end

%Get minimums
[n,m]=size(Mat);
Minimums=[];
for i=2:n-1
    for j=1:m
        if Mat(1+mod(i-2,n),j)>=Mat(i,j)&&Mat(1+mod(i,n),j)>=Mat(i,j)&& ...
            Mat(i,1+mod(j-2,m))>=Mat(i,j)&&Mat(i,1+mod(j,m))>=Mat(i,j)&& ...
            Mat(1+mod(i-2,n),1+mod(j-2,m))>=Mat(i,j)&&Mat(1+mod(i,n),1+mod(j,m))>=Mat(i,j)&& ...
            Mat(1+mod(i,n),1+mod(j-2,m))>=Mat(i,j)&&Mat(1+mod(i-2,n),1+mod(j,m))>=Mat(i,j)
        %disp('Minimum found!')
        Minimums=[Minimums;(i-1)*3,(j-1)*3,Mat(i,j)];
        end 
    end
end
load MinimumsForAu19.mat
Minimums(:,1)=180-Minimums(:,1)
Minimums=[Minimums;[Minimums(:,1) Minimums(:,2)+120 Minimums(:,3) Minimums(:,4) Minimums(:,5)];[Minimums(:,1) Minimums(:,2)+240 Minimums(:,3) Minimums(:,4) Minimums(:,5)]];
Minimums=[Minimums;[Minimums(:,1) -Minimums(:,2)+360 Minimums(:,3) Minimums(:,4) Minimums(:,5)]];
hold on
scatter3(Minimums(:,1),Minimums(:,2),Minimums(:,3)+.01,'ro')

AddingAtomDots_Au19_1
hold on
Minimums(:,1)=180-Minimums(:,1)
scatter3(Minimums(:,1),Minimums(:,2),Minimums(:,3)+.01,'ro')

figure
AddingAtomDots_HCOO_H
hold on
Minimums(:,1)=180-Minimums(:,1)
scatter3(Minimums(:,1),Minimums(:,2),Minimums(:,3)+.01,'ro')


