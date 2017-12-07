function PlotBarReactionsWithBoltzmanDisttribution(PESFull,options)
PESFull.Class=PESFull.Class.interpolatePESGrid;
KMCresults=PESFull.Results;

if nargin <= 1
    options=0;
end

Energies=zeros((length(PESFull.Mins)),1);
for i=1:length(PESFull.Mins)
    Energies(i) = getMEPEnergies(PESFull.Class,PESFull.Mins(i,1),PESFull.Mins(i,2));
end

T=298.15; %in Kelvin
kT=0.0000861731*T; %in eV

N_i=exp(-Energies/kT);
P_i=N_i/sum(N_i);

%Time on each site
count2=zeros(1,max(max(KMCresults.Locations)));
for i=2:2:length(KMCresults.Names)
    count2(KMCresults.Locations(1,i-1))=KMCresults.EventsIndividual(end,i-1)*KMCresults.ResidenceTimesIndividual(end,i-1)+count2(KMCresults.Locations(1,i-1));
end
figure
bar(count2)
hold on
P_i_normalized=P_i*sum(count2);
plot(1:length(P_i_normalized),P_i_normalized,'or')

title('Amount of time on each site')
if options=='log'
    set(gca,'YScale','log');
end
end