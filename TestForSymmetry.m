a=diffusionPES('17_1.dat');
b=translateByBeta(a,120);
c=translateByBeta(a,240);
a=a+b+c;

d=reflectAboutBetaZero(a);
a=a+d

plotGoogleMapsView(a)