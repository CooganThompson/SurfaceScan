function NewPoints=alignWithZAxis(Coords,iZAxisAtom,COM)
% // Aligns the COM and iZAxisAtom axis with the z-axis, and returns
% // transformated points
% // COM must be specified seperately because we coords might
% // contain other atoms that are not included in COM


% Debug
%{
iZAxisAtom=19;
Coords=cartesianCoords(...
    [  0.000000000000      0.000000000000      0.000000000000
     0.000000000000      0.000000000000      2.899320000000
      2.851434648216      0.000000000000      2.917087021442
      1.321396967650      1.996753123757      1.346197452945
        3.218814952978      0.254268897523      0.074451427246
       4.490899283932     -1.754950498045      1.600213554517
       5.788742041375     -3.755656049369      0.234382436752
        4.478802084000     -1.893411602728     -1.341817239027
        3.062336128796     -4.056713898457      0.069862659513
      1.587453924346     -2.047186654130      1.589482327128
     -1.309509366175     -2.111867066586      1.497651628454
      0.327124844118     -4.066052758281      0.217862676449
        1.650206712212     -2.266861732447     -1.599052226286
       3.086670231629     -0.079598156927     -2.766612504038
      1.573734568035      2.032375876605     -1.380134733604
       1.539887619161      1.776547299781     -4.117005593442
        0.148720558635     -0.259940834575     -2.903017735005
       -1.203774464556     -2.066691167309     -1.351708772543
         1.444553055159      2.756685642736     -5.502732889305]);
%}





% // Get the COM-Atom vector, and normalize it
Coords = cartesianCoords(Coords);
CenteredCoords = Coords.centerAbout(COM);
COMAtom_Vector = CenteredCoords.getRows(iZAxisAtom);
COMAtom_Vector = COMAtom_Vector.normalizeRows();

% // Z axis vector
% Can be any other vector which we wish to align the COM_Atom vector with
ZAxis=[0 0 1].';

% // Get the rotation matrix that aligns the COM-Atom vector with the z axis
% vector

GG = @(A,B) [ dot(A,B) -norm(cross(A,B)) 0;
    norm(cross(A,B)) dot(A,B)  0;
    0              0           1];
FFi = @(A,B) [ A (B-dot(A,B)*A)/norm(B-dot(A,B)*A) cross(B,A) ];
UU = @(Fi,G) Fi*G/Fi;


U = UU(FFi(COMAtom_Vector.coords(:),ZAxis(:)), ...
    GG(COMAtom_Vector.coords(:),ZAxis(:)));

NewPoints=(U*CenteredCoords.coords.').';
end