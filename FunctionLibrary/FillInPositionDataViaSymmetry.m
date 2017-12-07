function PositionData=FillInPositionDataViaSymmetry(PositionData,symmetry)
%Uses the point group symmetry operations (kinda) so you can graph the whole surface with only each unique part specified
disp('Be sure you defined alpha and beta well')
disp('Currently only good for C_n and C_nv symmetries, C_nh is untested')

MatForm=cell2mat(PositionData');

if strcmpi(symmetry,'Cs')
    symmetry='C1v';
end

if symmetry(1)=='C'||symmetry(1)=='c'
    n=str2double(symmetry(2));
    for i=2:n
        MatForm=translateByBeta(MatForm,360/n);
    end
end

if symmetry(3)=='v'
    MatForm=reflectAboutBetaZero(MatForm);
end

if symmetry(3)=='h'
    MatForm=reflectAboutAlphaNinety(MatForm);
end

%             if (symmetry(1)=='T'||symmetry(1)=='t')&&symmetry(2)=='d'
%             end
[Length,~]=size(MatForm);
for k=1:Length
    
    PositionData{k}=MatForm(k,:);
    
end

end

function MatForm=translateByBeta(MatForm,beta)

[Length,~]=size(MatForm);
for i=1:Length
    MatForm(Length+i,:)=MatForm(i,:);
    MatForm(Length+i,2)=mod(MatForm(Length+i,2)+beta,360);
    MatForm(Length+i,3)=MatForm(i,3)*cosd(beta)-MatForm(i,4)*sind(beta);
    MatForm(Length+i,4)=MatForm(i,4)*cosd(beta)+MatForm(i,3)*sind(beta);
end
MatForm=RemoveDuplicates(MatForm);

end

function MatForm=reflectAboutBetaZero(MatForm)

[Length,~]=size(MatForm);
for i=1:Length
    MatForm(Length+i,:)=MatForm(i,:);
    MatForm(Length+i,2)=360-MatForm(Length+i,2);
    MatForm(Length+i,4)=-MatForm(Length+i,4);
end
MatForm=RemoveDuplicates(MatForm);

end

function MatForm=reflectAboutAlphaNinety(MatForm)

[Length,~]=size(MatForm);
for i=1:Length
    MatForm(Length+i,:)=MatForm(i,:);
    MatForm(Length+i,1)=180-MatForm(Length+i,1);
    MatForm(Length+i,5)=-MatForm(Length+i,5);
end
MatForm=RemoveDuplicates(MatForm);

end

function MatForm=RemoveDuplicates(MatForm)
[Length,~]=size(MatForm);
i=1;
while i<Length
    testX=MatForm(i,1)==MatForm([i+1:end],1);
    testY=MatForm(i,2)==MatForm([i+1:end],2);
    testboth=testX.*testY;
    if sum(testboth)>0
        for j=find(testboth==1)
            MatForm(i+j,:)=[];
        end
    else
        i=i+1;
        [Length,~]=size(MatForm);
    end
end
end