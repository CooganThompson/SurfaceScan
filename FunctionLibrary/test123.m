
[Length,~]=size(MatForm);
for i=1:Length
    MatForm(Length+i,:)=MatForm(i,:);
    MatForm(Length+i,2)=360-MatForm(Length+i,2);
    MatForm(Length+i,4)=-MatForm(Length+i,4);
end
MatForm=RemoveDuplicates(MatForm);

