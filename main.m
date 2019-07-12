clear all
clc
BM=BrM(1000,0.0001);
n=size(BM,1);
PreQ1=zeros(n-1,6);
PreQ2=zeros(n-1,6);
BMS=zeros(n,3);
BMS(1,:)=[0,0,1];
SV=[0,1,0]; 
SN=zeros(1,3);
L=sqrt((BM(2,:)-BM(1,:))*(BM(2,:)-BM(1,:))');


[SN,BMS(2,:)]=SP( SV,BMS(1,:),L );

 for i=2:n-1
    PI=(BM(i,:)-BM(i-1,:))/(sqrt((BM(i,:)-BM(i-1,:))*(BM(i,:)-BM(i-1,:))'));
    PO=(BM(i+1,:)-BM(i,:))/(sqrt((BM(i+1,:)-BM(i,:))*(BM(i+1,:)-BM(i,:))'));
    [SV,PreQ1(i-1,:)]=SOV(BMS(i,:),PI,PO,SN);
    L=(sqrt((BM(i+1,:)-BM(i,:))*(BM(i+1,:)-BM(i,:))'));
    [SN,BMS(i+1,:),PreQ2(i-1,:)]=SP(SV,BMS(i,:),L);
 end

plot3(BMS(1:100,1),BMS(1:100,2),BMS(1:100,3));
for i=1:1000
SR1(i,:)=[0,cos(i*2*pi/1000),sin(i*2*pi/1000)];
SR2(i,:)=[cos(i*2*pi/1000),0,sin(i*2*pi/1000)];
end
%plot3(BMS(1:485,1),BMS(1:485,2),BMS(1:485,3));
% hold on
% plot3(SR1(:,1),SR1(:,2),SR1(:,3));
% hold on
% plot3(SR2(:,1),SR2(:,2),SR2(:,3));