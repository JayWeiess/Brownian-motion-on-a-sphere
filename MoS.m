function [ BMS ] = MoS( BM )
%������ɢ������BM, ȷ����Ӧ���ظ�ƽ��켣�˶�������켣��ɢ��BMS
%   �˴���ʾ��ϸ˵��

n=size(BM,1);
BMS=zeros(n,3);
BMS(1,:)=[-1,0,0];
SV=[0,0,1]; 
L=sqrt((BM(2,:)-BM(1,:))*(BM(2,:)-BM(1,:))');
[SN,BMS(2,:)]=SP( SV,BMS(1,:),L );
PI=zeros(1,2);
PO=zeros(1,2);
SV=zeros(1,3);
for i=2:n-1
    PI=(BM(i,:)-BM(i-1,:))/(sqrt((BM(i,:)-BM(i-1,:))*(BM(i,:)-BM(i-1,:))'));
    PO=(BM(i+1,:)-BM(i,:))/(sqrt((BM(i+1,:)-BM(i,:))*(BM(i+1,:)-BM(i,:))'));
    SV=SOV(BMS(i,:),PI,PO,SN);
    L=(sqrt((BM(i+1,:)-BM(i,:))*(BM(i+1,:)-BM(i,:))'));
    [SN,BMS(i+1,:)]=SP(SV,BMS(i,:),L);
end

end

