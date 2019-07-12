function [ SN,SNP,PreQ ] = SP( SV,BMS,L )
%% SP ������ʻ������ȷ����һ����������
% SN--��һ������켣��Բ�����ڵ�ƽ��ķ���
% SNP-��һ�������Ӧλ��
% SV--��һ�����淽��
% BMS-��ǰ����λ��
% L---��һ������
% PreQ-�洢һԪ���η��̵�ϵ��(p,q,l,Delta)
%%

%   �˴���ʾ��ϸ˵��
SN=zeros(1,3);
%�򻯵�ǰ����ļǺ�
x0=BMS(1,1);
y0=BMS(1,2);
z0=BMS(1,3);
%������һ������켣��Բ�����ڵ�ƽ��ķ���
SN(1,1)=y0*SV(1,3)-z0*SV(1,2);
SN(1,2)=z0*SV(1,1)-x0*SV(1,3);
SN(1,3)=x0*SV(1,2)-y0*SV(1,1);
a=SN(1,1);
b=SN(1,2);
c=SN(1,3);
p=1;
PreQ(1,6)=0;

if abs(SV(1,3))>0.57
%ȷ����Ԫһ�η��̵�ϵ��
q=2*cos(L)*(z0);
l=(a^2+b^2)*cos(L)^2-SV(1,3)^2;
u=q^2-4*p*l;
    if u>=0
    PreQ=[p,q,l,u,SV(1,3),1];
 %���������ʽ�ⷽ��
    z1=(q+sqrt(q^2-4*p*l))/(2*p);
    x1=(SV(1,1)*z1-b*cos(L))/SV(1,3);
    y1=(SV(1,2)*z1+a*cos(L))/SV(1,3);

    z2=(q-sqrt(q^2-4*p*l))/(2*p);
    x2=(SV(1,1)*z2-b*cos(L))/SV(1,3);
    y2=(SV(1,2)*z2+a*cos(L))/SV(1,3);


%-----------------------------------------------------------
    else
        if abs(SV(1,2))>0.57
 
    q=2*cos(L)*(y0);
    l=(a^2+c^2)*cos(L)^2-SV(1,2)^2;
    u=q^2-4*p*l;
    if u>=0
    PreQ=[p,q,l,q^2-4*p*l,SV(1,2),2];
%���������ʽ�ⷽ��
    y1=(q+sqrt(q^2-4*p*l))/(2*p);
    x1=(SV(1,1)*y1+c*cos(L))/SV(1,2);
    z1=(SV(1,3)*y1-a*cos(L))/SV(1,2);

    y2=(q-sqrt(q^2-4*p*l))/(2*p);
    x2=(SV(1,1)*y2+c*cos(L))/SV(1,2);
    z2=(SV(1,3)*y2-a*cos(L))/SV(1,2);
    
    
 %-----------------------------------------------------------
    else 
            q=2*cos(L)*(x0);
            l=(b^2+c^2)*cos(L)^2-SV(1,1)^2;
            PreQ=[p,q,l,q^2-4*p*l,SV(1,1),3];
%���������ʽ�ⷽ��
            x1=(q+sqrt(q^2-4*p*l))/(2*p);
            y1=(SV(1,2)*x1-c*cos(L))/SV(1,1);
            z1=(SV(1,3)*x1+b*cos(L))/SV(1,1);

            x2=(q-sqrt(q^2-4*p*l))/(2*p);
            y2=(SV(1,2)*x1-c*cos(L))/SV(1,1);
            z2=(SV(1,3)*x1+b*cos(L))/SV(1,1);
            
            
        
    end
        else
            q=2*cos(L)*(x0);
            l=(b^2+c^2)*cos(L)^2-SV(1,1)^2;
            PreQ=[p,q,l,q^2-4*p*l,SV(1,1),3];
%���������ʽ�ⷽ��
            x1=(q+sqrt(q^2-4*p*l))/(2*p);
            y1=(SV(1,2)*x1-c*cos(L))/SV(1,1);
            z1=(SV(1,3)*x1+b*cos(L))/SV(1,1);

            x2=(q-sqrt(q^2-4*p*l))/(2*p);
            y2=(SV(1,2)*x1-c*cos(L))/SV(1,1);
            z2=(SV(1,3)*x1+b*cos(L))/SV(1,1);
        end
    end
    else
        if abs(SV(1,2))>0.57
 
    q=2*cos(L)*(y0);
    l=(a^2+c^2)*cos(L)^2-SV(1,2)^2;
    u=q^2-4*p*l;
    if u>=0
    PreQ=[p,q,l,q^2-4*p*l,SV(1,2),2];
%���������ʽ�ⷽ��
    y1=(q+sqrt(q^2-4*p*l))/(2*p);
    x1=(SV(1,1)*y1+c*cos(L))/SV(1,2);
    z1=(SV(1,3)*y1-a*cos(L))/SV(1,2);

    y2=(q-sqrt(q^2-4*p*l))/(2*p);
    x2=(SV(1,1)*y2+c*cos(L))/SV(1,2);
    z2=(SV(1,3)*y2-a*cos(L))/SV(1,2);
    
    
 %-----------------------------------------------------------
    else 
            q=2*cos(L)*(x0);
            l=(b^2+c^2)*cos(L)^2-SV(1,1)^2;
            PreQ=[p,q,l,q^2-4*p*l,SV(1,1),3];
%���������ʽ�ⷽ��
            x1=(q+sqrt(q^2-4*p*l))/(2*p);
            y1=(SV(1,2)*x1-c*cos(L))/SV(1,1);
            z1=(SV(1,3)*x1+b*cos(L))/SV(1,1);

            x2=(q-sqrt(q^2-4*p*l))/(2*p);
            y2=(SV(1,2)*x1-c*cos(L))/SV(1,1);
            z2=(SV(1,3)*x1+b*cos(L))/SV(1,1);
            
            
        
    end
        else
            q=2*cos(L)*(x0);
            l=(b^2+c^2)*cos(L)^2-SV(1,1)^2;
            PreQ=[p,q,l,q^2-4*p*l,SV(1,1),3];
%���������ʽ�ⷽ��
            x1=(q+sqrt(q^2-4*p*l))/(2*p);
            y1=(SV(1,2)*x1-c*cos(L))/SV(1,1);
            z1=(SV(1,3)*x1+b*cos(L))/SV(1,1);

            x2=(q-sqrt(q^2-4*p*l))/(2*p);
            y2=(SV(1,2)*x1-c*cos(L))/SV(1,1);
            z2=(SV(1,3)*x1+b*cos(L))/SV(1,1);
        end
    end
        
if SV*[x1;y1;z1]>0
    SNP=[x1,y1,z1];
else
    SNP=[x2,y2,z2];
end

