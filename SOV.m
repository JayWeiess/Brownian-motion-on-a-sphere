function [ SV,PreQ,ABC ] = SOV( BMS,PI,PO,SN)
%% SOV 由平面驶出方向确定球面驶出方向
% BMS--当前球面位置
% PI---平面驶入方向
% PO---平面驶出方向
% SN---球面驶入法向
%% 通过解三元二次方程来得到相应的结果
% 简化记号
a=SN(1,1); b=SN(1,2); c=SN(1,3);
x0=BMS(1,1);
y0=BMS(1,2);
z0=BMS(1,3);

%外法向单位化
%d=sqrt(a^2+b^2+c^2);
SIN=PI(1,2)*PO(1,1)-PI(1,1)*PO(1,2);
s=y0*c-z0*b;
m=z0*a-x0*c;
l=x0*b-y0*a;
ABC=[a,b,c];
% 确定一元二次方程系数
p=1;
PreQ(1,6)=0;
% 由求根公式求解方程
if abs(l)>0.57

q=2*SIN*(c);
k=(x0^2+y0^2)*SIN^2-l^2;
PreQ=[p,q,k,q^2-4*p*k,l,1];

z1=(q-sqrt(q^2-4*p*k))/(2*p);
x1=(s*z1-y0*SIN)/l;
y1=(m*z1+x0*SIN)/l;

z2=(q+sqrt(q^2-4*p*k))/(2*p);
x2=(s*z2-y0*SIN)/l;
y2=(m*z2+x0*SIN)/l;

Oa=y0*z1-z0*y1;
Ob=z0*x1-x0*z1;
Oc=x0*y1-y0*x1;


%-----------------------------------------------------------
else if abs(m)>0.57
        q=2*SIN*b;
        k=(x0^2+z0^2)*SIN^2-m^2;
        PreQ=[p,q,k,q^2-4*p*k,m,2];
        
        y1=(q-sqrt(q^2-4*p*k))/(2*p);
        x1=(s*y1+z0*SIN)/m;
        z1=(l*y1-x0*SIN)/m;

        y2=(q+sqrt(q^2-4*p*k))/(2*p);
        x2=(s*y2+z0*SIN)/m;
        z2=(l*y2-x0*SIN)/m;
        
         Oa=y0*z1-z0*y1;
         Ob=z0*x1-x0*z1;
         Oc=x0*y1-y0*x1;
         
    else 
              q=2*SIN*a;
               k=(y0^2+z0^2)*SIN^2-s^2;
               PreQ=[p,q,k,q^2-4*p*k,s,3];
               
               x1=(q-sqrt(q^2-4*p*k))/(2*p);
               y1=(m*x1-z0*SIN)/s;
               z1=(l*x1+y0*SIN)/s;
            
               x2=(q+sqrt(q^2-4*p*k))/(2*p);
               y2=(m*x2-z0*SIN)/s;
               z2=(l*x2+y0*SIN)/s;     
               
               Oa=y0*z1-z0*y1;
               Ob=z0*x1-x0*z1;
               Oc=x0*y1-y0*x1;
               
        
    end


end
% 通过求球面驶出法向来确定解
COS=PI*PO';

if sign(SN*[Oa;Ob;Oc])==sign(COS)
    SV=[x1,y1,z1];
else
    SV=[x2,y2,z2];
end

