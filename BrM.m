function [ A ] = BrM(D,T)
%生成二维布朗运动散点
%   A返回(0,DT)时间内的布朗运动点,第一列为时间，第二列为位置
%  T单位运动时长
%  D总步数
C=randn(D,2).*sqrt(T);
A(1,:)=C(1,:);
for i=2:D
    
    A(i,:)=A(i-1,:)+C(i,:);
    

end

