function Tr=Trans()
clear
clc
gauss=0.774596669241483;
S=[-gauss -gauss -gauss    0   0    0    gauss gauss gauss];
T=[-gauss  0      gauss -gauss 0  gauss -gauss   0   gauss];

for i=1:9
    X(i,:)=[1 S(i) T(i) S(i)^2 S(i)*T(i) T(i)^2 (S(i)^2)*T(i) S(i)*(T(i)^2) (S(i)^2)*(T(i)^2)];
end
a=1;
S=[-a a  a  -a   0    a   0 -a ];
T=[-a -a a   a  -a    0   a  0 ];
for i=1:8
    Y(i,:)=[1 S(i) T(i) S(i)^2 S(i)*T(i) T(i)^2 (S(i)^2)*T(i) S(i)*(T(i)^2) (S(i)^2)*(T(i)^2)];
end

x=inv(X);
Tr=Y*x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%