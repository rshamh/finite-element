% function isopara4_2
clear
clc
XXX=xlsread('iso8.xlsx',1);
nx=XXX(:,1);
xx=XXX(:,2);
yy=XXX(:,3);
n=max(nx)
e=(size(nx,1)/8)


syms s t x1 x2 x3 x4 x5 x6 x7 x8 y1 y2 y3 y4 y5 y6 y7 y8 

N(1)=((1-s)*(1-t)*(-s-t-1))/4;
N(2)=((1+s)*(1-t)*(s-t-1))/4;
N(3)=((1+s)*(1+t)*(s+t-1))/4;
N(4)=((1-s)*(1+t)*(-s+t-1))/4;
N(5)=((1-t)*(1+s)*(1-s))/2;
N(6)=((1+s)*(1+t)*(1-t))/2;
N(7)=((1+t)*(1+s)*(1-s))/2;
N(8)=((1-s)*(1+t)*(1-t))/2;

Xst=0.25*((-1+s^2+s*t+t^2-(s^2)*t-s*(t^2))*x1+(-1+s^2-s*t+t^2-(s^2)*t+s*(t^2))*x2+(-1+s^2+s*t+t^2+(s^2)*t+s*(t^2))*x3+(-1+s^2-s*t+t^2+(s^2)*t-s*(t^2))*x4+(2-2*t-2*s^2+2*(s^2)*t)*x5+(2+2*s-2*t^2-2*(t^2)*s)*x6+(2+2*t-2*s^2-2*(s^2)*t)*x7+(2-2*s-2*t^2+2*(t^2)*s)*x8);
Yst=0.25*((-1+s^2+s*t+t^2-(s^2)*t-s*(t^2))*y1+(-1+s^2-s*t+t^2-(s^2)*t+s*(t^2))*y2+(-1+s^2+s*t+t^2+(s^2)*t+s*(t^2))*y3+(-1+s^2-s*t+t^2+(s^2)*t-s*(t^2))*y4+(2-2*t-2*s^2+2*(s^2)*t)*y5+(2+2*s-2*t^2-2*(t^2)*s)*y6+(2+2*t-2*s^2-2*(s^2)*t)*y7+(2-2*s-2*t^2+2*(t^2)*s)*y8);

a=diff(Yst,t)
b=diff(Yst,s)
c=diff(Xst,s)
d=diff(Xst,t)




% a=(1/4)*(y1*(s-1)+y2*(-1-s)+y3*(1+s)+y4*(1-s));
% b=(1/4)*(y1*(t-1)+y2*(1-t)+y3*(1+t)+y4*(-1-t));
% c=(1/4)*(x1*(t-1)+x2*(1-t)+x3*(1+t)+x4*(-1-t));
% d=(1/4)*(x1*(s-1)+x2*(-1-s)+x3*(1+s)+x4*(1-s));

B=vpa(zeros(3,16));

for i=1:8
dNs=diff(N(i),s);
dNt=diff(N(i),t);
B(1,2*i-1)=(a*dNs)-(b*dNt);
B(2,2*i-1)=0;
B(3,2*i-1)=c*(dNt)-d*(dNs);
B(1,2*i)=0;
B(2,2*i)=c*(dNt)-d*(dNs);
B(3,2*i)=a*(dNs)-b*(dNt);
end
B;

det_J=c*a-b*d;

JB=det_J*B;

E=input('enter E:\n');
v=input('enter v:\n');
tt=input('enter thickness:\n');
K=zeros(2*n,2*n);
D=(E/(1-v^2))*[1 v 0;v 1 0; 0 0 (1-v)/2];

%%%%%% Gauss Nodes
gauss=0.57735026918962;
S=[-gauss gauss gauss -gauss 0 gauss 0 -gauss];
T=[-gauss -gauss gauss gauss -gauss 0 gauss 0];

%%%%%% Stifness Matrix for each element
for e=1:e
ii=[xx(8*e-7) yy(8*e-7) nx(8*e-7)];
jj=[xx(8*e-6) yy(8*e-6) nx(8*e-6)];
mm=[xx(8*e-5) yy(8*e-5) nx(8*e-5)];
zz=[xx(8*e-4) yy(8*e-4) nx(8*e-4)];
ij=[xx(8*e-3) yy(8*e-3) nx(8*e-3)];
jm=[xx(8*e-2) yy(8*e-2) nx(8*e-2)];
mz=[xx(8*e-1) yy(8*e-1) nx(8*e-1)];
zi=[xx(8*e) yy(8*e) nx(8*e)];

x1=ii(1);
x2=jj(1);
x3=mm(1);
x4=zz(1);
x5=ij(1);
x6=jm(1);
x7=mz(1);
x8=zi(1);

y1=ii(2);
y2=jj(2);
y3=mm(2);
y4=zz(2);
y5=ij(2);
y6=jm(2);
y7=mz(2);
y8=zi(2);

% J=eval(J);
Ke=zeros(16,16);

for j=1:8
s=S(j);
t=T(j);
B=eval(JB);
Ke=(B'*D*B*tt)+Ke;
end
Ke;

%%%% Assembel Stiffness Matrix
g=1;
f=1;
for x=[2*ii(3)-1 2*ii(3) 2*jj(3)-1 2*jj(3) 2*mm(3)-1 2*mm(3) 2*zz(3)-1 2*zz(3) 2*ij(3)-1 2*ij(3) 2*jm(3)-1 2*jm(3) 2*mz(3)-1 2*mz(3) 2*zi(3)-1 2*zi(3) ]
    for y=[2*ii(3)-1 2*ii(3) 2*jj(3)-1 2*jj(3) 2*mm(3)-1 2*mm(3) 2*zz(3)-1 2*zz(3) 2*ij(3)-1 2*ij(3) 2*jm(3)-1 2*jm(3) 2*mz(3)-1 2*mz(3) 2*zi(3)-1 2*zi(3) ]
     K(x,y)=Ke(g,f)+K(x,y);
     f=f+1;
    end
    f=1;
    g=g+1;
end
end

K

pause
%-------- ^^^ End of Stiffness Matrix  ^^^ ---------%

%------------  make forces matrix -------------%
F=zeros(2*n,1);
numf=input('enter number of Forces:\n');

for q=1:numf
    f=input('enter node number for f:\n');
    F(2*f-1,1)=input('enter Fx value:\n');
    F(2*f,1)=input('enter Fy value:\n');
end

% --- for displacement Node's number:
dd=zeros(2*n,1);
for i=1:n
    dd(2*i-1,1)=i;
    dd(2*i)=i;
end
numR=input('enter number of R:\n');
kk_row=zeros(2*n-numR,2*n);
FF=zeros(2*n-numR,1);
ddd=zeros(2*n-numR,1);
qq=1;
r=input('[enter node numbers in Stiffness matrix for R:]\n');
s=0;
for q=1:2*n
   for i=1:numR
        if q==r(i)
            s=1;
        end
   end
   if s==0
          kk_row(qq,:)=K(q,:)+kk_row(qq,:);
          FF(qq,:)=F(q,:)+FF(qq,:);
          ddd(qq,:)=dd(q,:)+ddd(qq,:);
          qq=qq+1;
   end
   s=0;
end
F=FF;
% dd=ddd;
K=kk_row;
kk_col=zeros(2*n-numR,2*n-numR);
qq=1;
s=0;
for q=1:2*n
   for i=1:numR
        if q==r(i)
            s=1;
        end
   end    
   if s==0
       kk_col(:,qq)=K(:,q)+kk_col(:,qq);
       qq=qq+1;
   end
   s=0;
end

K=kk_col
F

%---------------   Start Gauss-Jordan Eq solution   -------------------
nn=size(K,1);
O=zeros(nn,nn+1);
for q=1:nn
    O(q,1:nn)=K(q,:)+O(q,nn);
    O(q,nn+1)=F(q,:)+O(q,nn+1);
end
O;
%------------------------
A=O;
nA=size(A,1);
for j=1:nA-1
for i=j+1:nA    
    A(i,:)=-(A(i,j)/A(j,j))*A(j,:)+A(i,:);          
end
end
for j=2:nA
for i=1:j-1
    A(i,:)=-(A(i,j)/A(j,j))*A(j,:)+A(i,:);
end
end
for i=1:nA
    A(i,:)=(1/A(i,i))*A(i,:);
end
d=A(:,nA+1)
% ------------ End Gauss-Jordan Eq solution   ------------------------
format long
displacements=zeros(nA,2);
displacements(:,1)=d;
displacements(:,2)=ddd;

displacements

%%%% ----------   Stress  ---------  %%%%
dTotal=zeros(2*n,2);
qq=1;
q=0;
for i=1:2*n

   for j=1:numR
       if i==r(j)
           dTotal(i,1)=0;
           q=1;
       end
   end
   
   if q==0
       dTotal(i,1)=d(qq);
       qq=qq+1;
   end
  q=0;   
end
dTotal(:,2)=dd;
dTotal
pause
for e=1:e
ii=[xx(4*e-3) yy(4*e-3) nx(4*e-3)];
jj=[xx(4*e-2) yy(4*e-2) nx(4*e-2)];
mm=[xx(4*e-1) yy(4*e-1) nx(4*e-1)];
zz=[xx(4*e) yy(4*e) nx(4*e)];
x1=ii(1);
x2=jj(1);
x3=mm(1);
x4=zz(1);
y1=ii(2);
y2=jj(2);
y3=mm(2);
y4=zz(2);

for j=1:4
s=S(j);
t=T(j);
B=eval(JB);
switch j
    case 1
        B1=B;
    case 2
        B2=B;
    case 3
        B3=B;
    case 4
        B4=B;
end
end
ne=zeros(8,1);
q=1;
for nd=[ii(3) jj(3) mm(3) zz(3)]
[a,b]=find(dTotal==nd);
a=dTotal(a);
ne([q q+1],1)=[a(1) a(2)];
% ne(q+1,1)=a(2)
q=q+2;
end
ne
S1=D*B1*ne
S2=D*B2*ne
S3=D*B3*ne
S4=D*B4*ne
pause
end

