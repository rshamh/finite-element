% function isopara8
clear
clc
E=xlsread('iso8.xlsx',1,'B2');
v=xlsread('iso8.xlsx',1,'B3');
tt=xlsread('iso8.xlsx',1,'B4');
n=xlsread('iso8.xlsx',1,'B5');
e=xlsread('iso8.xlsx',1,'B6');

nodes=xlsread('iso8.xlsx',1,'D2:D1000');
cordX=xlsread('iso8.xlsx',1,'E2:E1000');
cordY=xlsread('iso8.xlsx',1,'F2:F1000');
e_cords=xlsread('iso8.xlsx',1,'O2:V1000');

RRR(:,[2 3])=xlsread('iso8.xlsx',1,'G2:H1000');
RRR(:,1)=xlsread('iso8.xlsx',1,'D2:D1000');

Fp=xlsread('iso8.xlsx',1,'J2:L1000');

all_nodes=zeros(1,8*e);
ggg=1;

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

a=diff(Yst,t);
b=diff(Yst,s);
c=diff(Xst,s);
d=diff(Xst,t);


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

det_J=abs(c*a-b*d);

JB=(1/det_J)*B;


K=zeros(2*n,2*n);
D=(E/(1-v^2))*[1 v 0;v 1 0; 0 0 (1-v)/2];

%%%%%% Gauss Nodes
gauss=0.774596669241483;
S=[-gauss -gauss -gauss    0   0    0    gauss gauss gauss];
T=[-gauss  0      gauss -gauss 0  gauss -gauss   0   gauss];
W=[0.308642 0.493827 0.308642 0.493827 0.790123 0.493827 0.308642 0.493827 0.308642];

%%%%%% Stifness Matrix for each element
for e=1:e
    B_element=zeros(3,16,9);
    
ii=[cordX(e_cords(e,1)) cordY(e_cords(e,1)) nodes(e_cords(e,1))];
jj=[cordX(e_cords(e,2)) cordY(e_cords(e,2)) nodes(e_cords(e,2))];
mm=[cordX(e_cords(e,3)) cordY(e_cords(e,3)) nodes(e_cords(e,3))];
zz=[cordX(e_cords(e,4)) cordY(e_cords(e,4)) nodes(e_cords(e,4))];
ij=[cordX(e_cords(e,5)) cordY(e_cords(e,5)) nodes(e_cords(e,5))];
jm=[cordX(e_cords(e,6)) cordY(e_cords(e,6)) nodes(e_cords(e,6))];
mz=[cordX(e_cords(e,7)) cordY(e_cords(e,7)) nodes(e_cords(e,7))];
zi=[cordX(e_cords(e,8)) cordY(e_cords(e,8)) nodes(e_cords(e,8))];

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


Ke=zeros(16,16);

for j=1:9
s=S(j);
t=T(j);
w=W(j);
B=eval(JB);
det_JJ=eval(det_J);
B_element(:,:,j)=B;
Ke=(B'*D*B*tt*w*det_JJ)+Ke;
end
Ke;
B_Total{e}=B_element;
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
% shomareye tamami gereh ha be tartibe estefade:
for i=[ii(3) jj(3) mm(3) zz(3) ij(3) jm(3) mz(3) zi(3)]
    all_nodes(ggg)=i;
    ggg=ggg+1;
end

end
% tedade nodhaye bekar rafte baraye miangin girie tanesh ha:
num_nodes=zeros(n,2);
for i=1:n
    O=find(all_nodes==i);
    num_nodes(i,1)=i;
    num_nodes(i,2)=size(O,2);
end

K


%-------- ^^^ End of Stiffness Matrix  ^^^ ---------%

%------------  make forces matrix -------------%
F=zeros(2*n,1);
numf=size(Fp,1);

for q=1:numf
    F(2*Fp(q,1)-1,1)=Fp(q,2);
    F(2*Fp(q,1),1)=Fp(q,3);
end

% --- for displacement Node's number:
dd=zeros(2*n,1);
for i=1:n
    dd(2*i-1,1)=i;
    dd(2*i)=i;
end




g=1;
for i=1:n
   if RRR(i,2)==0
       R(g)=2*RRR(i,1)-1;
       g=g+1; 
   end
   if RRR(i,3)==0
       R(g)=2*RRR(i,1);
       g=g+1; 
   end
end
numR=size(R,2);
kk_row=zeros(2*n-numR,2*n);
FF=zeros(2*n-numR,1);
ddd=zeros(2*n-numR,1);

qq=1;


s=0;
for q=1:2*n
   for i=1:numR
        if q==R(i)
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
K=kk_row;
kk_col=zeros(2*n-numR,2*n-numR);
qq=1;
s=0;
for q=1:2*n
   for i=1:numR
        if q==R(i)
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


d_total=zeros(2*n,1);
R_T=1:2*n;
R_T(R)=[];
for i=1:2*n-numR
    d_total(R_T(i))=d(i);
end

Sx=zeros(n,1);
Sy=zeros(n,1);
Sxy=zeros(n,1);

% %%%% ----------   Stress  ---------  %%%%

T=Trans();
       
for e=1:e
   d_element=zeros(16,1);

   ii=e_cords(e,1);
   jj=e_cords(e,2);
   mm=e_cords(e,3);
   zz=e_cords(e,4);
   ij=e_cords(e,5);
   jm=e_cords(e,6);
   mz=e_cords(e,7);
   zi=e_cords(e,8);
   
      g=1;
   for i=[2*ii-1 2*ii 2*jj-1 2*jj 2*mm-1 2*mm 2*zz-1 2*zz 2*ij-1 2*ij 2*jm-1 2*jm 2*mz-1 2*mz 2*zi-1 2*zi ]
      d_element(g)=d_total(i);
      g=g+1;
   end
   
   BB=B_Total{e};
   
   for i=1:9    
     SS(:,:,i)=D*BB(:,:,i)*d_element;
   end
   
      % SSx: sigma x           SSy: sigma y       SSt: tanesh boreshi
   
   for i=1:9
       SSx(i,1)=SS(1,:,i);
       SSy(i,1)=SS(2,:,i);
       SSxy(i,1)=SS(3,:,i);
   end

   Sxe=T*SSx;
   Sye=T*SSy;
   Sxye=T*SSxy;
   
   g=1;
   for i=[ii jj mm zz ij jm mz zi]
       Sx(i)=Sxe(g)+Sx(i);
       g=g+1;
   end
 
      g=1;
   for i=[ii jj mm zz ij jm mz zi]
       Sy(i)=Sye(g)+Sy(i);
       g=g+1;
   end
   
      g=1;
   for i=[ii jj mm zz ij jm mz zi]
       Sxy(i)=Sxye(g)+Sxy(i);
       g=g+1;
   end   
end

Sx_Final=zeros(n,2);
Sy_Final=zeros(n,2);
Sxy_Final=zeros(n,2);


for i=1:n
    Sx_Final(i,1)=i;
    Sy_Final(i,1)=i;
    Sxy_Final(i,1)=i;
    
    Sx_Final(i,2)=Sx(i)/num_nodes(i,2);
    Sy_Final(i,2)=Sy(i)/num_nodes(i,2);
    Sxy_Final(i,2)=Sxy(i)/num_nodes(i,2);
end
  
Sx_Final;
Sy_Final;
Sxy_Final;

fprintf('\n Displacements: \n\n');
for i=1:2*n
    node=ceil(i/2);
    if rem(i,2)==1
        fprintf('%d \t Ux \t %6.5f \n',node,d_total(i)) ;
    end
    if rem(i,2)==0
        fprintf('%d \t Uy \t %6.5f \n\n',node,d_total(i)) ;
    end    
end

fprintf('\n ******************************************** \n');


fprintf('Sress X: \n\n')
for i=1:n
    fprintf('%d \t %6.5f \n',i,Sx_Final(i,2))
end
disp('***************************************')
fprintf('Stres Y: \n\n')
for i=1:n
    fprintf('%d \t %6.5f \n',i,Sy_Final(i,2))
end
disp('***************************************')
fprintf('Stres XY: \n\n')
for i=1:n
    fprintf('%d \t %6.5f \n',i,Sxy_Final(i,2))
end
       

Eresult{1,1}=['num of node'];
Eresult{1,2}=['Ux'];
Eresult{1,3}=['Uy'];
Eresult{1,4}=['Sigma X'];
Eresult{1,5}=['Sigma Y'];
Eresult{1,6}=['Tou XY'];

for i=1:n
    Eresult{i+1,1}=i;
    Eresult{i+1,2}=d_total(2*i-1);
    Eresult{i+1,3}=d_total(2*i);
    Eresult{i+1,4}=Sx_Final(i,2);
    Eresult{i+1,5}=Sy_Final(i,2);
    Eresult{i+1,6}=Sxy_Final(i,2);
end
xlswrite('result.xlsx',Eresult,'Res','A1');

