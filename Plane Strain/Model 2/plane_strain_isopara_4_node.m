% function isopara4_2
clear
clc

E=xlsread('iso.xlsx',1,'B6');
v=xlsread('iso.xlsx',1,'B7');
tt=xlsread('iso.xlsx',1,'B8');
n=xlsread('iso.xlsx',1,'B9');
e=xlsread('iso.xlsx',1,'B10');

RRR=xlsread('iso.xlsx',1,'H6:J1000');
numR=size(RRR,1);

nodes=xlsread('iso.xlsx',1,'D6:D1000');
Xcords=xlsread('iso.xlsx',1,'E6:E1000');
Ycords=xlsread('iso.xlsx',1,'F6:F1000');

e_cords=xlsread('iso.xlsx',1,'N6:Q1000');
e_cords=e_cords';
%bare noghtei:
Fp=xlsread('iso.xlsx',1,'S6:U1000');
 
D=(E/((1+v)*(1-2*v)))*[1-v v 0;v 1-v 0; 0 0 (1-2*v)/2];



 %%%%%% Gauss 2x2
gauss=0.57735026918962;
S=[-gauss gauss gauss -gauss];
T=[-gauss -gauss gauss gauss];

K=zeros(2*n,2*n);
F=zeros(2*n,1);

all_nodes=zeros(1,4*e);
ggg=1;

for e=1:e
    Ke=zeros(8,8);
    B_element=zeros(3,8,4);
    
   ii=[ Xcords(e_cords(1,e)) Ycords(e_cords(1,e)) e_cords(1,e) ];
   jj=[ Xcords(e_cords(2,e)) Ycords(e_cords(2,e)) e_cords(2,e) ];
   mm=[ Xcords(e_cords(3,e)) Ycords(e_cords(3,e)) e_cords(3,e) ];
   zz=[ Xcords(e_cords(4,e)) Ycords(e_cords(4,e)) e_cords(4,e) ];
   
   x1=ii(1);
   x2=jj(1);
   x3=mm(1);
   x4=zz(1);
   y1=ii(2);
   y2=jj(2);
   y3=mm(2);
   y4=zz(2);

  Fe=zeros(8,1); 
for j=1:4
    s=S(j);
    t=T(j);
   
dXs=0.25*((t-1)*x1+(1-t)*x2+(1+t)*x3+(-1-t)*x4);
dXt=0.25*((s-1)*x1+(-1-s)*x2+(1+s)*x3+(1-s)*x4);
dYs=0.25*((t-1)*y1+(1-t)*y2+(1+t)*y3+(-1-t)*y4);
dYt=0.25*((s-1)*y1+(-1-s)*y2+(1+s)*y3+(1-s)*y4);

N=0.25*[(1-s)*(1-t) (1+s)*(1-t) (1+s)*(1+t) (1-s)*(1+t)];

dNs=[0.25*(t-1) 0.25*(1-t) 0.25*(t+1) 0.25*(-t-1)];
dNt=[0.25*(s-1) 0.25*(-s-1) 0.25*(s+1) 0.25*(1-s)];

det_J=abs((dXs*dYt)-(dXt*dYs));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%          B & stiffness Matrix for Elements  -------------------------------------
    B=zeros(3,8);
    for i=1:4
       a=(dYt*dNs(i)-dYs*dNt(i))/det_J;
       b=(-dXt*dNs(i)+dXs*dNt(i))/det_J;
       B([1,3],2*i-1)=[a,b];
       B([2,3],2*i)=[b,a];
    end    
    B_element(:,:,j)=B;
    Ke=((B')*D*B*det_J)+Ke;



end
B_Total{e}=B_element;



%%%% Assembel Stiffness Matrix    -------------------------------------
g=1;
f=1;
for x=[2*ii(3)-1 2*ii(3) 2*jj(3)-1 2*jj(3) 2*mm(3)-1 2*mm(3) 2*zz(3)-1 2*zz(3)]
    for y=[2*ii(3)-1 2*ii(3) 2*jj(3)-1 2*jj(3) 2*mm(3)-1 2*mm(3) 2*zz(3)-1 2*zz(3)]
     K(x,y)=Ke(g,f)+K(x,y);
     f=f+1;
    end
    f=1;
    g=g+1;
end

% shomareye tamami gereh ha be tartibe estefade:
for i=[ii(3) jj(3) mm(3) zz(3)]
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
for i=1:numR
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

K(R,:)=[];
K(:,R)=[];
F(R,:)=[];


d=K\F;


d_total=zeros(2*n,1);
R_T=1:2*n;
R_T(R)=[];

for i=1:2*n-numR
    d_total(R_T(i))=d(i);
end


Sx=zeros(n,1);
Sy=zeros(n,1);
Sxy=zeros(n,1);
Sz=zeros(n,1);




T=[1+(sqrt(3)/2) -0.5 -0.5 1-(sqrt(3)/2) ;
    -0.5 1-(sqrt(3)/2) 1+(sqrt(3)/2) -0.5 ;
    1-(sqrt(3)/2) -0.5 -0.5 1+(sqrt(3)/2) ;
    -0.5 1+(sqrt(3)/2) 1-(sqrt(3)/2) -0.5 ];


for e=1:e
    d_element=zeros(8,1);
    
   ii=e_cords(1,e);
   jj=e_cords(2,e);
   mm=e_cords(3,e);
   zz=e_cords(4,e);
   
   g=1;
   for i=[2*ii-1 2*ii 2*jj-1 2*jj 2*mm-1 2*mm  2*zz-1 2*zz]
       d_element(g)=d_total(i);
       g=g+1;
   end
   
   BB=B_Total{e};
   
   for i=1:4    
   SS(:,:,i)=D*BB(:,:,i)*d_element;
   ee(:,:,i)=BB(:,:,i)*d_element;   % ee: strain(x,y,yz)
   end
   % SSx: sigma x           SSy: sigma y       SSt: tanesh boreshi     
   
   for i=1:4
       SSx(i,1)=SS(1,:,i);
       SSy(i,1)=SS(2,:,i);
       SSxy(i,1)=SS(3,:,i);
       
       eex(i,1)=ee(1,:,i);
       eey(i,1)=ee(2,:,i);
   end
   
   % SSz: Sigma z
   SSz=(E/(1+v))*(v/(1-2*v))*(eex+eey);
   
   
   Sxe=T*SSx;
   Sye=T*SSy;
   Sxye=T*SSxy;
   
   Sze=T*SSz; 
  
   g=1;
   for i=[ii jj mm zz]
       Sx(i)=Sxe(g)+Sx(i);
       g=g+1;
   end
 
      g=1;
   for i=[ii jj mm zz]
       Sy(i)=Sye(g)+Sy(i);
       g=g+1;
   end
   
      g=1;
   for i=[ii jj mm zz]
       Sxy(i)=Sxye(g)+Sxy(i);
       g=g+1;
   end
   
         g=1;
   for i=[ii jj mm zz]
       Sz(i)=Sze(g)+Sz(i);
       g=g+1;
   end
   
end


Sx_Final=zeros(n,2);
Sy_Final=zeros(n,2);
Sxy_Final=zeros(n,2);

Sz_Final=zeros(n,2);

for i=1:n
    Sx_Final(i,1)=i;
    Sy_Final(i,1)=i;
    Sxy_Final(i,1)=i;
    Sz_Final(i,1)=i;
    
    Sx_Final(i,2)=Sx(i)/num_nodes(i,2);
    Sy_Final(i,2)=Sy(i)/num_nodes(i,2);
    Sxy_Final(i,2)=Sxy(i)/num_nodes(i,2);
    Sz_Final(i,2)=Sz(i)/num_nodes(i,2);
end




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
fprintf('Stress X: \n\n')
for i=1:n
    fprintf('%d \t %6.5f \n',i,Sx_Final(i,2))
end
disp('***************************************')
fprintf('Stress Y: \n\n')
for i=1:n
    fprintf('%d \t %6.5f \n',i,Sy_Final(i,2))
end
disp('***************************************')
fprintf('Tau XY: \n\n')
for i=1:n
    fprintf('%d \t %6.5f \n',i,Sxy_Final(i,2))
end
disp('***************************************')
fprintf('Stress Z: \n\n')
for i=1:n
    fprintf('%d \t %6.5f \n',i,Sz_Final(i,2))
end







