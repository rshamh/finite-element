clear
clc

inp=xlsread('Q4.xlsx',1,'A2:F2');
E=inp(1); v=inp(2); thickness=inp(3); n=inp(4); e=inp(5);

Coords=xlsread('Q4.xlsx',1,'C7:E1000');
nodes=zeros(n,1);
X=zeros(n,1);
Y=zeros(n,1);

for i=1:n
   nodes(i)=Coords(i,1);
   X(i)=Coords(i,2);
   Y(i)=Coords(i,3);
end

e_coords=xlsread('Q4.xlsx',1,'I7:L1000');

Fp=xlsread('Q4.xlsx',1,'T7:V1000');

D=(E/(1-v^2))*[1 v  0     ;
               v 1  0     ;
               0 0 (1-v)/2];        
           
Suport=xlsread('Q4.xlsx',1,'O7:Q1000');
numS=size(Suport,1);

Gauss=xlsread('Q4.xlsx',1,'Y7:AA10');

K=zeros(2*n,2*n);
F=zeros(2*n,1);

for e=1:e
    Ke=zeros(8,8);
    B=zeros(3,8,4);
    
    P1=e_coords(e,1);
    P2=e_coords(e,2);
    P3=e_coords(e,3);
    P4=e_coords(e,4);
    
    for i=1:4      
          s(i)=Gauss(i,2);
          t(i)=Gauss(i,3);
                
          Ns=[ 0.25*(t(i)-1)  0.25*(1-t(i)) 0.25*(t(i)+1)  0.25*(-t(i)-1)];         
          Nt=[ 0.25*(s(i)-1)   0.25*(-s(i)-1) 0.25*(s(i)+1)  0.25*(1-s(i))];
              
          Xs(i) =0.25*((t(i)-1)*X(P1) + (1-t(i))*X(P2)  + (1+t(i))*X(P3) + (-1-t(i))*X(P4));
          Xt(i) =0.25*((s(i)-1)*X(P1) + (-1-s(i))*X(P2) + (1+s(i))*X(P3) + (1-s(i))*X(P4));
          Ys(i) =0.25*((t(i)-1)*Y(P1) + (1-t(i))*Y(P2)  + (1+t(i))*Y(P3) + (-1-t(i))*Y(P4));
          Yt(i) =0.25*((s(i)-1)*Y(P1) + (-1-s(i))*Y(P2) + (1+s(i))*Y(P3) + (1-s(i))*Y(P4));
          
          det_j=((Xs(i)*Yt(i))-(Ys(i)*Xt(i)));
          
          B(1,1,i)=Ns(1)*Yt(i)-Ys(i)*Nt(1);
          B(3,1,i)=Xs(i)*Nt(1)-Ns(1)*Xt(i);
          B(2,2,i)=Xs(i)*Nt(1)-Ns(1)*Xt(i);
          B(3,2,i)=Ns(1)*Yt(i)-Ys(i)*Nt(1);
 
          B(1,3,i)=Ns(2)*Yt(i)-Ys(i)*Nt(2);
          B(3,3,i)=Xs(i)*Nt(2)-Ns(2)*Xt(i);
          B(2,4,i)=Xs(i)*Nt(2)-Ns(2)*Xt(i);
          B(3,4,i)=Ns(2)*Yt(i)-Ys(i)*Nt(2);
          
          B(1,5,i)=Ns(3)*Yt(i)-Ys(i)*Nt(3);
          B(3,5,i)=Xs(i)*Nt(3)-Ns(3)*Xt(i);
          B(2,6,i)=Xs(i)*Nt(3)-Ns(3)*Xt(i);
          B(3,6,i)=Ns(3)*Yt(i)-Ys(i)*Nt(3);         
          
          B(1,7,i)=Ns(4)*Yt(i)-Ys(i)*Nt(4);
          B(3,7,i)=Xs(i)*Nt(4)-Ns(4)*Xt(i);
          B(2,8,i)=Xs(i)*Nt(4)-Ns(4)*Xt(i);
          B(3,8,i)=Ns(4)*Yt(i)-Ys(i)*Nt(4); 
          
          B(:,:,i)=B(:,:,i)/det_j;
          
          Ke=((B(:,:,i)')*D*B(:,:,i)*det_j*thickness)+Ke;
    end
    B_T{e}=B;
    
H=[2*P1-1  2*P1  2*P2-1  2*P2  2*P3-1  2*P3  2*P4-1  2*P4];    
    for x=1:8 
        for y=1:8
            K(H(x),H(y))=Ke(x,y)+K(H(x),H(y));
        end
    end
end


numf=size(Fp,1);
for q=1:numf
    F(2*Fp(q,1)-1,1)=Fp(q,2);
    F(2*Fp(q,1),1)=Fp(q,3);
end


R=[];
for i=1:size(Suport,1)
    Rlabel=Suport(i,1);
    res=Suport(i,2:3);
    for j=1:2
        if res(j)==0
            R1=2*(Rlabel-1)+j;
            R=[R R1];
        end
    end
end
   
K(R,:)=[];
K(:,R)=[];
F(R,:)=[];

d=K\F;

d_total=zeros(2*n,1);
R_T=1:2*n;
R_T(R)=[];
numR=size(R,2);
for i=1:2*n-numR
    d_total(R_T(i))=d(i);
end

T =[ 1.8660   -0.5000   -0.5000    0.1340; -0.5000    0.1340    1.8660   -0.5000; 0.1340   -0.5000   -0.5000    1.8660; -0.5000    1.8660    0.1340   -0.5000];

for e=1:e
    d_element=zeros(8,1);
      
   g=1;
   for i=[2*e_coords(e,1)-1 2*e_coords(e,1) 2*e_coords(e,2)-1 2*e_coords(e,2) 2*e_coords(e,3)-1 2*e_coords(e,3)  2*e_coords(e,4)-1 2*e_coords(e,4)]
       d_element(g)=d_total(i);
       g=g+1;
   end
   
   B=B_T{e};
   
   for i=1:4    
   Stresses(:,:,i)=D*B(:,:,i)*d_element;
   end
   
   for i=1:4
       Sx(i,1)=Stresses(1,:,i);
       Sy(i,1)=Stresses(2,:,i);
       Sxy(i,1)=Stresses(3,:,i);
   end
  
   disp('*************')
   fprintf('Element %d \n',e)
   disp('*************')
   Sigma_x=T*Sx
   Sigma_y=T*Sy
   tau_xy=T*Sxy
   disp('---------------------')
end
disp('  ');
fprintf('Displacement Of Node\n');
disp('  ');
disp('Node             U             V');
disp('  ');
 for i=1:n
     fprintf('%3g',i); fprintf('= %13G',d_total(2*i-1)); fprintf(' %13G\n',d_total(2*i));
     
 end
fprintf(' -------------------------------\n');
disp('  ');
          