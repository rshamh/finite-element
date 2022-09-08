clear all
clc

E=xlsread('A3.xlsx',1,'C2');
v=xlsread('A3.xlsx',1,'C3');
n=xlsread('A3.xlsx',1,'C4');
e=xlsread('A3.xlsx',1,'C5');

Pz=xlsread('A3.xlsx',1,'C6');
Pr=xlsread('A3.xlsx',1,'C7');

coords=xlsread('A3.xlsx',1,'E3:G1000');
e_coords=xlsread('A3.xlsx',1,'I3:L1000');

Forces=xlsread('A3.xlsx',1,'N3:Q1000');

Reactions=xlsread('A3.xlsx',1,'S3:U1000');

D1=(1-v)/(1-2*v);
D2=(v)/(1-2*v);
D=(E/(1+v))*[D1 D2 D2 0
             D2 D1 D2 0
             D2 D2 D1 0
             0  0  0  0.5];
B=zeros(4,6,e);   
K=zeros(2*n,2*n);
F=zeros(2*n,1);
for e=1:e
    Fe=zeros(6,1);
    
   ii=e_coords(e,2);
   jj=e_coords(e,3);
   kk=e_coords(e,4);
   
   Ri=coords(ii,2);
   Zi=coords(ii,3);
   Rj=coords(jj,2);
   Zj=coords(jj,3);
   Rk=coords(kk,2);
   Zk=coords(kk,3);
   
   Lij=abs(Ri-Rj);
   Ljk=abs(Rj-Rk);
   Lik=abs(Ri-Rk);
%%%%%%%%%%   
A=abs(0.5*det([Ri Zi 1; Rj Zj 1; Rk Zk 1]));  
%%%%%%%%%%

ai=Rj*Zk-Rk*Zj;
aj=Rk*Zi-Ri*Zk;
ak=Ri*Zj-Rj*Zi;

bi=Zj-Zk;
bj=Zk-Zi;
bk=Zi-Zj;

ci=Rk-Rj;
cj=Ri-Rk;
ck=Rj-Ri;

r_bar=(Ri+Rj+Rk)/3;
z_bar=(Zi+Zj+Zk)/3;

Ni=(1/(2*A))*(ai+bi*r_bar+ci*z_bar);
Nj=(1/(2*A))*(aj+bj*r_bar+cj*z_bar);
Nk=(1/(2*A))*(ak+bk*r_bar+ck*z_bar);

% 
% B(:,:,e)=(1/(2*A))*[bi 0 bj 0 bk 0
%           2*A*Ni/r_bar 0 2*A*Nj/r_bar 0 2*A*Nk/r_bar 0
%           0 ci  0 cj 0 ck
%           ci bi cj bj ck bk];
      
B(:,:,e)=(1/(2*A))*[bi 0 bj 0 bk 0
                    0 ci 0  cj 0 ck
                   (ai/r_bar)+bi+(ci*z_bar/r_bar) 0 (aj/r_bar)+bj+(cj*z_bar/r_bar) 0 (ak/r_bar)+bk+(ck*z_bar/r_bar) 0
                    ci bi cj bj ck bk];

      
      
Ke=2*pi*A*r_bar*B(:,:,e)'*D*B(:,:,e);

if Forces(e,3)==1
    Fe=(2*pi*Lij*(1/6))*[(2*Ri+Rj)*Pr ; (2*Ri+Rj)*Pz ; (Ri+2*Rj)*Pr ; (Ri+2*Rj)*Pz; 0; 0];
end

g1=1;
for x=[2*ii-1 2*ii 2*jj-1 2*jj 2*kk-1 2*kk]
    g2=1;
    F(x)=Fe(g1)+F(x);
    for y=[2*ii-1 2*ii 2*jj-1 2*jj 2*kk-1 2*kk]
        K(x,y)=Ke(g1,g2)+K(x,y);
        g2=g2+1;
    end
    g1=g1+1;
end



end
num_R=size(Reactions,1);

for i=1:num_R
    i1=Reactions(i,1);
    j1=2*i1-1;
    j2=2*i1;
     if Reactions(i,2)==1
        K(j1,:)=0;
        K(:,j1)=0;
        K(j1,j1)=1;
    end
     if Reactions(i,3)==1
        K(j2,:)=0;
        K(:,j2)=0;
        K(j2,j2)=1;
     end
end



d=K\F;



for e=1:e
    d_element=zeros(6,1);
    
    g3=1;
    for x=[2*ii-1 2*ii 2*jj-1 2*jj 2*kk-1 2*kk];
        d_element(g3)=d(x);
        g3=g3+1;
    end
    
    S(:,:,e)=D*B(:,:,e)*d_element;
    
end

disp(' ')
disp(' *************************')
disp('  ')
fprintf('Displacements:\n')
disp('  ')

fprintf('Node \t Ur \t    Uz \n')
for i=1:n
    Ur(i)=d(2*i-1);
    Uz(i)=d(2*i);
    fprintf('%d \t %f \t %f \n',i, d(2*i-1), d(2*i))
end

disp('  ')
fprintf('_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-\n')
disp(' ')

fprintf('Stress:\n\n')
fprintf('Element       Sigma r        Sigma z       Sigma teta      tau rz \n\n')
for i=1:e
   fprintf('%d            %g        %g       %g       %g \n\n',i , S(1,:,i) , S(2,:,i) , S(3,:,i), S(4,:,i))
end
  
S_r=zeros(n,1);
S_z=zeros(n,1);
S_teta=zeros(n,1);
S_rz=zeros(n,1);

for e=1:e
    ii=e_coords(e,2);
    jj=e_coords(e,3);
    kk=e_coords(e,4);
    
    S_r(ii)=S(1,:,e)+S_r(ii);
    S_r(jj)=S(1,:,e)+S_r(jj);
    S_r(kk)=S(1,:,e)+S_r(kk);

    S_z(ii)=S(2,:,e)+S_z(ii);
    S_z(jj)=S(2,:,e)+S_z(jj);
    S_z(kk)=S(2,:,e)+S_z(kk);
    
    S_teta(ii)=S(3,:,e)+S_teta(ii);
    S_teta(jj)=S(3,:,e)+S_teta(jj);
    S_teta(kk)=S(3,:,e)+S_teta(kk);
    
    S_rz(ii)=S(4,:,e)+S_rz(ii);
    S_rz(jj)=S(4,:,e)+S_rz(jj);
    S_rz(kk)=S(4,:,e)+S_rz(kk);
end

%miangin giri:
num_n=zeros(n,1);
for n=1:n
    num_n(n)=size(find(e_coords(:,[2 3 4])==n),2);
end
 
ave_r=zeros(n,1);
ave_z=zeros(n,1);
ave_teta=zeros(n,1);
ave_rz=zeros(n,1);

for n=1:n
    ave_r(n)=S_r(n)/num_n(n);
    ave_z(n)=S_z(n)/num_n(n);
    ave_teta(n)=S_teta(n)/num_n(n);
    ave_rz(n)=S_rz(n)/num_n(n);
end

%mohasebat jahate rasme nemoodar:

[meshX,meshY]=meshgrid([0 7.5 12.5 17.5],[0 -5 -10 -15]);
Z_r=zeros(4,4);
Z_z=zeros(4,4);
Z_teta=zeros(4,4);
Z_rz=zeros(4,4);
Z_Ur=zeros(4,4);
Z_Uz=zeros(4,4);

for n=1:n
    Z_r(n)=ave_r(n);
    Z_z(n)=ave_z(n);
    Z_teta(n)=ave_teta(n);
    Z_rz(n)=ave_rz(n);
    
    Z_Ur(n)=Ur(n);
    Z_Uz(n)=Uz(n);
end
    Z_r=Z_r';
    Z_z=Z_z';
    Z_teta=Z_teta';
    Z_rz=Z_rz';

    
    % rasmee nemoodarha:
figure
contourf(meshX,meshY,Z_r)
title('Sigma r')
xlabel('r')
ylabel('z')
hold on
scatter(coords(:,2),coords(:,3))
colorbar

figure
contourf(meshX,meshY,Z_z)
title('Sigma z')
xlabel('r')
ylabel('z')
hold on
scatter(coords(:,2),coords(:,3))
colorbar

figure
contourf(meshX,meshY,Z_teta)
title('Sigma teta')
xlabel('r')
ylabel('z')
hold on
scatter(coords(:,2),coords(:,3))
colorbar

figure
contourf(meshX,meshY,Z_rz)
title('tau rz')
xlabel('r')
ylabel('z')
hold on
scatter(coords(:,2),coords(:,3))
colorbar

figure
contourf(meshX,meshY,Z_Ur)
title('Displacement Ur')
xlabel('r')
ylabel('z')
hold on
scatter(coords(:,2),coords(:,3))
colorbar

figure
contourf(meshX,meshY,Z_Uz)
title('Displacement Uz')
xlabel('r')
ylabel('z')
hold on
scatter(coords(:,2),coords(:,3))
colorbar