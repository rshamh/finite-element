clear all
clc
fprintf('                                ******** 3D Frame Matrix Analysis ********\n\n');
fprintf('Program is running...\n');

%read from Excel:
E_Total=xlsread('frame.xlsx',1,'C6:C50');
A_Total=xlsread('frame.xlsx',1,'D6:D50');

Iy_Total=xlsread('frame.xlsx',1,'E6:E50');
Iz_Total=xlsread('frame.xlsx',1,'F6:F50');
J_Total=xlsread('frame.xlsx',1,'G6:G50');
G_Total=xlsread('frame.xlsx',1,'H6:H50');

nodes=xlsread('frame.xlsx',1,'A5');
e=xlsread('frame.xlsx',1,'A7');

%input cordinates:
coords=xlsread('frame.xlsx',1,'J6:Q50');

l_Total=xlsread('frame.xlsx',1,'R6:R50');
m_Total=xlsread('frame.xlsx',1,'S6:S50');
n_Total=xlsread('frame.xlsx',1,'T6:T50');
D_Total=xlsread('frame.xlsx',1,'U6:U50');

Reactions=xlsread('frame.xlsx',1,'W6:AC50');
Forces=xlsread('frame.xlsx',1,'AH6:AN50');

L_Total=xlsread('frame.xlsx',1,'V6:V100');

%Calculate K for each Element:
F=zeros(6*nodes,1);
K=zeros(6*nodes,6*nodes);

for e=1:e
    node(1)=coords(e,4);
    node(2)=coords(e,8);    

    Iy=Iy_Total(e);
    Iz=Iz_Total(e);
    J=J_Total(e);
    G=G_Total(e);
    l=l_Total(e);
    m=m_Total(e);
    n=n_Total(e);
    D=D_Total(e);
    
    L=L_Total(e);
    A=A_Total(e);
    E=E_Total(e);
    
Ke(:,:,e)=         [(A*E)/L           0                           0                 0               0                 0                 -(A*E)/L             0                       0                   0               0                  0             ;
                    0                (12*E*Iz)/(L^3)              0                 0               0                (6*E*Iz)/(L^2)       0                -(12*E*Iz)/(L^3)          0                   0               0                 (6*E*Iz)/(L^2) ;
                    0                 0                          (12*E*Iy)/(L^3)    0             -(6*E*Iy)/(L^2)     0                   0                  0                     -(12*E*Iy)/(L^3)      0             -(6*E*Iy)/(L^2)      0             ;
                    0                 0                           0                (G*J)/L          0                 0                   0                  0                       0                 -(G*J)/L          0                  0             ;
                    0                 0                         -(6*E*Iy)/(L^2)     0              (4*E*Iy)/L         0                   0                  0                      (6*E*Iy)/(L^2)       0              (2*E*Iy)/L          0             ;
                    0                (6*E*Iz)/(L^2)               0                 0               0                (4*E*Iz)/(L)         0                -(6*E*Iz)/(L^2)           0                   0               0                 (2*E*Iz)/(L)   ;
                    
                  -(A*E)/L            0                           0                 0               0                 0                  (A*E)/L             0                       0                   0               0                  0             ;
                    0               -(12*E*Iz)/(L^3)              0                 0               0               -(6*E*Iz)/(L^2)       0                 (12*E*Iz)/(L^3)          0                   0               0                -(6*E*Iz)/(L^2) ;
                    0                 0                         -(12*E*Iy)/(L^3)    0              (6*E*Iy)/(L^2)     0                   0                  0                      (12*E*Iy)/(L^3)      0              (6*E*Iy)/(L^2)      0             ;
                    0                 0                           0               -(G*J)/L          0                 0                   0                  0                       0                  (G*J)/L          0                  0             ;
                    0                 0                         -(6*E*Iy)/(L^2)     0              (2*E*Iy)/L         0                   0                  0                      (6*E*Iy)/(L^2)       0              (4*E*Iy)/L          0             ;
                    0                (6*E*Iz)/(L^2)               0                 0               0                (2*E*Iz)/(L)         0                -(6*E*Iz)/(L^2)           0                   0               0                 (4*E*Iz)/(L)   ];
                    
                    
                    

  if ~(n==1 && l==0 && m==0) 
      if ~(n==-1 && l==0 && m==0)
T(:,:,e)=           [l           m       n       0       0       0       0       0       0       0       0       0;
                   -(m/D)        l/D     0       0       0       0       0       0       0       0       0       0;
                   -(l*n)/D    -(m*n)/D  D       0       0       0       0       0       0       0       0       0;
                     0           0       0       l       m       n       0       0       0       0       0       0;
                     0           0       0     -(m/D)    l/D     0       0       0       0       0       0       0;
                     0           0       0     -(l*n)/D -(m*n)/D D       0       0       0       0       0       0;
                     0           0       0       0       0       0       l       m       n       0       0       0;
                     0           0       0       0       0       0     -(m/D)    l/D     0       0       0       0;
                     0           0       0       0       0       0     -(l*n)/D -(m*n)/D D       0       0       0;
                     0           0       0       0       0       0       0       0       0       l       m       n;
                     0           0       0       0       0       0       0       0       0     -(m/D)    l/D     0;
                     0           0       0       0       0       0       0       0       0     -(l*n)/D -(m*n)/D D];
      end
      
  elseif (n==1 && l==0 && m==0)
      
T(:,:,e)=           [0           0       1       0       0       0       0       0       0       0       0       0;
                     0           1       0       0       0       0       0       0       0       0       0       0;
                    -1           0       0       0       0       0       0       0       0       0       0       0;
                     0           0       0       0       0       1       0       0       0       0       0       0;
                     0           0       0       0       1       0       0       0       0       0       0       0;
                     0           0       0      -1       0       0       0       0       0       0       0       0;
                     0           0       0       0       0       0       0       0       1       0       0       0;
                     0           0       0       0       0       0       0       1       0       0       0       0;
                     0           0       0       0       0       0      -1       0       0       0       0       0;
                     0           0       0       0       0       0       0       0       0       0       0       1;
                     0           0       0       0       0       0       0       0       0       0       1       0;
                     0           0       0       0       0       0       0       0       0      -1       0       0];
      
  elseif (n==-1 && l==0 && m==0)
 T(:,:,e)=          [0           0       -1      0       0       0       0       0       0       0       0       0;
                     0           1       0       0       0       0       0       0       0       0       0       0;
                     1           0       0       0       0       0       0       0       0       0       0       0;
                     0           0       0       0       0      -1       0       0       0       0       0       0;
                     0           0       0       0       1       0       0       0       0       0       0       0;
                     0           0       0       1       0       0       0       0       0       0       0       0;
                     0           0       0       0       0       0       0       0      -1       0       0       0;
                     0           0       0       0       0       0       0       1       0       0       0       0;
                     0           0       0       0       0       0       1       0       0       0       0       0;
                     0           0       0       0       0       0       0       0       0       0       0      -1;
                     0           0       0       0       0       0       0       0       0       0       1       0;
                     0           0       0       0       0       0       0       0       0       1       0       0];
  end
                     
                     
                     
                     
                     
KE(:,:,e)=T(:,:,e)'*Ke(:,:,e)*T(:,:,e);

%assembling:
H=[6*node(1)-5 6*node(1)-4 6*node(1)-3 6*node(1)-2 6*node(1)-1 6*node(1) 6*node(2)-5 6*node(2)-4 6*node(2)-3 6*node(2)-2 6*node(2)-1 6*node(2)];
for x=1:12
    for y=1:12
        K(H(x),H(y))=KE(x,y,e)+K(H(x),H(y));
    end
end
end
KK=K;


% Boundry Conditions:
R=[];
numR=size(Reactions,1);
for i=1:numR
   if Reactions(i,2)==0
       R=[R,6*Reactions(i,1)-5]; 
   end
   if Reactions(i,3)==0
       R=[R,6*Reactions(i,1)-4];
   end
   if Reactions(i,4)==0
       R=[R,6*Reactions(i,1)-3];
   end
   if Reactions(i,5)==0
       R=[R,6*Reactions(i,1)-2];
   end
   if Reactions(i,6)==0
       R=[R,6*Reactions(i,1)-1];
   end
   if Reactions(i,7)==0
       R=[R,6*Reactions(i,1)];
   end
   
end


%Forces Matrix:
numF=size(Forces,1);
for i=1:numF
    
    F(6*Forces(i,1)-5,1)=Forces(i,2);
    F(6*Forces(i,1)-4,1)=Forces(i,3);
    F(6*Forces(i,1)-3,1)=Forces(i,4);
    F(6*Forces(i,1)-2,1)=Forces(i,5);
    F(6*Forces(i,1)-1,1)=Forces(i,6);
    F(6*Forces(i,1),1)=Forces(i,7);
    
end

K(R,:)=[];
K(:,R)=[];
F(R,:)=[];
%%%%%%%%%%%%

%---------------   Start Gauss-Jordan Eq solution   -------------------
d=K\F;
% ------------ End Gauss-Jordan Eq solution   ------------------------

%%%%%%%%%%%%
d_Total=zeros(6*nodes,1);
H=1:6*nodes;
H(R)=[];

%Calculate Total Displacement Matrix:
for i=1:size(H,2)
    d_Total(H(i))=d(i);
end

%Calculate Displacement Matrix for each Elements:
d_element=zeros(12,1);
for e=1:e
    node(1)=coords(e,4);
    node(2)=coords(e,8);
    H=[6*node(1)-5 6*node(1)-4 6*node(1)-3 6*node(1)-2 6*node(1)-1 6*node(1) 6*node(2)-5 6*node(2)-4 6*node(2)-3 6*node(2)-2 6*node(2)-1 6*node(2)];
    for i=1:12
        d_element(i,1,e)=d_Total(H(i));
    end    
    
    %mohasebeye matrise f har eleman dar dastage mahalie khodash
    f(:,1,e)=Ke(:,:,e)*T(:,:,e)*d_element(:,1,e);
    
end

fprintf('.........................................................................................................\n\n');
fprintf('\t\t\t\t\t\t\t***********************************************\n');
fprintf('\t\t\t\t\t\t\t***                                         ***\n');
fprintf('\t\t\t\t\t\t\t***                Results                  ***\n');
fprintf('\t\t\t\t\t\t\t***                                         ***\n');
fprintf('\t\t\t\t\t\t\t***********************************************\n');


fprintf('\nDisplacements:\n\n')
fprintf('Node num.\t\tUX\t\t\tUY\t\t\tUZ\t\t    RX\t\t    RY\t\t    RZ\n');
fprintf('---------\t --------\t --------\t --------\t --------\t --------\t --------\n');
for i=1:nodes
    fprintf('\t%d\t\t %7.6f\t %7.6f\t %7.6f\t %7.6f\t %7.6f\t %7.6f\n',i,d_Total(6*i-5,1),d_Total(6*i-4,1),d_Total(6*i-3,1),d_Total(6*i-2,1),d_Total(6*i-1,1),d_Total(6*i,1))
end
fprintf('\n***********************************************************************\n')
fprintf('\nElements Forces:\n\n')
for e=1:e
    fprintf('Element %d\n',e)
    fprintf('F%d X:   %4.3f \nF%d Y:   %4.3f \nF%d Z:   %4.3f \nM%d X:   %4.3f \nM%d Y:   %4.3f \nM%d Z:   %4.3f \n\nF%d X:   %4.3f \nF%d Y:   %4.3f \nF%d Z:   %4.3f \nM%d X:   %4.3f \nM%d Y:   %4.3f \nM%d Z:   %4.3f \n\n********************** \n'...
        ,coords(e,4),f(1,:,e),coords(e,4),f(2,:,e),coords(e,4),f(3,:,e)...
        ,coords(e,4),f(4,:,e),coords(e,4),f(5,:,e),coords(e,4),f(6,:,e)...
        ,coords(e,8),f(7,:,e),coords(e,8),f(8,:,e),coords(e,8),f(9,:,e)...
        ,coords(e,8),f(10,:,e),coords(e,8),f(11,:,e),coords(e,8),f(12,:,e))
end

FF=KK*d_Total;