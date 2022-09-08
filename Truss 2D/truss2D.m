clear
clc

E=xlsread('truss2D.xlsx',1,'B5');
% I=xlsread('frame.xlsx',1,'C5');
% A=xlsread('frame.xlsx',1,'D5');
A_Total=xlsread('truss2D.xlsx',1,'D6:D50');
% n=input('enter the number of nodes:\n');
n=xlsread('truss2D.xlsx',1,'B8');
% e=input('enter the number of elements:\n');
e=xlsread('truss2D.xlsx',1,'C8');
% R0=input('enter the number of support reactions:\n');

%input cordinates X Y n for each node:
coords=xlsread('truss2D.xlsx',1,'G6:L50');
Angles=xlsread('truss2D.xlsx',1,'M6:M50');
%%%%%%%%%%
R0=xlsread('truss2D.xlsx',1,'D8');
Reactions=xlsread('truss2D.xlsx',1,'O6:Q50');

Forces=xlsread('truss2D.xlsx',1,'W6:Y50');

I_Total=xlsread('truss2D.xlsx',1,'E6:E50');
Lenghts=zeros(e,1);
for e=1:e
    Lenghts(e)=sqrt((coords(e,1)-coords(e,4))^2+(coords(e,2)-coords(e,5))^2);
end

K=zeros(2*n,2*n);
F=zeros(2*n,1);
for e=1:e
    node1=coords(e,3);
    node2=coords(e,6);    

    A=A_Total(e);
    I=I_Total(e);
    L=Lenghts(e);
    a=Angles(e);


Ke(:,:,e)=        (A*E/L)*[1    0   -1  0;
                           0    0   0   0;
                           -1   0   1   0;
                           0    0   0   0];
                       
T(:,:,e)=[cosd(a) sind(a)   0           0;
         -sind(a) cosd(a)   0           0;
          0        0        cosd(a)     sind(a);
          0        0        -sind(a)    cosd(a)];


KE(:,:,e)=T(:,:,e)'*Ke(:,:,e)*T(:,:,e);

%assembling:
H=[2*node1-1 2*node1 2*node2-1 2*node2];
for x=1:4
    for y=1:4
        K(H(x),H(y))=KE(x,y,e)+K(H(x),H(y));
    end
end
end
KK=K;
% Boundry Conditions:
g=1;
numR=size(Reactions,1);
for i=1:numR
   if Reactions(i,2)==0
       R(g)=2*Reactions(i,1)-1;
       g=g+1; 
   end
   if Reactions(i,3)==0
       R(g)=2*Reactions(i,1);
       g=g+1; 
   end
end
K(R,:)=[];
K(:,R)=[];


%Forces Matrix:
numF=size(Forces,1);
for i=1:numF
    F(2*Forces(i,1)-1,1)=Forces(i,2);
    F(2*Forces(i,1),1)=Forces(i,3);
end
F(R,:)=[];
%%%%%%%%%%%%
d=K\F;
%%%%%%%%%%%%
d_Total=zeros(2*n,1);
H=1:2*n;
H(R)=[];
%be dast avardane matrise jabejaiye kol - 9*1
for i=1:size(H,2)
    d_Total(H(i))=d(i);
end

%be das avardane matrise d dar har eleman
d_element=zeros(4,1);
for e=1:e
    node1=coords(e,3);
    node2=coords(e,6);
    H=[2*node1-1 2*node1 2*node2-1 2*node2];
    for i=1:4
        d_element(i,1,e)=d_Total(H(i));
    end    
end

%mohasebeye matrise f har eleman dar dastage mahalie khodash
for e=1:e
    f(:,1,e)=T(:,:,e)*Ke(:,:,e)*T(:,:,e)'*d_element(:,1,e);
end

fprintf('\nDisplacements:\n\n')
fprintf('Node num.\t\t X\t\t\t Y\t\t \n');
for i=1:n
    fprintf('\t%d\t\t %7.6f\t %7.6f\t \n',i,d_Total(2*i-1,1),d_Total(2*i,1))
end
fprintf('\n***********************************************************************\n')
fprintf('\nElements Forces:\n\n')
for e=1:e
    fprintf('Element %d\n',e)
    fprintf('F%d X:   %4.3f \nF%d Y:   %4.3f \n\nF%d X:   %4.3f \nF%d Y:   %4.3f \n\n********************** \n'...
        ,coords(e,3),f(1,:,e),coords(e,3),f(2,:,e)...
        ,coords(e,6),f(3,:,e),coords(e,6),f(4,:,e))
end

FF=KK*d_Total