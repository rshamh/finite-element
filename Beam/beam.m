clear
clc

E=xlsread('beam1.xlsx',1,'B5');
% I=xlsread('frame.xlsx',1,'C5');
% A=xlsread('frame.xlsx',1,'D5');
A_Total=xlsread('beam1.xlsx',1,'D6:D50');
% n=input('enter the number of nodes:\n');
n=xlsread('beam1.xlsx',1,'B8');
% e=input('enter the number of elements:\n');
e=xlsread('beam1.xlsx',1,'C8');
% R0=input('enter the number of support reactions:\n');

%input cordinates X Y n for each node:
coords=xlsread('beam1.xlsx',1,'G6:L50');
Angles=xlsread('beam1.xlsx',1,'M6:M50');
%%%%%%%%%%
R0=xlsread('beam1.xlsx',1,'D8');
Reactions=xlsread('beam1.xlsx',1,'O6:Q50');

Forces=xlsread('beam1.xlsx',1,'T6:V50');

I_Total=xlsread('beam1.xlsx',1,'E6:E50');
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


Ke(:,:,e)=        (2*E*I/(L^3))*[6          3*L         -6        3*L;
                                 3*L        2*(L^2)     -3*L      L^2;
                                 -6         -3*L        6         -3*L;
                                 3*L        L^2         -3*L      2*(L^2)];
                             

%assembling:
H=[2*node1-1 2*node1 2*node2-1 2*node2];
for x=1:4
    for y=1:4
        K(H(x),H(y))=Ke(x,y,e)+K(H(x),H(y));
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
%be dast avardane matrise jabejaiye kol - 4*1
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

for e=1:e
    f(:,1,e)=Ke(:,:,e)'*d_element(:,1,e);
end


fprintf('\nDisplacements:\n\n')
fprintf('Node num.\t\t Y\t\t   Teta\t\t \n');
for i=1:n
    fprintf('\t%d\t\t %7.6f\t %7.6f\t \n',i,d_Total(2*i-1,1),d_Total(2*i,1))
end
fprintf('\n***********************************************************************\n')
fprintf('\nElements Forces:\n\n')
for e=1:e
    fprintf('Element %d\n',e)
    fprintf('V%d  :   %4.3f \nM%d  :   %4.3f \n\nV%d  :   %4.3f \nM%d  :   %4.3f\n********************** \n'...
        ,coords(e,3),f(1,:,e),coords(e,3),f(2,:,e)...
        ,coords(e,6),f(3,:,e),coords(e,6),f(4,:,e))
end

FF=KK*d_Total