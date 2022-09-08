function CST
n=input('enter nom. of nodes:\n');
% v=input('enter v:\n');
% E=input('enter E:\n');
% t=input('enter t:\n');
v=0.25;
E=30*10e6;
t=1;
e=input('enter number of elements:\n');
K=zeros(2*n,2*n);
for e=1:e
i=input('enter [Xi Yi node]:\n');
j=input('enter [Xj Yj node]:\n');
m=input('enter [Xm Ym node]:\n');

Bi=j(2)-m(2);
Bj=m(2)-i(2);
Bm=i(2)-j(2);
yi=m(1)-j(1);
yj=i(1)-m(1);
ym=j(1)-i(1);

AA=[1 i(1) i(2);1 j(1) j(2); 1 m(1) m(2)];
% A2=det(AA);
A2=((j(1)*m(2)-m(1)*j(2))-(i(1)*m(2)-m(1)*i(2))+(i(1)*j(2)-j(1)*i(2)));


B=(1/A2)*[Bi 0 Bj 0 Bm 0;0 yi 0 yj 0 ym;yi Bi yj Bj ym Bm];


D=(E/(1-v^2))*[1 v 0;v 1 0;0 0 (1-v)/2];

BT=zeros(6,3);

for l=1:3
for ll=1:6
BT(ll,l)=B(l,ll);
end
end


KK=t*(A2/2)*(BT)*D*B;
g=1;
f=1;
for x=[2*i(3)-1 2*i(3) 2*j(3)-1 2*j(3) 2*m(3)-1 2*m(3)]
    for y=[2*i(3)-1 2*i(3) 2*j(3)-1 2*j(3) 2*m(3)-1 2*m(3)]
     K(x,y)=KK(g,f)+K(x,y);
     f=f+1;
    end
    f=1;
    g=g+1;
end


end
K

%-------- ^^^ End of Stiffness Matrix  ^^^ ---------%

%------------  make forces matrix -------------%
F=zeros(2*n,1);
numf=input('enter number of Forces:\n');

for q=1:numf
    f=input('enter node number for f:\n');
    F(2*f-1,1)=input('enter Fx value:\n');
    F(2*f,1)=input('enter Fy value:\n');
end

% --- for number of displacement:
dd=zeros(2*n,1);
for i=1:2*n
    dd(2*i-1,1)=i;
    dd(2*i)=i;
end
%--------------------------------

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
dd=ddd;
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

%--------------------------------------

nn=size(K,1);
O=zeros(nn,nn+1);
for q=1:nn
    O(q,1:nn)=K(q,:)+O(q,nn);
    O(q,nn+1)=F(q,:)+O(q,nn+1);
end
O

A=O;


%---------------------------------------------
n=size(A,1);

for j=1:n-1
for i=j+1:n

    
    A(i,:)=-(A(i,j)/A(j,j))*A(j,:)+A(i,:);
    

            
end
end
for j=2:n
for i=1:j-1
    A(i,:)=-(A(i,j)/A(j,j))*A(j,:)+A(i,:);
end
end
for i=1:n
    A(i,:)=(1/A(i,i))*A(i,:);
end
d=A(:,n+1)


displacements=zeros(n,2);
displacements(:,1)=d;
displacements(:,2)=dd;

displacements

















 