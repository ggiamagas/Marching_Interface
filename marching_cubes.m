function [P,C,S,UI,NUV,NU,Bar,E,Time]=marching_cubes(phi,psi,x,y,z,u,v,w,Nx,Ny,Nz)
%Marching Cubes main output is P C and Time
%P: Matrix of all interfacial points determined by means of linear interpolation. 
%Each point occupies a row of the matrix while each column corresponds to the x y and z coordinate respectively
%C: Connectivity Matrix determines which points (rows of P matrix) will be triangulated i.e. connected with straight lines
%S: Surfactant Concentration on the interface
%The way the points are stored in P simplifies the connectivity, every three rows represent the vertices of one triangle
%NU: Matrix containing all the Unit Normals to all the Triangles. Using quiver with vector end point the coordinates 
%and vector start point any point on the triangle draws the unit normals at the specified point
%NUV: Unit Normals on all triangle edges computed through the gradient of Phi
%UI: The velocity vector interpolated on all triangle vertices.
%Bar: Matrix containing all the Triangle Barycentra(typically Centroids). These points are selected as the starting points of the unit normals
%Time: The elapsed time of the algorithm. Two essential findings in terms of speeding up: 
%A. P and C must be global matrices over the whole grid and rendered only at the end rather
%than at each cube. B. Matlab is very much (!) faster in columnwise rather than rowise reading. In other words,
%the growing index of the matrices should preferably be the column number. Therefore the transpose of P and C is 
%initially saved and at the end their proper form for the triangulation is recovered via a simple tranponse.

T=tri_matrix(); %Call Triangulation Matrix
isolevel=0;    %Choose the isolevel

tStart=tic;
q=1;
r=1;
p=1;
f=1;
for i=1:Nx-1   %Loop over a total of (Nx-1)*(Ny-1)*(Nz-1) Cubes 
    tic
for j=1:Ny-1
for k=1:Nz-1

%Phi
phi0=phi(i,k+1,j+1);
phi1=phi(i+1,k+1,j+1);
phi2=phi(i+1,k+1,j);
phi3=phi(i,k+1,j);

phi4=phi(i,k,j+1);
phi5=phi(i+1,k,j+1);
phi6=phi(i+1,k,j);
phi7=phi(i,k,j);

Phi=[phi0,phi1,phi2,phi3,phi4,phi5,phi6,phi7];


%Velocities
u0=u(i,k+1,j+1);
u1=u(i+1,k+1,j+1);
u2=u(i+1,k+1,j);
u3=u(i,k+1,j);

u4=u(i,k,j+1);
u5=u(i+1,k,j+1);
u6=u(i+1,k,j);
u7=u(i,k,j);

U=[u0,u1,u2,u3,u4,u5,u6,u7];


v0=v(i,k+1,j+1);
v1=v(i+1,k+1,j+1);
v2=v(i+1,k+1,j);
v3=v(i,k+1,j);

v4=v(i,k,j+1);
v5=v(i+1,k,j+1);
v6=v(i+1,k,j);
v7=v(i,k,j);

V=[v0,v1,v2,v3,v4,v5,v6,v7];


w0=w(i,k+1,j+1);
w1=w(i+1,k+1,j+1);
w2=w(i+1,k+1,j);
w3=w(i,k+1,j);

w4=w(i,k,j+1);
w5=w(i+1,k,j+1);
w6=w(i+1,k,j);
w7=w(i,k,j);

W=[w0,w1,w2,w3,w4,w5,w6,w7];


%Gradients
if (i==1)
    Gx0=(phi(i+1,k+1,j+1)-phi(i,k+1,j+1))/(x(i+1)-x(i));
    Gx3=(phi(i+1,k+1,j)-phi(i,k+1,j))/(x(i+1)-x(i));
    Gx4=(phi(i+1,k,j+1)-phi(i,k,j+1))/(x(i+1)-x(i));
    Gx7=(phi(i+1,k,j)-phi(i,k,j))/(x(i+1)-x(i));    
else
    Gx0=(phi(i+1,k+1,j+1)-phi(i-1,k+1,j+1))/(2*(x(i+1)-x(i)));
    Gx3=(phi(i+1,k+1,j)-phi(i-1,k+1,j))/(2*(x(i+1)-x(i)));
    Gx4=(phi(i+1,k,j+1)-phi(i-1,k,j+1))/(2*(x(i+1)-x(i)));
    Gx7=(phi(i+1,k,j)-phi(i-1,k,j))/(2*(x(i+1)-x(i)));  
end

if (i==Nx-1)
    Gx1=(phi(i+1,k+1,j+1)-phi(i,k+1,j+1))/(x(i+1)-x(i));
    Gx2=(phi(i+1,k+1,j)-phi(i,k+1,j))/(x(i+1)-x(i));
    Gx5=(phi(i+1,k,j+1)-phi(i,k,j+1))/(x(i+1)-x(i));
    Gx6=(phi(i+1,k,j)-phi(i,k,j))/(x(i+1)-x(i));
else
    Gx1=(phi(i+2,k+1,j+1)-phi(i,k+1,j+1))/(2*(x(i+2)-x(i+1)));
    Gx2=(phi(i+2,k+1,j)-phi(i,k+1,j))/(2*(x(i+2)-x(i+1)));
    Gx5=(phi(i+2,k,j+1)-phi(i,k,j+1))/(2*(x(i+2)-x(i+1)));
    Gx6=(phi(i+2,k,j)-phi(i,k,j))/(2*(x(i+2)-x(i+1)));
end

if (j==1)
    Gy2=(phi(i+1,k+1,j+1)-phi(i+1,k+1,j))/(y(j+1)-y(j));
    Gy3=(phi(i,k+1,j+1)-phi(i,k+1,j))/(y(j+1)-y(j));
    Gy6=(phi(i+1,k,j+1)-phi(i+1,k,j))/(y(j+1)-y(j));
    Gy7=(phi(i,k,j+1)-phi(i,k,j))/(y(j+1)-y(j));      
else
    Gy2=(phi(i+1,k+1,j+1)-phi(i+1,k+1,j-1))/(2*(y(j+1)-y(j)));
    Gy3=(phi(i,k+1,j+1)-phi(i,k+1,j-1))/(2*(y(j+1)-y(j)));
    Gy6=(phi(i+1,k,j+1)-phi(i+1,k,j-1))/(2*(y(j+1)-y(j)));
    Gy7=(phi(i,k,j+1)-phi(i,k,j-1))/(2*(y(j+1)-y(j))); 
end

if (j==Ny-1)
    Gy0=(phi(i,k+1,j+1)-phi(i,k+1,j))/(y(j+1)-y(j));
    Gy1=(phi(i+1,k+1,j+1)-phi(i+1,k+1,j))/(y(j+1)-y(j));
    Gy4=(phi(i,k,j+1)-phi(i,k,j))/(y(j+1)-y(j));
    Gy5=(phi(i+1,k,j+1)-phi(i+1,k,j))/(y(j+1)-y(j)); 
else
    Gy0=(phi(i,k+1,j+2)-phi(i,k+1,j))/(2*(y(j+2)-y(j+1)));
    Gy1=(phi(i+1,k+1,j+2)-phi(i+1,k+1,j))/(2*(y(j+2)-y(j+1)));
    Gy4=(phi(i,k,j+2)-phi(i,k,j))/(2*(y(j+2)-y(j+1)));
    Gy5=(phi(i+1,k,j+2)-phi(i+1,k,j))/(2*(y(j+2)-y(j+1)));
end

if (k==1)
    Gz4=(phi(i,k+1,j+1)-phi(i,k,j+1))/(z(k+1)-z(k));
    Gz5=(phi(i+1,k+1,j+1)-phi(i+1,k,j+1))/(z(k+1)-z(k));
    Gz6=(phi(i+1,k+1,j)-phi(i+1,k,j))/(z(k+1)-z(k));
    Gz7=(phi(i,k+1,j)-phi(i,k,j))/(z(k+1)-z(k));
else
    Gz4=(phi(i,k+1,j+1)-phi(i,k-1,j+1))/(z(k+1)-z(k-1));
    Gz5=(phi(i+1,k+1,j+1)-phi(i+1,k-1,j+1))/(z(k+1)-z(k-1));
    Gz6=(phi(i+1,k+1,j)-phi(i+1,k-1,j))/(z(k+1)-z(k-1));
    Gz7=(phi(i,k+1,j)-phi(i,k-1,j))/(z(k+1)-z(k-1));
end

if (k==Nz-1)
    Gz0=(phi(i,k+1,j+1)-phi(i,k,j+1))/(z(k+1)-z(k));
    Gz1=(phi(i+1,k+1,j+1)-phi(i+1,k,j+1))/(z(k+1)-z(k));
    Gz2=(phi(i+1,k+1,j)-phi(i+1,k,j))/(z(k+1)-z(k));
    Gz3=(phi(i,k+1,j)-phi(i,k,j))/(z(k+1)-z(k));   
else
    Gz0=(phi(i,k+2,j+1)-phi(i,k,j+1))/(z(k+2)-z(k));
    Gz1=(phi(i+1,k+2,j+1)-phi(i+1,k,j+1))/(z(k+2)-z(k));
    Gz2=(phi(i+1,k+2,j)-phi(i+1,k,j))/(z(k+2)-z(k));
    Gz3=(phi(i,k+2,j)-phi(i,k,j))/(z(k+2)-z(k));    
end

Gx=[Gx0,Gx1,Gx2,Gx3,Gx4,Gx5,Gx6,Gx7];
Gy=[Gy0,Gy1,Gy2,Gy3,Gy4,Gy5,Gy6,Gy7];
Gz=[Gz0,Gz1,Gz2,Gz3,Gz4,Gz5,Gz6,Gz7];
G=[Gx;Gy;Gz];
L=sqrt(Gx.^2+Gy.^2+Gz.^2);

%Cube Vertex Unit Normals 
VN=zeros(3,8);
for o=1:8
   VN(:,o)=G(:,o)/L(o);
end


%Psi 
psi0=psi(i,k+1,j+1);
psi1=psi(i+1,k+1,j+1);
psi2=psi(i+1,k+1,j);
psi3=psi(i,k+1,j);

psi4=psi(i,k,j+1);
psi5=psi(i+1,k,j+1);
psi6=psi(i+1,k,j);
psi7=psi(i,k,j);

Psi=[psi0,psi1,psi2,psi3,psi4,psi5,psi6,psi7];


%Coordinates
x0=x(i);
y0=y(j+1);
z0=z(k+1);
x1=x(i+1);
y1=y(j+1);
z1=z(k+1);
x2=x(i+1);
y2=y(j);
z2=z(k+1);
x3=x(i);
y3=y(j);
z3=z(k+1);

x4=x(i);
y4=y(j+1);
z4=z(k);
x5=x(i+1);
y5=y(j+1);
z5=z(k);
x6=x(i+1);
y6=y(j);
z6=z(k);
x7=x(i);
y7=y(j);
z7=z(k);

X=[x0,x1,x2,x3,x4,x5,x6,x7];
Y=[y0,y1,y2,y3,y4,y5,y6,y7];
Z=[z0,z1,z2,z3,z4,z5,z6,z7];

%Bin for Binary
for m=1:8 
    if Phi(m)<isolevel
        Bin(m)=1;
    else
        Bin(m)=0;
    end
end



%Cubil ID
cubeindex=1*Bin(1)+2*Bin(2)+4*Bin(3)+8*Bin(4)+16*Bin(5)+32*Bin(6)+64*Bin(7)+128*Bin(8); %Cube ID based on its vertices 
Tri=T(cubeindex+1,:);


%Interpolation
for l=1:+3:13
    if Tri(l)~=-1
        [P1,S1,IN1,U1x,U1y,U1z]=interpolator(Tri(l),isolevel,Phi,VN,Psi,U,V,W,X,Y,Z);
        [P2,S2,IN2,U2x,U2y,U2z]=interpolator(Tri(l+1),isolevel,Phi,VN,Psi,U,V,W,X,Y,Z);
        [P3,S3,IN3,U3x,U3y,U3z]=interpolator(Tri(l+2),isolevel,Phi,VN,Psi,U,V,W,X,Y,Z);
        P(:,q)=P1;            %Points (this is the tranpose of the actual P!)
        P(:,q+1)=P2;
        P(:,q+2)=P3;
        C(:,r) = [q q+1 q+2]; %Connectivity (this is the tranpose of the actual C!)
        
        AB=[P2(1)-P1(1);P2(2)-P1(2);P2(3)-P1(3)]; %P1P2
        BC=[P3(1)-P2(1);P3(2)-P2(2);P3(3)-P2(3)]; %P2P3
        CA=[P1(1)-P3(1);P1(2)-P3(2);P1(3)-P3(3)]; %P3P1
        
        %Triangle normal [or cross(BC,CA) or cross(AB,CA)]
        n=cross(AB,BC); 
        NU(:,r)=n/norm(n);
        %O(:,r)=[P2(1),P2(2),P2(3)];      
        
        %Barycentra
        Barx=(P1(1)+P2(1)+P3(1))/3;
        Bary=(P1(2)+P2(2)+P3(2))/3;
        Barz=(P1(3)+P2(3)+P3(3))/3;
        Bar(:,r)=[Barx,Bary,Barz];
        
        %Triangle surface
        E(r)=1/2*(norm(cross(AB,CA))); 
        
        %Interpolated unit normals
        NUV(:,q)=IN1;
        NUV(:,q+1)=IN2;
        NUV(:,q+2)=IN3;  
        
        %Interpolated velocity vector
        UI(:,q)=[U1x U1y U1z];
        UI(:,q+1)=[U2x U2y U2z];
        UI(:,q+2)=[U3x U3y U3z];
        
        %Mean surfactant concentration on the triangle        
        SM=(S1+S2+S3)/3;
        S(q)=SM;            
        S(q+1)=SM;
        S(q+2)=SM;
                  
        q=q+3;
        r=r+1;   
    end
   
end
 
end
end
    time(i)=toc %Elapsed time per x-step (total Nx-1 steps) (show on command window helps keep track of the running speed)
end

Time=toc(tStart); %Overall elapsed time

end
