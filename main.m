function []=main()
clear;
close all;
clc;

cd ./input

%Grid Specifications
Lx=4*pi;
Ly=2*pi;
h=1; 
Nx=512;
Nz=513;
Ny=256;
xp = fopen('x.dat','r');
zp = fopen('z.dat','r');
yp= fopen('y.dat','r');
x=fread(xp,'double','ieee-le');
z=fread(zp,'double','ieee-le');
y=fread(yp,'double','ieee-le');
fclose(xp);
fclose(zp);
fclose(yp);


%Phase and Velocity Fields Reading
l=0;
for m=1
  l=l+1;
  filename = sprintf('phi.dat', m);
  uname=sprintf('u.dat', m);
  vname=sprintf('v.dat', m);
  wname=sprintf('w.dat', m);
  p = fopen(filename,'r');
  up=fopen(uname,'r');
  vp=fopen(vname,'r');
  wp=fopen(wname,'r');  
  A=zeros(Nx,Nz,Ny);
  B=zeros(Nx,Nz,Ny);
  C=zeros(Nx,Nz,Ny);
  D=zeros(Nx,Nz,Ny);
  phir=fread(p,'double','ieee-le');
  ur=fread(up,'double','ieee-le');
  vr=fread(vp,'double','ieee-le');
  wr=fread(wp,'double','ieee-le');
  A=reshape(phir,[Nx Nz Ny]);
  B=reshape(ur,[Nx Nz Ny]);
  C=reshape(vr,[Nx Nz Ny]);
  D=reshape(wr,[Nx Nz Ny]);
  phi(:,:,:,l)=A(:,:,:);
  u(:,:,:,l)=B(:,:,:);
  v(:,:,:,l)=C(:,:,:);
  w(:,:,:,l)=D(:,:,:);
  fclose(p);
  fclose(up);
  fclose(vp);
  fclose(wp);
end
lmax=l;

cd ..

%Surfactant Field Reading
psi=rand(Nx,Nz,Ny,lmax);

%Call to the main function 
[P,C,S,U,NUV,NU,Bar,E,Time]=marching_cubes(phi,psi,x,y,z,u,v,w,Nx,Ny,Nz);
                                      

disp('Elapsed time')
T=Time  %Show Overall Elapsed Time
disp('Interface Surface Area')
E=sum(E)  %Show Total Interface Surface Area

%Drawing
figures(P,C,S,U,NUV,NU,Bar,E,Lx,Ly,h);


end