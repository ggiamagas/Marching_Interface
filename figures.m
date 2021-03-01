function []=figures(P,C,S,U,NUV,NU,Bar,E,Lx,Ly,h)
%Contains various visualization options and interface rendering 
%A.Simple isosurface rendering with colormap indicating height. The triangles can either be set to be shown or not. If they are the interface
%will look black due to the millions of lines due to the very fine %discretization.
%B.The barycentra of all the triangles used for triangulation give another kind of surface view
%C.Unit normals calculated at every triangle i.e. flat subsurface by means of taking the cross product of its edges vectors (the differences
%between the triangle vertices and dividing by its magnitude (which keep in mind that is also the surface area of the parallelogram defined by the parent vectors). 
%The normals are then rendered with quiver by spacifying the beginning of each vector to be the barycentra of the triangle. Note that the three
%coordinates of the unit normal are simply containing the information about the vectors direction i.e. any vector being vertical to the parent
%vectors. The start of the vector is an arbitrary selection.
%D.Surfactant concentration map. A colormap created to show the surfactant concentration on the interface by using a matrix of 
%all the mean surfactant concentrations over each triangle which has been computed as the average between the vertices surfactant value which in
%their turn are found by means of linear interpolation from the grid vertices on the triangle at the same moment that these become known
%through the marching cubes interpolations.



%Render the Isosurface!
figure;
TR=triangulation(C',P');   %Triangulation (Notice the Transpose of C and P Matrices !!!)
SURFACE=trisurf(TR);             %Render the 3D isosurface
%S=trisurf(TR,'FaceAlpha',0.5); %Create a semi-tranparent (phantom!) surface
% set(S,'LineStyle','none'); %Remove the thousands of triangle edges that hide the colormap due to the grid density
%title('Interface Deformation At High Surface Tension (We=0.1)')
title('Interface Deformation At Low Surface Tension (We=1)')
colormap;
c = colorbar;
ylabel(c, 'Deformation amplitude (A/h)')
Max_Elevation=max(P(3,:));
Min_Elevation=min(P(3,:));
% caxis([0,0.6]);
gu= gca; gu.XAxis.Visible = 'off'; gu.YAxis.Visible = 'off'; gu.ZAxis.Visible = 'off'; grid off
xlim([0 Lx])
ylim([0 Ly])
zlim([-h h])
daspect([1 1 1])
xlabel('x')
ylabel('y')
zlabel('z')  
%zoom on



%Render the Isosurface With Triangulation Visible!
figure;
TR=triangulation(C',P');   %Triangulation (Notice the Transpose of C and P Matrices !!!)
SURFACE=trisurf(TR);             %Render the 3D isosurface
%S=trisurf(TR,'FaceAlpha',0.5); %Create a semi-tranparent (phantom!) surface
%set(S,'LineStyle','none'); %Remove the thousands of triangle edges that hide the colormap due to the grid density
%title('Interface Deformation At High Surface Tension (We=0.1)')
title('Interface Deformation At Low Surface Tension (We=1)')
colormap;
c = colorbar;
ylabel(c, 'Deformation amplitude (A/h)')
max(P(3,:));
min(P(3,:));
%caxis([-0.5,0.6]);
xlim([0 Lx])
ylim([0 Ly])
zlim([-h h])
daspect([1 1 1])
xlabel('x')
ylabel('y')
zlabel('z') 
grid on



%Triangle Barycentra!
% figure;
% Barycentra=Bar';
% plot3(Barycentra(:,1),Barycentra(:,2),Barycentra(:,3),'o')
% title('Triangulation Barycentra')
% xlim([0 Lx])
% ylim([0 Ly])
% zlim([-h h])
% daspect([1 1 1])
% xlabel('x')
% ylabel('y')
% zlabel('z')  


%Unit Normals!
figure;
N=NU';
Barycentra=Bar';
quiver3(Barycentra(:,1),Barycentra(:,2),Barycentra(:,3),N(:,1),N(:,2),N(:,3));
title('Unit Normals At The Barycentra')
xlim([0 Lx])
ylim([0 Ly])
zlim([-h h])
daspect([1 1 1])
xlabel('x')
ylabel('y')
zlabel('z')  



%Render the Isosurface With Unit Normals!
figure;
TR=triangulation(C',P');   %Triangulation 
SURFACE=trisurf(TR);             %Render the 3D isosurface
hold on
N=NU';
Barycentra=Bar';
quiver3(Barycentra(:,1),Barycentra(:,2),Barycentra(:,3),N(:,1),N(:,2),N(:,3));
%set(S,'LineStyle','none'); %Remove the thousands of triangle edges that hide the colormap due to the grid density
title('Interface Rendering Plus Unit Normals')
xlim([0 Lx])
ylim([0 Ly])
zlim([-h h])
daspect([1 1 1])
xlabel('x')
ylabel('y')
zlabel('z')  
% zoom on


%Computed Unit Normals!
figure;
%TR=triangulation(C',P');   %Triangulation (Notice the Transpose of C and P Matrices !!!)
%S=trisurf(TR);             %Render the 3D isosurface
%hold on
NUV=NUV';
O=P';
quiver3(O(:,1),O(:,2),O(:,3),NUV(:,1),NUV(:,2),NUV(:,3));
%set(S,'LineStyle','none'); %Remove the thousands of triangle edges that hide the colormap due to the grid density
title('Gradient Computed Unit Normals')
xlim([0 Lx])
ylim([0 Ly])
zlim([-h h])
daspect([1 1 1])
xlabel('x')
ylabel('y')
zlabel('z')  


%Interpolated Velocity Field
figure;
TR=triangulation(C',P');   %Triangulation 
SURFACE=trisurf(TR);             %Render the 3D isosurface
U=U';
O=P';
quiver3(O(:,1),O(:,2),O(:,3),U(:,1),U(:,2),U(:,3));
%set(S,'LineStyle','none'); %Remove the thousands of triangle edges that hide the colormap due to the grid density
title('Velocity Field On The Interface')
xlim([0 Lx])
ylim([0 Ly])
zlim([-h h])
daspect([1 1 1])
xlabel('x')
ylabel('y')
zlabel('z')  


%Interface with Surfactant Concentration Map!
figure;
TR=triangulation(C',P');   %Triangulation (Notice the Transpose of C and P Matrices !!!)
SURFACE=trisurf(TR);             %Render the 3D isosurface
SURFACE.CData=S;
set(SURFACE,'LineStyle','none'); %Remove the thousands of triangle edges that hide the colormap due to the grid density
title('Interface With Surfactant Concentration Rendering')
xlim([0 Lx])
ylim([0 Ly])
zlim([-h h])
daspect([1 1 1])
xlabel('x')
ylabel('y')
zlabel('z')  



end