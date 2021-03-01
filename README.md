# Marching_Interface

A tool for the extraction of isosurfaces from a scalar 3D field on rectirilinear grids 

This code is a tool that is used for the extraction of the interface between two fluid phases in a 3D channel flow. The kernel of the code is based on the marching cubes algorithm. The main idea lies in the construction of an interface as a contour of a scalar field based on a preselected isolevel. In the examined case the scalar field is the phase field and the isolevel is zero which is the point where according to phase field method is where the first phase ends and the second begins. A 3D scanning of the domain takes place looking at the values of the phase field at each grid cube (note that in the case of a non-uniform Chebyshev descritization we are dealing with parallelograms instead of cubes, however the algorithm still performs for adequately fine descritazations) and all the locations where the field obtains a value equal to the isolevel are determined by means of linear interpolation. Based on the connectivity of these points a triangulation takes place and the resulting 3D shape representing the interface is rendered. If avalaible the velocity and surfactant fields can also be scanned and interpolated on the interface.
Extensions of the main output include the determination of unit normal vectors, based on both the estimation of cross-products between triangle edges and the interpolation of the phase field divergence on the interface, the creation of surfactant colormaps and the projection of velocity vectors on the interface. The computational time rises with every one of these extensions. 


Subroutines:

main.m: grid specifications size and x,y,z coordinate files reading. Reading phase field and also velocity and surfactant fields. A call to marching_cubes.m. The display of the overall time and the total interface area after the return of the subroutine. A call to figures.m where all the plots are created. 

marching_cubes.m: a triple for loop where working on a single cube (parallelogram) per step the points on the edges where	the interface passes from are retrieved through a call to the interpolator.m. The connectivity of these points is then figured by a call to a triangulation matrix. In this way the different triangles which will represent the interfacece can later be drawn.

interpolator.m: based on the values of the phase field on the vertexes of each cube it estimates the location of the intersection of the cube edge with the interface 

triangulation.m: the triangulation matrix

figures.m: all plots

An input folder should be created inside the folder where the .m files are. This folder should be named "input" and it should contain the fields that are going to be read, namely: phi.dat, u.dat, v.dat, w.dat and also the grid coordinate arrays x.dat, y.dat and z.dat.
