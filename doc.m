%% Spherpinski
%
% Function to compute, display, and save a Sierpinski
% sphere at any iteration number / depth level.
%
% Author & support : nicolas.douillet (at) free.fr, 2017-2020.
%
%% Syntax
%
% Spherpinski;
%
% Spherpinski(nb_it);
%
% Spherpinski(nb_it, option_fill);
%
% Spherpinski(nb_it, option_fill, option_display);
%
% [V, T] = Spherpinski(nb_it, option_fill, option_display);
%
%% Description
%
% Spherpinski computes and display the 3-Sierpinski
% sphere of unit radius, based on the regular octahedron.
%
% Spherpinski(nb_it) computes and display the nb_it Sierpinski
% sphere of unit radius, based on the regular octahedron.
%
% Spherpinski(nb_it, option_fill) fills the octahedron remaining empty
% quadrants with small Sierpinski spherical pieces when option_fill is set to 'double',
% and doesn't -default behaviour- when it is set to 'simple'.
%
% Spherpinski(nb_it, option_fill, option_display) displays it when
% option_display is set to logical *true/1 (default), and doesn't
% when it is set to  logical false/0.
%
% [V,T] = Spherpinski(nb_it, option_fill, option_display) saves the resulting
% vertex coordinates in the array V, and the triangulation in the array T.
%
%% See also
%
% <https://fr.mathworks.com/matlabcentral/fileexchange/73432-n-level-sierpinski-ball Sierpinski_ball> |
% <https://fr.mathworks.com/matlabcentral/fileexchange/79152-sierpinski-octahedron Sierpinski_octahedron>
%
%% Input arguments
%
% - nb_it : positive integer scalar double, the number of iterations / depth level.
%
% - option fill : character string in the set {*'simple,'double','SIMPLE','DOUBLE'}.
%
% - option_display : either logical, *true/false or numeric *1/0.
%
%% Output arguments
%
%        [ |  |  |]
% - V = [Vx Vy Vz], real matrix double, the vertex coordinates. Size(V) = [nb_vertices,3].
%        [ |  |  |]
%
%        [ |  |  |]
% - T = [T1 T2 T3], positive integer matrix double, the triangulation. Size(T) = [nb_triangles,3].
%        [ |  |  |]
%
%% Example #1
% Computes and displays the simple Sierpinski sphere at iteration 3

Spherpinski;

%% Example #2
% Computes, displays, and saves the double Sierpinski sphere at iteration 7

[V,T] = Spherpinski(7,'double',true);