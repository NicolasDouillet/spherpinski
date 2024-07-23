function [V, T] = Spherpinski(nb_it, option_fill, option_display)
%% Spherpinski : function to compute and display the Sierpinski
% sphere at any iteration number / depth level.
%
% Author : nicolas.douillet9 (at) gmail.com, 2017-2024.
%
%
% Syntax
%
% Spherpinski;
% Spherpinski(nb_it);
% Spherpinski(nb_it, option_fill);
% Spherpinski(nb_it, option_fill, option_display);
% [V, T] = Spherpinski(nb_it, option_fill, option_display);
%
%
% Description
%
% Spherpinski computes and display the 3-Sierpinski
% sphere of unit radius, based on the regular octahedron.
%
% Spherpinski(nb_it) computes and display the nb_it Sierpinski
% sphere of unit radius, based on the regular octahedron.
%
% Spherpinski(nb_it, option_fill) fills the sphere remaining empty
% quadrants with small Sierpinski spherical pieces when option_fill is set to 'double',
% and doesn't -default behaviour- when it is set to 'simple'.
%
% Spherpinski(nb_it, option_fill, option_display) displays it when
% option_display is set to logical *true/1 (default), and doesn't
% when it is set to  logical false/0.
%
% [V, T] = Spherpinski(nb_it, option_fill, option_display) saves the resulting
% vertex coordinates in the array V, and the triangulation in the array T.
%
%
% Input arguments
%
% - nb_it : positive integer scalar double, the number of iterations / depth level.
% - option fill : character string in the set {*'simple,'double','SIMPLE','DOUBLE'}.
% - option_display : either logical, *true/false or numeric *1/0.
%
%
% Output arguments
%
%       [ |  |  |]
% - V = [Vx Vy Vz], real matrix double, the vertex coordinates. Size(V) = [nb_vertices,3].
%       [ |  |  |]
%
%       [ |  |  |]
% - T = [T1 T2 T3], positive integer matrix double, the triangulation. Size(T) = [nb_triangles,3].
%       [ |  |  |]
%
%
% Example #1
%
% Computes and displays the simple Sierpinski sphere at iteration 3
%
% Spherpinski;
%
%
% Example #2
%
% Computes, displays, and saves the double Sierpinski sphere at iteration 7
%
% [V,T] = Spherpinski(7,'double',true);


%% Input parsing
assert(nargin < 4,'Too many input arguments.');

if ~nargin
    nb_it = 3;
    option_fill = 'simple';
    option_display = true;
elseif nargin > 0
    assert(isnumeric(nb_it) && nb_it == floor(nb_it) && nb_it >= 0,'nb_it parameter value must be numeric positive or null integer.');
    if nargin > 1
        assert(ischar(option_fill) && (strcmpi(option_fill,'simple') || strcmpi(option_fill,'double')),...
            'option_fill must be a character string in the set {''simple'',''double'',''SIMPLE'',''DOUBLE''}.');
        if nargin > 2
            assert(islogical(option_display) || isnumeric(option_display),'option_display parameter type must be either logical or numeric.');
        else
            option_display = true;
        end
    else
        option_fill = 'simple';
        option_display = true;
    end
end

warning('on');
if option_display && nb_it > 7
    warning('%s facets to display ! Make sure your graphic card has enough memory.',num2str(4^nb_it));
end
warning('off');


%% Body
% Basis vectors
I = [1 0 0]';
J = [0 1 0]';
K = [0 0 1]';

nb_max_it = 9;
sample_step = 2^(nb_max_it-nb_it);

% Create a meshed triangle
[V1, T_array_1] = sample_triangle(I,J,K,sample_step);
[V2, T_array_2] = sample_triangle(0.5*(J+K),0.5*(J-I),0.5*(K-I),sample_step);

% Middle edges vertices
edg_idx1 = 1 + sample_step/2;
edg_idx3_vect = cumsum(edg_idx1:1+sample_step);
edg_idx3 = edg_idx3_vect(end);
edg_idx2 = edg_idx3 - sample_step/2;

V_array_1 = V1;
V_array_2 = V2;

% Iterate over nb_it
p = 0;

while p ~= nb_it
    
    new_V_array_1 = repmat(V_array_1, [1 1 3]);
    new_V_array_2 = repmat(V_array_2, [1 1 3]);
    
    for j = 1:size(V_array_1,3) % Loops on current nb Sierpinski triangles        
            
            new_V_array_1(:,:,3*(j-1)+1) = sample_triangle(V_array_1(1,:,j)',V_array_1(edg_idx1,:,j)',V_array_1(edg_idx2,:,j)',sample_step);
            new_V_array_1(:,:,3*(j-1)+2) = sample_triangle(V_array_1(1 + sample_step,:,j)',V_array_1(edg_idx3,:,j)',V_array_1(edg_idx1,:,j)',sample_step);
            new_V_array_1(:,:,3*(j-1)+3) = sample_triangle(V_array_1(end,:,j)',V_array_1(edg_idx2,:,j)',V_array_1(edg_idx3,:,j)',sample_step);
            
    end
    
    for j = 1:size(V_array_2,3)
            
            new_V_array_2(:,:,3*(j-1)+1) = sample_triangle(V_array_2(1,:,j)',V_array_2(edg_idx1,:,j)',V_array_2(edg_idx2,:,j)',sample_step);
            new_V_array_2(:,:,3*(j-1)+2) = sample_triangle(V_array_2(1 + sample_step,:,j)',V_array_2(edg_idx3,:,j)',V_array_2(edg_idx1,:,j)',sample_step);
            new_V_array_2(:,:,3*(j-1)+3) = sample_triangle(V_array_2(end,:,j)',V_array_2(edg_idx2,:,j)',V_array_2(edg_idx3,:,j)',sample_step);
            
    end
    
    V_array_1 = new_V_array_1;
    V_array_2 = new_V_array_2;
    
    p = p+1;
    
end

V1 = V_array_1(:,:,1);
T1 = T_array_1(:,:,1);

for k = 1:size(V_array_1,3)
    
    T1 = cat(1,T1,T_array_1(:,:,1)+size(V1,1));
    V1 = cat(1,V1,V_array_1(:,:,k));  
    
end

V2 = V_array_2(:,:,1);
T2 = T_array_2(:,:,1);

for k = 2:size(V_array_2,3)
    
    T2 = cat(1,T2,T_array_2(:,:,1)+size(V2,1));
    V2 = cat(1,V2,V_array_2(:,:,k));  
    
end

% Project on the unitary sphere
D1 = sqrt(sum(V1.^2,2));
D2 = sqrt(sum(V2.^2,2));
V1 = V1./D1;

% Simple / double Spherpinski option
if strcmpi(option_fill,'simple')
    
    V = V1;
    T = T1;
    
elseif strcmpi(option_fill,'double')
    
    V2 = V2./D2;
    V = cat(1,V1,V2);
    T = cat(1,T1,T2+size(V1,1));
    
end

% Perform rotations such that the resulting
% Sierpinski sphere is based on the regular octahedron

RzV = ([-1 0 0; 0 -1 0; 0 0 1]*V')';
T = cat(1,T,T+size(V,1));
V = cat(1,V,RzV);

RxV = ([1 0 0; 0 -1 0; 0 0 -1]*V')';
T = cat(1,T,T+size(V,1));
V = cat(1,V,RxV);

% Remove duplicated vertices
[V,T] = remove_duplicated_vertices(V,T);

% % Ellipsoid option
% V(:,2) = 0.25*(1+sqrt(5))^2*V(:,2);
% V(:,1) = 0.5*(1+sqrt(5))*V(:,1);

% Remove duplicated triangles
T = unique(sort(T,2),'rows','stable');

% Display
if option_display
    
    figure;
    set(gcf,'Color',[0 0 0]), set(gca,'Color',[0 0 0]);
    trisurf(T,V(:,1),V(:,2),V(:,3),'EdgeColor',[0 1 0]), shading interp, hold on;
    colormap([0 1 0]);
    axis equal, axis tight, axis off;
    camlight left;
    
end


end % Spherpinski


%% sample_triangle subfunction
function [V, T] = sample_triangle(V1, V2, V3, nbstep)

% Create sampling grid
global Ndim;

Ndim = size(V1,1);

% (V1V2, V1V3) base
u = (V2 - V1);
v = (V3 - V1);

V = zeros(sum(1:nbstep+1),Ndim);

nu = u / norm(u);
nv = v / norm(v);
stepu = norm(u) / nbstep;
stepv = norm(v) / nbstep;
k = 1;

% Sampling & vertices generation
for m = 0:nbstep
    
    for n = 0:nbstep
        
        if m+n <= nbstep % in (V1,V2,V3) triangle conditions ; indices # nb segments
            
            % translation vector
            tv = m*stepu*nu + n*stepv*nv;
            V(k,:) = (V1 + tv)';
            k = k+1;
            
        end
        
    end
    
end

% Index triplets list construction
T = zeros(nbstep^2,3);
row_length = 1 + nbstep;
cum_row_length = row_length;
row_idx = 1;
p = 1;

while p <= nbstep^2 && row_length > 1
    
     i = p;
    
    if p < 2 % "right" triangle serie only
        
        while (i < cum_row_length)
            
            T(row_idx,:) = [i i+1 i+row_length];
            row_idx = row_idx + 1;
            i = i +1;
            
        end
        
        row_length = row_length - 1;
        cum_row_length = cum_row_length + row_length;
        p = p + row_length+1;
        
    else
        
        % Since p >= 2
        while i < cum_row_length % both triangle series
            
            T(row_idx,:) = [i i+1 i+row_length];
            row_idx = row_idx + 1;            
            T(row_idx,:) = [i i-row_length i+1]; % + upside-down triangles serie
            row_idx = row_idx + 1;            
            i = i +1;
            
        end
        
        row_length = row_length - 1;
        cum_row_length = cum_row_length + row_length;
        p = p + row_length+1;
        
    end
    
end

T = sort(T,2);
T = unique(T,'rows','stable');

end % sample_triangle


%% Remove duplicated vertices subfunction
function [V_out, T_out] = remove_duplicated_vertices(V_in, T_in)

tol = 1e4*eps;
[V_out,~,n] = uniquetol(V_in,tol,'ByRows',true);
T_out = n(T_in);

end % remove_duplicated_vertices