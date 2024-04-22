function [V,dbox,duples] = dCHEBVAND(deg,X,dbox)

%--------------------------------------------------------------------------
% Object:
% This routine computes the Chebyshev-Vandermonde matrix for degree "deg"
% on a d-dimensional point cloud "X".
%
% The Chebyshev basis is the tensorial Chebyshev basis of total degree
% "deg", shifted on the hyperrectangle defined by "dbox".
%
% If "dbox" is not provided, the routine sets that variable to define the
% smaller "hyperrectangle" (box) with sides parallel to the cartesian
% axes and containing the pointset "X".
%--------------------------------------------------------------------------
% Input:
% deg: polynomial degree;
% X: d-column array of "m" points cloud (matrix "m x d");
% * dbox: variable that defines the smallest hyperectangle with sides
%     parallel to the axis, containing the domain.
%     If "dbox" is not provided, it defines the smaller "hyperrectangle", 
%     with sides parallel to the cartesian axes, containing the pointset 
%     "X". 
%     It is a matrix with dimension "2 x d", where "d" is the dimension  
%     of the space in which it is embedded the domain. 
%     For instance, for a 2-sphere, it is "d=3", while for a 2 dimensional
%     polygon it is "d=2".
%     As example, the set "[-1,1] x [0,1]" is described as "[-1 0; 1 1]".
% 
% Note: the variables with an asterisk "*" are not mandatory and can be 
% also set as empty matrix.
%--------------------------------------------------------------------------
% Output:
% V: shifted Chebyshev-Vandermonde matrix for degree "deg" on the pointset
%    "X", relatively to "dbox".
% dbox: variable that defines the hyperrectangle with sides parallel to the
%     axis, containing the domain.
% duples: polynomial degree in each variable.
%--------------------------------------------------------------------------
% Data:
% The original routine has been written by M. Dessole, F. Marcuzzi and
% M. Vianello on 11/06/2020.
% It has been modified on:
% 22/10/2020 by A. Sommariva;
% 29/10/2020 by M. Vianello;
% 05/11/2020 by A. Sommariva.
%--------------------------------------------------------------------------



% ........................... Function body ...............................



% ...... troubleshooting ......

% box containing the cloud
if nargin < 3, dbox=[]; end
if isempty(dbox) 
    a=min(X); b=max(X); dbox=[a;b];
else
    a=dbox(1,:); b=dbox(2,:);
end



% ..... main code below .....

% d-uples of indices with sum less or equal to "deg" graded lexicographical 
% order
d=size(X,2); 
N = nchoosek(deg+d,d); duples = zeros(N,d);
for i=2:N
    duples(i,:) = mono_next_grlex(d,duples(i-1,:)); 
end

% mapping the mesh in the hypercube "[-1,1]^d"
map = zeros(size(X));
for i=1:d
    map(:,i)=(2*X(:,i)-b(i)-a(i))/(b(i)-a(i)); 
end

% Chebyshev-Vandermonde matrix on the mesh
T=chebpolys(deg,map(:,1));
V=T(:,duples(:,1)+1);
for i=2:d
    T=chebpolys(deg,map(:,i)); 
    V=V.*T(:,duples(:,i)+1); 
end










function T=chebpolys(deg,x)

%--------------------------------------------------------------------------
% Object:
% This routine computes the Chebyshev-Vandermonde matrix on the real line
% by recurrence.
%--------------------------------------------------------------------------
% Input:
% deg: maximum polynomial degree
% x: 1-column array of abscissas
%--------------------------------------------------------------------------
% Output:
% T: Chebyshev-Vandermonde matrix at x, T(i,j+1)=T_j(x_i), j=0,...,deg.
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

T=zeros(length(x),deg+1);
t0=ones(length(x),1); T(:,1)=t0;
t1=x; T(:,2)=t1;

for j=2:deg
    t2=2*x.*t1-t0; 
    T(:,j+1)=t2;
    t0=t1; 
    t1=t2;
end



