

function xyw = cub_circsect(n,omega,r1,r2)

%--------------------------------------------------------------------------
% Object:
% The routine computes the nodes and weights of a product gaussian
% formula on a circular annular sector centered at the origin
% with angles in [-omega,omega]
%--------------------------------------------------------------------------
% Input:
% n: algebraic degree of exactness
% omega: half-length of the angular interval, 0<omega<=pi
% r1,r2: internal and external radius, 0<=r1<r2
%--------------------------------------------------------------------------
% Output:
% xyw: (ceil((n+2)/2) x (n+1)) x 3 array of (xnodes,ynodes,weights)
%--------------------------------------------------------------------------
% Required routines:
% 1. r_jacobi.m (www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html)
% 2. gauss.m (www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html)
% 3. quad_trig.m 
%--------------------------------------------------------------------------
% Written by Gaspare Da Fies and Marco Vianello, University of Padova
% Date: November 8, 2011.
% Last update: January 4. 2020.
%--------------------------------------------------------------------------

% trigonometric gaussian formula on the arc
tw=quad_trig(n,-omega,omega);

% algebraic gaussian formula on the radial segments
ab=r_jacobi(ceil((n+2)/2),0,0);
xw=gauss(ceil((n+2)/2),ab);
xw(:,1)=xw(:,1)*(r2-r1)/2+(r2+r1)/2;
xw(:,2)=xw(:,2)*(r2-r1)/2;

% creating the polar grid
[r,theta]=meshgrid(xw(:,1),tw(:,1));
[w1,w2]=meshgrid(xw(:,2),tw(:,2));

% nodal cartesian coordinates and weights
xyw(:,1)=r(:).*cos(theta(:));
xyw(:,2)=r(:).*sin(theta(:));
xyw(:,3)=r(:).*w1(:).*w2(:);



