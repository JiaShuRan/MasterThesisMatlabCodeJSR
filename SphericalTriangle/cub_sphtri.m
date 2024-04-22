
function xyzw = cub_sphtri(n,P1,P2,P3,pos)

%--------------------------------------------------------------------------
% Object:
% The routine computes a cubature formula "xyzw" with (numerical) algebraic 
% degree of precision "n" on the spherical triangle with vertices "P1",  
% "P2", "P3", by a 2D formula of degree of precision "m+n" on the xy  
% projection of the spherical triangle rotated to put the barycenter at the
% north pole.
%
% The parameter "m" depends on the size of the circumradius and it is "2" 
% for "small" "r", becoming bigger as "r" approaches "1".
%
% IMPORTANT: This routine requires Chebfun.
%--------------------------------------------------------------------------
% Input: 
% n: algebraic degree of precision of the rule;
% P1,P2,P3: column arrays of the spherical triangle vertices coords
% pos: positive weights for pos=1, possible neg weights (faster) for pos=0
%--------------------------------------------------------------------------
% Output:
% xywz : 4-column array of nodes cartesian coords and weights
%--------------------------------------------------------------------------
% Routines called:
% 1. compute_m (attached)
% 2. cub_circsect (external)
% 3. dCATCH (external)
%--------------------------------------------------------------------------
% Data:
% The original routine has been written by M. Vianello on May 2019.
% Last Update: 01/01/2021 by A. Sommariva.
%--------------------------------------------------------------------------
%% Copyright (C) 2021
%% Alvise Sommariva, Marco Vianello.
%%
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 2 of the License, or
%% (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
%%
%% Authors: 
%% Alvise Sommariva, Marco Vianello.
%%
%% Date: January 01, 2021
%--------------------------------------------------------------------------


% ............................ troubleshooting ............................

% working with vertices as column vectors
if size(P1,1) == 1, P1=P1'; end
if size(P2,1) == 1, P2=P2'; end
if size(P3,1) == 1, P3=P3'; end

% ------------------------------- main code -------------------------------
% We scale the problem in the unit-sphere with center [0 0 0].
% In other words, we make up a cubature rule on the unit sphere (by
% contraction w.r.t. the unit sphere and map back the so obtained rule to
% the starting sphere with radius R).

% ... barycenter ...
RR=norm(P1);
vert0=[P1 P2 P3]/RR; 
CC=1/3*(P1+P2+P3); CC=CC/norm(CC);

% .................... rotation matrix to the north pole ..................

[az,el,r] = cart2sph(CC(1),CC(2),CC(3));
phi=az; theta=pi/2-el;
cp=cos(phi); sp=sin(phi); ct=cos(theta); st=sin(theta);
R1=[ct 0 -st; 0 1 0; st 0 ct]; R2=[cp sp 0; -sp cp 0; 0 0 1]; 
R=R1*R2; invR=R';

% ............... vertices of the triangle at the North Pole ..............

vert1=R*vert0; 

% ...... computing "m" needed in determining the degree of precision ......
m=compute_m(vert1');

% .................... determining quadrature rule ........................  

xyzw=[]; % nodes on the sphere
vert1=[vert1 vert1(:,1)];

for i=1:3
    
    % affine transformation matrix 
    P=vert1(:,i); Q=vert1(:,i+1);
    om=acos(P'*Q);
    xi=cos(om/2); eta=sin(om/2);
    M=[xi eta 0 0; 0 0 xi eta; xi -eta 0 0; 0 0 xi -eta];
    h=[P(1); P(2); Q(1); Q(2)];
    u=M\h;
    T=[u(1) u(2); u(3) u(4)];
    
    % nodes and weights for the xy-projection of the rotated spher.triangle
    nw=cub_circsect(m+n,om/2,0,1); nod2=T*nw(:,1:2)';
    
    % spherical triangle nodes and weights on North Pole spher.triangle
    x=nod2(1,:); y=nod2(2,:); z=sqrt(1-x.^2-y.^2); 
    
    weights=abs(det(T))*nw(:,3)./z';
    
    % inverse rotation to original spherical triangle
    nodes=invR*[x; y; z];
    
    % cubature rule update
    xyzw=[xyzw;[nodes' weights]];
    
end

% Exporting results to the sphere with radius "R".
X=RR*xyzw(:,1); Y=RR*xyzw(:,2); Z=RR*xyzw(:,3); W=RR*xyzw(:,4);
xyzw=[X Y Z W];






function m=compute_m(vert)

%--------------------------------------------------------------------------
% Object:
% It computes the "m" positive integer value depending from the sph. tri.
% with vertices "vert", so that a WAM over the sph. triangle with degree
% equal to "n" can be obtained via WAM of degree "m+n" on its projection on
% the xy-plane.
%--------------------------------------------------------------------------
% Input:
% vert: points defining the sph. triangle (the k-th point is described 
%    by the k-th row.
%--------------------------------------------------------------------------
% Output:
% m: positive integer value depending from the sph.triangle with vertices 
% "vert", so that a WAM over the sph. triangle with degree equal to "n" can
% be obtained via cubature of degree "m+n" on its projection on the 
% xy-plane.
%--------------------------------------------------------------------------
% Data:
% First version: 13/01/2021 by A. Sommariva and M. Vianello.
%--------------------------------------------------------------------------

r=sqrt( (vert(:,1)).^2 + (vert(:,2)).^2 ); 
r0=max(r);

intvf=[0,r0]; F=@(r) sqrt(1-r); f=chebfun(F,intvf);
m=max(1,(length(f)-1));








