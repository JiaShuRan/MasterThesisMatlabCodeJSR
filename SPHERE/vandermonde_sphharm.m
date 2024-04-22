
function [V,degs]=vandermonde_sphharm(deg,nodes)

%--------------------------------------------------------------------------
% Object:
% This routine computes the spherical harmonics Vandermonde matrix of 
% degree "deg" on "nodes", whose cartesian coordinates are stored in 
% "nodes".
%--------------------------------------------------------------------------
% Input:
% deg: degree of the Vandermonde matrix
% nodes: matrix N x 3, describing N points of the sphere via cartesian
%      coordinates.
%--------------------------------------------------------------------------
% Output:
% V: Vandermonde matrix of degree "deg" on "nodes", whose cartesian  
%  coordinates are stored in "nodes". If "N" is the cardinality of the
%  nodes, then V has dimension "N x ((deg+1)^2)".
%--------------------------------------------------------------------------


z_phi=cart_2_cyl(nodes); z=z_phi(:,1); phi=z_phi(:,2);
[V,degs]=vandermonde_sph_harm(deg,z,phi);



function [VAN,degs]=vandermonde_sph_harm(deg,z,phi)

% create a vandermonde matrix of spherical harmonics at the points (z,phi)
% where we are using the cylindrical coordinates on the sphere that is
% x=cos(phi)*sqrt(1-z^2)
% y=sin(phi)*sqrt(1-z^2)
% z=z
%
% programs not usually lie in standard Matlab versions
% spharmrd.m
% from Paul Leopardi
% plkd.m (in spharmrd)
% from Paul Leopardi

cont=0;
degs=[];
for i=0:1:deg
    [VAN(cont+1:cont+2*i+1,:)]=spharmrd(i,z,phi);
    cont=cont+2*i+1;
    degsL=i*ones(1,2*i+1);
    degs=[degs degsL];
end
VAN=VAN';


function [YL,DYLdt,DYLdp] = spharmrd(L,z,phi)

%   Return a matrix YL containing the real spherical harmonics
%   Y_{L,-L}:Y_{L,L} of degree L at the points (z, phi),
%   with (z,phi) representing points on the sphere S^2 in cylindrical
%   coordinates.
%   Optionally return matrices DYLdt and DYLdp containing
%   partial derivatives with respect to theta and phi
%   where z == cos(theta).
%
% -- Input arguments --
% L      : Degree, Integer >= 0
% z      : Row vector of altitudes of points where functions are to be evaluated
% phi    : Row vector of longitudes of points where functions are to be evaluated
% -- Output arguments --
% YL     : Matrix of size (2*L+1)*prod(size(z)) giving the function values at z
% DYLdt  : Matrix of size (2*L+1)*prod(size(z)) giving theta derivatives at z
% DYLdp  : Matrix of size (2*L+1)*prod(size(z)) giving phi derivatives at z
%        : If called with only one output argument the derivatives will not be calculated
%
%   Layout of the matrix YL is:
%   YL(1,:)   = Y_{L,-L}(z,phi) = scaled P_L^L(z)*sin(L*phi)
%   ...
%   YL(L,:)   = Y_{L,-1}(z,phi) = scaled P_L^1(z)*sin(phi)
%   YL(L+1,:) = Y_{L,0} (z,phi) = scaled P_L(z)
%   YL(L+2,:) = Y_{L,1} (z,phi) = scaled P_L^1(z)*cos(phi)
%   ...
%   YL(N,:)   = Y_{L,L} (z,phi) = scaled P_L^1(z)*cos(L*phi)
%
%   The matrices DYLdt and DYLdp have appropriate partial
%   derivatives in the corresponding positions.
%
%   See http://mathworld.wolfram.com/SphericalHarmonic.html
%   for definition.
%   See also http://mathworld.wolfram.com/LegendrePolynomial.html
%   References: 
%   A,A&R: Andrews, Askey & Roy, "Special Functions",
%   Cambridge, 1999, pp302-303, 456-457.
%   A&S: Abramowitz & Stegun, "Handbook of Mathematical Functions",
%   pp333-334,778
%
%   [YL,DYLdt,DYLdp] = spharmrd(L,z,phi)

%   $Revision: 2.2  $ Paul Leopardi 2003-08-22
%   Improve comments and tidy up.
%   $Revision: 2.1  $ Paul Leopardi 2003-08-21
%   Factor PL calculation back out into plkd
%   $Revision: 2.0  $ Paul Leopardi 2003-08-20
%   Thoroughly debug
%   $Revision: 1.5  $ Paul Leopardi 2002-08-02
%   Incorporate diff_PL
%   $Revision: 1.4  $ Paul Leopardi 2002-08-01
%   Optimized spharmrd based on spharmvr and diff_PL
%   $Revision: 1.3  $ Paul Leopardi 2002-07-11
%   Optimized spharmvr based on spharmvz
%   $Revision: 1.2  $ Paul Leopardi 2002-07-10
%   Optimized spharmpl based on spharmx
%   $Revision: 1.1  $ Paul Leopardi 2002-07-09
%   Harmonize with spharmxc
%   $Revision: 1.0  $ Paul Leopardi 2002-07-02
%   Adapted by Paul Leopardi 2002-07-02
%   for UNSW School of Mathematics
%   Based on SPHARM2:
%   Denise L. Chen  9-1-93
%   Copyright 1984-2000 The MathWorks, Inc.
%   $Revision: 5.10 $  $Date: 2000/06/01 03:46:27 $

if nargin == 0
  help spharmrd
  return;
end

CalcDt = nargout > 1;
CalcDp = nargout > 2;

% Ensure input consists of row vectors
%
if size(z,2) == 1
    z = z';
end
if size(phi,2) == 1
    phi = phi';
end

% Get associated Legendre functions and derivatives.
%
if CalcDt
    [PL,DPLdt] = plkd(L,z);
else
    PL = plkd(L,z);
end

k = 1:L;
N = 2*L+1;
lz = prod(size(z));

YL = zeros(N,lz);
scaling = sqrt(N/(4*pi));
YL(1+L,:)        = scaling * PL(1,:);
if CalcDt
    DYLdt = zeros(N,lz);
    DYLdt(1+L,:) = scaling * DPLdt(1,:);
end
if CalcDp
    DYLdp = zeros(N,lz);
end
if L > 0
    cp = (cumprod(L:2*L).*[1,1,cumprod(L-1:-1:1)])';
    scaling = diag(sqrt(N./(2*pi*cp(k+1))));
    scaledPL  = scaling * PL(1+k,:);
    kkp = k'*phi;

    % sin comes first, in reverse order
    %
    sk = L:-1:1;

    % cos comes after P_L, in forward order
    %
    ck = L+2:N;

    YL(sk,:)        = scaledPL  .* sin(kkp);
    YL(ck,:)        = scaledPL  .* cos(kkp);
    
    % If required, use Leibniz rul to calculate derivatives
    %
    if CalcDt
        scaledDPLdt = scaling * DPLdt(1+k,:);
        DYLdt(sk,:) = scaledDPLdt .* sin(kkp);
        DYLdt(ck,:) = scaledDPLdt .* cos(kkp);
    end
    if CalcDp
        DYLdp(sk,:) = scaledPL  .* ( diag(k) * cos(kkp));
        DYLdp(ck,:) = scaledPL  .* (-diag(k) * sin(kkp));
    end
end










function [PL,DPLdt] = plkd(L,z)

%   Return a matrix PL containing the associated Legendre functions
%   P_L^0:P_L^L of degree L at the points z
%   Optionally return a matrix DPLdt containing the
%   derivative with respect to theta where z == cos(theta).
%
% -- Input arguments --
% L      : Degree, Integer >= 0
% z      : Row vector of points in [-1, 1] where functions are to be evaluated
% -- Output arguments --
% PL     : Matrix of size (L+1)*prod(size(z)) giving the function values at z
% DPLdt  : Matrix of size (L+1)*prod(size(z)) giving derivatives at z
%        : If called with only one output argument the derivative will not be calculated
%
%   See http://mathworld.wolfram.com/LegendrePolynomial.html
%   References:
%   A&S: Abramowitz & Stegun, "Handbook of Mathematical Functions".
%   PL: Paul Leopardi, "Computing derivatives of associated Legendre functions using Matlab".
%   S: Szego, "Orthogonal Polynomials".
%
%   [PL,DPLdt] = plkd(L,z)

%   $Revision: 1.2  $ Paul Leopardi 2003-08-26
%   Clean up comments
%   $Revision: 1.1  $ Paul Leopardi 2003-08-22
%   Simplify expressions and remove the calls to ultra
%   $Revision: 1.0  $ Paul Leopardi 2003-08-20
%   Based on spharmrd
%   for UNSW School of Mathematics
%   Based on SPHARM2:
%   Denise L. Chen  9-1-93
%   Copyright 1984-2000 The MathWorks, Inc.
%   $Revision: 5.10 $  $Date: 2000/06/01 03:46:27 $

if nargin == 0
  help plkd
  return;
end

CalcDt = nargout > 1;
if size(z,2) == 1
    z = z';
end
k = 1:L;
lz = prod(size(z));
PL = legendre(L,z);
if CalcDt
    DPLdt = zeros(L+1,lz);
    if L > 0
        pole = abs(z) == 1;
        if ~isempty(z(pole))

            % The following special cases for the derivatives at -1 and at 1
            % are defined in terms of limits. 
            % The formula for order 1 is PL Theorem 4.
            % It can be derived from S 4.21.7, p63, S 4.1.1, p58 and S 4.1.3, p59.
            % For all other orders, the derivative is 0, as per PL Theorems 3 and 5.
            %
            DPLdt(2,pole) = -L*(L+1)/2*PL(1,pole);
        end
        if ~isempty(z(~pole))
            oz = ones(size(z));
            ok = ones(size(k))';
            sint = sqrt(1-z.*z);

            % The following general formula for the derivative is PL Theorem 1.
            % It can be derived from A&S 8.5.2, p334, assuming -1 < z < 1.
            %
            DPLdt(1+k,~pole) = (((k-L-1).*(k+L))'*oz(~pole)).*PL(k,~pole) - ...
                (k'*oz(~pole)).*PL(1+k,~pole).*(ok*(z(~pole)./sint(~pole)));

            % The following special case is PL Theorem 2.
            %
            DPLdt(1,~pole) = PL(2,~pole);
        end
    end
end





function z_phi=cart_2_cyl(nodes)

%--------------------------------------------------------------------------
% PURPOSE:
%-------------
% This routine mimics Matlab built-in routine "cart2sph" but
%
% theta: lies in [0,pi] x [0,2*pi] (and not in [-pi/2,pi/2]) x [-pi,pi]).
%--------------------------------------------------------------------------
% INPUT:
%-------------
%
% nodes: n x 3 MATRIX, CONTAINING THE CARTESIAN COORDINATES OF THE POINTS.
%
%--------------------------------------------------------------------------
% OUTPUT:
%-------------
%
% z_phi: n x 3 MATRIX, CONTAINING THE CYLINDRICAL COORDINATES OF THE POINTS.
%
%--------------------------------------------------------------------------
% NOTE:
%-------------
%
% We used the following spherical coordinates:
%
% x_1=cos(theta_2)*sin(theta_1)
% x_2=sin(theta_2)*sin(theta_1)
% x_3=cos(theta_1)
%
% where theta_1=theta(1), theta_2=theta(2).
%
%--------------------------------------------------------------------------
% EXAMPLE 1:
%-------------
%
% >> x =[6.651394320702556e-01  7.394874417699735e-01 1.036718832169940e-01];
% >> theta=cart_2_sph(x);
% >> xx=sph_2_cart(theta)
%
% xx =
%
%      6.651394320702555e-01     7.394874417699734e-01     1.036718832169940e-01
%
% >> x =[6.651394320702556e-01  7.394874417699735e-01 1.036718832169940e-01];
% >> z_phi=cart_2_cyl(x)
%
% z_phi =
%
%      1.036718832169940e-01     8.382796045410373e-01
%
% >> theta=[acos(1.036718832169940e-01) 8.382796045410373e-01];
% >> xx=sph_2_cart(theta)
%
% xx =
%
%      6.651394320702555e-01     7.394874417699734e-01     1.036718832169940e-01
%
% >>
%
%--------------------------------------------------------------------------
% EXAMPLE 2:
%-------------
%
% >>z_phi=cart_2_cyl([1 0 0])
%
% z_phi =
%
%      0     0
%
% >>
%
%--------------------------------------------------------------------------
% EXAMPLE 3:
%-------------
%
% >> z_phi=cart_2_cyl([0 0 1])
%
% z_phi =
%
%      1     0
%
% >>
%
%--------------------------------------------------------------------------
%%
%% Copyright (C) 2013. Alvise Sommariva, Marco Vianello.
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
%%
%%          Alvise Sommariva <alvise@euler.math.unipd.it>
%%          Marco Vianello   <marcov@euler.math.unipd.it>
%%
%% Date: August 9, 2016
%%
%-------------------------------------------------------------------------

theta=cart2sph_sc_ss_c(nodes);
z_phi=[nodes(:,3) theta(:,2)];



function theta=cart2sph_sc_ss_c(nodes)

%--------------------------------------------------------------------------
% PURPOSE:
%-------------
% This routine mimics Matlab built-in routine "cart2sph" but
%
% theta: lies in [0,pi] x [0,2*pi] (and not in [-pi/2,pi/2]) x [-pi,pi]).
%--------------------------------------------------------------------------
% INPUT:
%-------------
%
% nodes: n x 3 MATRIX, CONTAINING THE CARTESIAN COORDINATES OF THE POINTS.
%        POINTS MUST BE STORED AS ROW VECTORS.
%
%--------------------------------------------------------------------------
% OUTPUT:
%-------------
%
% theta: n x 3 MATRIX, CONTAINING THE SPHERICAL COORDINATES OF THE POINTS.
%
%--------------------------------------------------------------------------
% NOTE:
%-------------
%
% We used the following spherical coordinates:
%
% x_1=sin(theta_1)*cos(theta_2)
% x_2=sin(theta_1)*sin(theta_2)
% x_3=cos(theta_1)
%
% where theta_1=theta(1), theta_2=theta(2).
%
%--------------------------------------------------------------------------
% EXAMPLES:
%-------------
% >> % EXAMPLE 1.
% >> theta=cart2sph_sc_ss_c([1 0 0])
%
% theta =
%
%      1.570796326794897e+00                   0
%
% >> pts=sph2cart_sc_ss_c(theta)
%
% pts =
%
%      1.000000000000000e+00                   0     6.123233995736766e-17
%
% >> % EXAMPLE 2.
% >> x =[6.651394320702556e-01 7.394874417699735e-01 1.036718832169940e-01];
% >> theta=cart2sph_sc_ss_c(x)
%
% theta =
%
%      1.466937831133291e+00   8.382796045410373e-01
%
% >> pts=sph2cart_sc_ss_c(theta)
%
% pts =
%
%      6.651394320702555e-01  7.394874417699734e-01  1.036718832169940e-01
%
%--------------------------------------------------------------------------
%%
%% Copyright (C) 2013. Alvise Sommariva, Marco Vianello.
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
%%
%%          Alvise Sommariva <alvise@euler.math.unipd.it>
%%          Marco Vianello   <marcov@euler.math.unipd.it>
%%
%% Date: May 19, 2013
%%
%-------------------------------------------------------------------------

[x_theta,x_phi,r]=cart2sph(nodes(:,1),nodes(:,2),nodes(:,3));

for i=1:1:size(x_theta)
    if x_theta(i) < 0
        x_theta(i)=2*pi+x_theta(i);
    end
end

for i=1:1:size(x_phi)
    if x_phi(i)>=0
        x_phi(i)=-(x_phi(i)-pi/2);
    end
    if x_phi(i)<0
        x_phi(i)=(-x_phi(i)+pi/2);
    end
end

theta=[x_phi x_theta];
