
function XYZW=cub_cube(ade,fam)

%--------------------------------------------------------------------------
% OBJECT.
%--------------------------------------------------------------------------
% Non-tensor Clenshaw-Curtis like formula for the Chebyshev measure 
%
%           w(x,y,z)=(1/pi^3)*sqrt(1-x^2)*sqrt(1-y^2)*sqrt(1-z^2)
% 
% in the unit-cube [-1,1]^3.
%--------------------------------------------------------------------------
% INPUT
%--------------------------------------------------------------------------
% ade: requested alg. degree of precision of the rule for integration of 
%               int(f(x,y,z)*w(x,y,z))    x,y,z in [-1,1]^3;
% 
% fam: parameter that allows to choose one of the 4 sub-grids below
%     (1 by default).
%--------------------------------------------------------------------------
% OUTPUT
%--------------------------------------------------------------------------
% nodes: cubature nodes;
% weights: cubature weights.
%--------------------------------------------------------------------------
% AUTHORS
%--------------------------------------------------------------------------
% Stefano De Marchi, Marco Vianello
% University of Padova
%--------------------------------------------------------------------------
% PAPER
%--------------------------------------------------------------------------
% NEW CUBATURE FORMULAE AND HYPERINTERPOLATION IN THREE VARIABLES
% S. DE MARCHI, M. VIANELLO and Y. XU
% BIT Numerical Mathematics (2006), DOI:10.1007/s10543-000-0000-x
%--------------------------------------------------------------------------
% First version: May 2008
% Revision: August 2022
% Modified by Alvise Sommariva.
%--------------------------------------------------------------------------

if nargin < 1, ade=5; end
if nargin < 2, fam=1; end

% .......................... Nodes generation .............................

n=floor(ade/2); % in the paper setting, the formula had ade="2*n+1".

k=0:n+1; c=cos(k*pi/(n+1));   % Chebyshev-Lobatto points
[X,Y,Z]=ndgrid(c);     % 3d grid of Chebyshev points in the cube
Sc=[X(:), Y(:), Z(:)]; % Chebyshev grid arranged as array with 3 coord.

% Here we extract the sub-vector E=EVEN and O=ODD of the Chebsyhev points
E=c(1:2:length(c)); O=c(2:2:length(c));

%--------------------------------------------------------------------------
% Generation of the nodes sub-grids.
%
% We have 4 possible families:
% 1) EEE OOO
% 2) EOO OEE
% 3) OEE EOO
% 4) EOE OEO
%--------------------------------------------------------------------------

switch(fam)
    case 1 % Family choosen in the lasso paper (default).
        [X1,Y1,Z1]=meshgrid(E,E,E);
        [X2,Y2,Z2]=meshgrid(O,O,O);
    case 2
        [X1,Y1,Z1]=meshgrid(E,E,O);
        [X2,Y2,Z2]=meshgrid(O,O,E);

    case 3
        [X1,Y1,Z1]=meshgrid(E,O,O);
        [X2,Y2,Z2]=meshgrid(O,E,E);
    case 4
        [X1,Y1,Z1]=meshgrid(O,O,O);
        [X2,Y2,Z2]=meshgrid(E,E,E);
end

EEE=[X1(:),Y1(:),Z1(:)]; OOO=[X2(:),Y2(:),Z2(:)]; XYZ=[EEE; OOO];

% .......................... Weights generation ...........................

X=XYZ(:,1); Y=XYZ(:,2); Z=XYZ(:,3);

Xtest= (X == -1) | (X == 1);
Ytest= (Y == -1) | (Y == 1);
Ztest= (Z == -1) | (Z == 1);

Wfactor=2.^(-(Xtest+Ytest+Ztest));
W0=(4/(n+1)^3)*ones(size(X));
W=Wfactor.*W0;

% .......................... Nodes/weights storage ........................

XYZW=[XYZ W];
