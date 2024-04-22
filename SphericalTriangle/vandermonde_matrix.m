function [V,dbox,degs] = vandermonde_matrix(deg,X,dbox,domain_structure,...
    dimpoly)

%--------------------------------------------------------------------------
% Object:
% This routine computes the Vandermonde matrix for degree "deg"
% on a d-dimensional point cloud "X" as determined by the variable 
% "domain_structure.domain".
% 
% If the domain is unknown, the Chebyshev basis we use the tensorial  
% Chebyshev basis of total degree "deg", shifted on the hyperrectangle 
% defined by "dbox".
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
% * domain_structure: structure defining the domain,
%    domain_struct.domain: string with domain name
%    domain_struct.parms: domain parameters.
% * dimpoly: dimension of polynomial space.
% Note: the variables with an asterisk "*" are not mandatory and can be 
% also set as empty matrix.
%--------------------------------------------------------------------------
% Output:
% V: Vandermonde matrix for degree "deg" on the pointset "X".
% dbox: variable that defines the hyperrectangle with sides parallel to the
%     axis, containing the domain.
%--------------------------------------------------------------------------
% Data:
% The original routine has been written by M. Dessole, F. Marcuzzi and
% M. Vianello on 11/06/2020.
% It has been modified on:
% 22/10/2020 by A. Sommariva;
% 29/10/2020 by M. Vianello;
% 05/11/2020 by A. Sommariva.
% 13/11/2020 by A. Sommariva.
%--------------------------------------------------------------------------



% ........................... Function body ...............................
if nargin < 5, dimpoly=[]; end
if nargin < 4, domain_structure.domain='generic'; end
if isempty(domain_structure), domain_structure.domain='generic'; end

domain=domain_structure.domain;

switch domain
    case {'sphere','spherical-triangle','spherical-polygon',...
            'spherical-rectangle'}
        %fprintf('\n \t * vandermonde_sphharm');
        [V,degs]=vandermonde_sphharm(deg,X);
        dbox=[-1 1; -1 1; -1 1]';
    case 'circle'
        theta=angle(X(:,1)+i*X(:,2));
        V=[ones(size(X,1),1)];
        for k=1:deg, VL=[cos(k*theta); sin(k*theta)]; V=[V; VL]; end
        dbox=[-1 1; -1 1]';
    case 'unit-disk'
        %fprintf('\n \t * vandermonde_logan_shepp');
        [V,degs]=vandermonde_logan_shepp(deg,X);
        % V=vandermonde_zernike(deg,X);
        dbox=[-1 1; -1 1]';
    case 'unit-simplex'
        %fprintf('\n \t * vandermonde_dubiner');
        V=vandermonde_dubiner(deg,X);
        dbox=[0 1; 0 1]';
    otherwise
        %fprintf('\n \t * dCHEBVAND');
        [V,dbox,degs] = dCHEBVAND(deg,X,dbox);
       
end