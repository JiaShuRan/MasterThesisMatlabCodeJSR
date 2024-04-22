function [coeff,R,jvec,dbox,degs] = dHYPERFIT2(deg,nodes,w,f,jvec,dbox,...
    domain_structure,dimpoly)

%--------------------------------------------------------------------------
% Object:
% This routine computes the coefficients of the weighted-regression
% polynomial for total-degree "deg", w.r.t. a "w"-orthogonal basis on a
% d-dim point cloud "nodes".
%--------------------------------------------------------------------------
% Input: 
% deg: total degree of the algebraic regression polynomial;
% nodes: d-column array of point coordinates;
% w: 1-column array of nonnegative weights, or nonnegative scalar in
%    case of equal weights (in case of doubts set it as "[]");
% f: 1-column array of sampled function values at "nodes";
% * jvec: vector of column indices that selects a polynomial basis;
% * dbox: variable that defines a hyperrectangle with sides parallel to the
%    axis, containing the domain (or pointset "nodes" in the discrete 
%    case).
%    If "dbox" is not provided, it is the smaller "hyperrectangle", with
%    sides parallel to the cartesian axes, containing the pointset "nodes".
%    It is a matrix with dimension "2 x d", where "d" is the dimension of
%    the space in which it is embedded the domain.
%    For instance, for a 2-sphere, it is "d=3", for a 2 dimensional
%    polygon it is "d=2".
%    As example, the set "[-1,1] x [0,1]" is described as
%                          "dbox=[-1 0; 1 1]".
% * domain_structure: structure defining the domain,
%    domain_struct.domain: string with domain name
%    domain_struct.parms: domain parameters.
%
% Note: the variables with an asterisk "*" are not mandatory and can be
% also set as empty matrix.
%--------------------------------------------------------------------------
% Output:
% coeff: 1-column array of weighted regression coefficients;
% R: triangular factor in the economy size QR decomposition
%                   diag(sqrt(w))*C(:,jvec)=Q*R
%   where "C=dCHEBVAND(n,Y)" (or a particular basis if for the domain some 
%   orthogonal basis are available);
% jvec: vector of column indices, selects a polynomial basis;
% dbox: variable that defines the d-dim box where to adapt the
%       product-Chebyshev basis.
%--------------------------------------------------------------------------
% Dates:
% Written on 26/07/2020 by M. Dessole, F. Marcuzzi, M. Vianello.
%
% Modified by:
% 29/10/2020: M. Vianello;
% 03/11/2020: A. Sommariva.
% 04/11/2020: M. Dessole and M. Vianello;
% 05/11/2020: A. Sommariva.
%--------------------------------------------------------------------------

% ..... troubleshooting .....

if nargin < 8, dimpoly=[]; end
if nargin < 7, domain_structure.domain='generic'; end
if isempty(domain_structure), domain_structure.domain='generic'; end
if nargin < 6, dbox=[]; end
if nargin < 5, jvec=[]; end
if nargin < 4, w=[]; end

% .............................  main code ................................

domain=domain_structure.domain;

switch domain
        
    % domain with special basis available
    case {'sphere','circle','unit-simplex','unit-disk'} 
        [V,dbox,degs] = vandermonde_matrix(deg,nodes,dbox,domain_structure);
        coeff=(f.*w)'*V; coeff=coeff';
        R=[]; jvec=1:length(degs);
       
    case {'spherical-triangle','spherical-polygon','spherical-rectangle'}
        dimpoly=(deg+1)^2;
        reorths=1;
        [V,degs]=vandermonde_sphharm(deg,nodes);
        [~,jvec,Q,R,dbox]=dORTHVAND2(deg,nodes,w,jvec,V,dbox,...
            dimpoly,reorths);
        coeff=Q'*(sqrt(w).*f);

    otherwise  % generic case
        [~,jvec,Q,R,dbox,degs]=dORTHVAND2(deg,nodes,w,jvec,[],dbox,dimpoly);
        coeff=Q'*(sqrt(w).*f);

end

