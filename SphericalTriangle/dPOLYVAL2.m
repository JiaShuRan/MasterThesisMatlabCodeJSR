function pval = dPOLYVAL2(deg,coeff,Z,R,jvec,dbox,domain_structure,dimpoly)

%--------------------------------------------------------------------------
% Object:
% This routine evaluates the weighted-regression polynomial of degree
% "deg", defined through its coefficient vector "coeff", at the d-dim
% target points "Z".
%--------------------------------------------------------------------------
% Input:
% deg: polynomial degree
% coeff: 1-column array of of weighted-regression coefficients;
% Z: d-column array of target point coordinates;
% R: invertible triangular factor provided by dPOLYFIT;
% jvec: vector of column indices, selects a polynomial basis;
% dbox: d-dim box where to adapt the product Chebyshev basis.
% dimpoly: dimension of polynomial space.
%--------------------------------------------------------------------------
% Output;
% pval: evaluation of the weighted-regression polynomial at the target
%    points
%--------------------------------------------------------------------------
% Dates:
% Written on 26/07/2020 by M. Dessole, F. Marcuzzi, M. Vianello.
%
% Modified by:
% 29/10/2020: M. Vianello;
% 03/11/2020: A. Sommariva.
% 04/11/2020: M. Dessole and M. Vianello;
% 05/11/2020: A. Sommariva.
% 14/11/2020: A. Sommariva.
%--------------------------------------------------------------------------
% .........................  Function Body ................................

% ..... troubleshooting .....
if nargin < 8, dimpoly=[]; end
if nargin < 7, domain_structure.domain='generic'; end
if nargin < 6, dbox=[]; end
if nargin < 5, jvec=[]; end

% ..... main code below .....
[V,dbox] = vandermonde_matrix(deg,Z,dbox,domain_structure,dimpoly);

domain=domain_structure.domain;

switch domain
    case {'sphere','circle','unit-simplex','unit-disk'}
        pval=V*coeff;
    otherwise
        L=length(R);
        for k=1:L
            jvecL=jvec{k};
            RL=R{k};
            V=V(:,jvecL)/RL;
        end
        pval = V*coeff;
end


