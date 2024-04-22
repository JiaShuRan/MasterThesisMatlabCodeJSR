function XW=quad_gausslegendre(n)

%--------------------------------------------------------------------------
% Object:
% Gauss Legendre rule on the unit interval [-1,1].
% * Important: ade=2*n-1
%--------------------------------------------------------------------------
% Input:
% n: cardinality of the Gaussian rule
%--------------------------------------------------------------------------
% Output:
% XW: n x 2 matrix, where XW(:,1) are the nodes and XW(:,2) are the
%    weights.
%--------------------------------------------------------------------------

% Gaussian rule in [-1,1]: ade=2*n -> n+1 points.
ab=r_jacobi(n+1,0,0); XW=gauss(n+1,ab);