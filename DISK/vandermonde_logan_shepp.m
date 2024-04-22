
function [V,degs]=vandermonde_logan_shepp(n,pts)

%--------------------------------------------------------------------------
% Object:
% Vandermonde matrix, relatively to the Logan-Shepp basis on the disk, of
% total degree at most "n" and points "pts".
% Notice that all the moments are null with the exception of the first that
% is theoretically "sqrt(pi)".
%--------------------------------------------------------------------------
% Input:
% n: total degree of the basis.
% pts: points in which the basis is considered.
%--------------------------------------------------------------------------
% Output:
% V: Logan-Shepp Vandermonde matrix relatively to the Logan-Shepp basis on
%    the disk, of total degree at most "n" and points "pts".
%    If "dim" is the dimension of the polynomial space and "M" is the
%    cardinality of "pts" then the dimension of "V" is "dim x M".
% degs: vector whose i-th component is degree of the i-th polynomial.
%--------------------------------------------------------------------------

V=[]; degs=[];
for ii=0:n
    Vloc=logan_shepp(pts,(0:ii)',ii); V=[V Vloc];
    degs=[degs; ii*ones(ii+1,1)];
end









function fx=logan_shepp(pts,j,n)

%--------------------------------------------------------------------------
% Input:
% pts: points in which the basis is considered, i.e. pts(k,:)=(x_k,y_k).
% j: M x 1 vector whose components are integers with j(k) <= n
% n: Logan-Shepp degree
%--------------------------------------------------------------------------
% Output:
% fx: evaluation of the j-th Logan-Shepp polynomial of order "n" at "pts".
%     It is a matrix of dimension "card(pts) x (n+1)".
%--------------------------------------------------------------------------

normalized_ls=0; % 0: normalized, 1: not normalized.

x=pts(:,1); y=pts(:,2);

M=size(j,1);
N=size(pts,1);

x=repmat(x,M,1); y=repmat(y,M,1);
j=repmat(j,1,N); j=j'; j=j(:);

args= x.*cos( j*pi/(n+1) ) + y.*sin( j*pi/(n+1) );

if normalized_ls == 0
    fx=(1/sqrt(pi))*chebyshev_second_type(args,n);
else
    fx=chebyshev_second_type(args,n);
end

fx=reshape(fx,N,M);









function fx=chebyshev_second_type(x,n)

%--------------------------------------------------------------------------
% Object:
% Evaluation of Chebyshev polynomials of second-type of degree "n" at "x".
%--------------------------------------------------------------------------
% Input:
% x: points in which the basis is considered.
% n: degree of the Chebyshev polynomials of second-type.
%--------------------------------------------------------------------------
% Output:
% fx: evaluation of Chebyshev polynomials of 2nd-type of degree "n" at "x".
%--------------------------------------------------------------------------

un_2=ones(size(x)); un_1=2*x;

switch n
    case 0
        fx=un_2;
    case 1
        fx=un_1;
    otherwise
        for index=2:n, un=2*x.*un_1-un_2; un_2=un_1; un_1=un; end
        fx=un;
end


