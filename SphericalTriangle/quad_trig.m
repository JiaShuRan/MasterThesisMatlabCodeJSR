
function tw=quad_trig(n,alpha,beta)

%--------------------------------------------------------------------------
% Object:
% Computation of trigonometric gaussian rules on the unit arc from "alpha"
% to "beta".
%--------------------------------------------------------------------------
% Inputs:
% n    : the rule computes computes the n+1 angles and weights of a 
%        trigonometric gaussian quadrature formula on [alpha,beta], with 
%                          0 < beta-alpha <= 2*pi;
% alpha, beta: arc angles for trigonometric gaussian rules on the unit arc
%       from "alpha" to "beta".
%--------------------------------------------------------------------------
% Outputs
% tw   : (n+1) x 2 matrix, where the first column contains the nodes, while
%        the second one contains the weights.
%--------------------------------------------------------------------------
% First version: May 18, 2013 by G. Da Fies, A. Sommariva, M. Vianello.
% Successive versions have been made by G. Meurant, A. Sommariva and
% M. Vianello.
% Last release: January 04, 2020.
%--------------------------------------------------------------------------

% .......................... troubleshooting  .............................

if nargin < 1, n=10; end
if nargin < 2, beta=pi; alpha=-beta; end
if nargin < 3, beta=alpha; alpha=-beta; end

n=n+1;
omega=(beta-alpha)/2;
ab = r_subchebyshev(n,omega);
xw_symm_eigw = SymmMw(n,ab);
tw=quad_trig_conversion(xw_symm_eigw,omega);
tw(:,1)=tw(:,1)+(beta+alpha)/2;









function ab=r_subchebyshev(n,omega)

%--------------------------------------------------------------------------
% Object:
% Recurrence coeffs for the monic OPS w.r.t. the weight function
%         w(x)=2*sin(omega/2)/sqrt(1-sin^2(omega/2)*x^2)
% by the modified Chebyshev algorithm.
% The reference angle of the rule is [-omega,omega].
%--------------------------------------------------------------------------
% Inputs:
% n     : number of points.
% omega : arc angle.
%--------------------------------------------------------------------------
% Output:
% ab   : three terms recursion
%--------------------------------------------------------------------------
% First version: May 18, 2013 by G. Da Fies, A. Sommariva, M. Vianello.
% Successive versions have been made by G. Meurant, A. Sommariva and
% M. Vianello.
% Last release: January 04, 2020.
%--------------------------------------------------------------------------

N=n; n=n-1;

% modified Chebyshev moments by recurrence
if rem(N,2) == 1, NN=N+1; nn=n+1; else NN=N; nn=n; end
mom=fast_moments_computation(omega,2*nn+1);

% recurrence coeffs of the monic Chebyshev polynomials
abm(:,1)=zeros(2*nn+1,1);
abm(:,2)=0.25*ones(2*nn+1,1); abm(1,2)=pi; abm(2,2)=0.5;

% recurrence coeffs for the monic OPS w.r.t. the weight function
ab = fast_chebyshev(NN,mom,abm);









function x = tridisolve(a,b,c,d)

%--------------------------------------------------------------------------
% Object:
% Solution of tridiagonal system of equations.
%--------------------------------------------------------------------------
% From Cleve Moler's Matlab suite
% http://www.mathworks.it/moler/ncmfilelist.html
%--------------------------------------------------------------------------
% x = TRIDISOLVE(a,b,c,d) solves the system of linear equations
%
%     b(1)*x(1) + c(1)*x(2) = d(1),
%     a(j-1)*x(j-1) + b(j)*x(j) + c(j)*x(j+1) = d(j), j = 2:n-1,
%     a(n-1)*x(n-1) + b(n)*x(n) = d(n).
%
% The algorithm does not use pivoting, so the results might be inaccurate 
% if abs(b) is much smaller than abs(a)+abs(c).
%
% More robust, but slower, alternatives with pivoting are:
%     x = T\d where T = diag(a,-1) + diag(b,0) + diag(c,1)
%     x = S\d where S = spdiags([[a; 0] b [0; c]],[-1 0 1],n,n)
%--------------------------------------------------------------------------
% Optimized version by G. Meurant.
%--------------------------------------------------------------------------

x = d;
n = length(x);
bi = zeros(n,1);

for j = 1:n-1
    bi(j) = 1 / b(j);
    mu = a(j) * bi(j);
    b(j+1) = b(j+1) - mu * c(j);
    x(j+1) = x(j+1) - mu * x(j);
end

x(n) = x(n) / b(n);
for j = n-1:-1:1
    x(j) = (x(j) - c(j) * x(j+1)) * bi(j);
end






function ab=fast_chebyshev(N,mom,abm)

%--------------------------------------------------------------------------
% Object:
% Modified Chebyshev algorithm, that works only for the subperiodic weight 
% function.
%--------------------------------------------------------------------------
% From Gautschi's code (simplified)
% Mar 2012
%--------------------------------------------------------------------------
% Optimized version by G. Meurant.
%--------------------------------------------------------------------------

ab = zeros(N,2);
sig = zeros(N+1,2*N);

ab(1,2) = mom(1);

sig(1,1:2*N) = 0;
sig(2,:) = mom(1:2*N);

for n = 3:N+1
    for m = n-1:2*N-n+2
        sig(n,m) = sig(n-1,m+1) + abm(m,2) * sig(n-1,m-1) - ...
            ab(n-2,2) * sig(n-2,m);
    end
    
    ab(n-1,2) = sig(n,n-1) / sig(n-1,n-2);
end









function mom=fast_moments_computation(omega,n)

%--------------------------------------------------------------------------
% Object: 
%--------------------------------------------------------------------------
% Inputs:
%--------------------------------------------------------------------------
% Outputs:
%--------------------------------------------------------------------------
% Authors G. Meurant and A. Sommariva
% June 2012
%--------------------------------------------------------------------------

mom=zeros(1,n+1);
mom(1)=2*omega; % FIRST MOMENT.

if(n>=2)
    
    if(omega<=1/4*pi)
        l=10;
    elseif(omega<=1/2*pi)
        l=20;
    elseif(omega<=3/4*pi)
        l=40;
    else
        if omega == pi
            l=2*ceil(10*pi);
        else
            l=2*ceil(10*pi/(pi-omega));
        end
    end
    
    
    temp=(2:2:n+2*l-2); % AUXILIAR VECTORS.
    temp2=temp.^2-1;
    
    dl=1/4 -1./(4*(temp-1)); % DIAGONALS.
    dc=1/2 -1/sin(omega/2)^2 -1./(2*temp2);
    du=1/4 +1./(4*(temp+1));
    
    d=4*cos(omega/2)/sin(omega/2)./temp2'; % COMPUTING KNOWN TERM.
    d(end)=d(end);                         % PUT LAST MOMENT NULL
    
    z=tridisolve(dl(2:end),dc,du(1:end-1),d); % SOLVE SYSTEM.
    mom(3:2:n+1)=z(1:floor(n/2)); % SET ODD MOMENTS.
    
end

mom=mom';

normalized = 0;

if normalized == 0
    M=length(mom);
    kk=2.^(-((1:2:M)-2))'; kk(1)=1;
    v=ones(M,1);
    v(1:2:M)=kk;
    mom=v.*mom;
end









function xw=SymmMw(N,ab)

%--------------------------------------------------------------------------
% Object:
% Computation of the nodes and weights for a symmetric weight function
% this version uses the reduced matrix and eig and computation of weights 
% with the 3-term recurrence.
%--------------------------------------------------------------------------
% Input:
% N : cardinality of the rule
% ab: 3-term recurrence for the orthogonal polynomials same as in OPQ
%     ab(1,2) is the 0th moment.
%--------------------------------------------------------------------------
% Output:
% xw : xw(:,1) nodes, xw(:,2) weights of the quadrature rule
%--------------------------------------------------------------------------
% Reference paper:
% Fast variants of the Golub and Welsch algorithm for symmetric
% weight functions by G. Meurant and A. Sommariva (2012)
%--------------------------------------------------------------------------
% Data:
% Written by G. Meurant and A. Sommariva on June 2012
%--------------------------------------------------------------------------

N0 = size(ab,1);
if N0 < N
    error('SymmMw: input array ab is too short')
end

na = norm(ab(:,1));
if na > 0
    error('SymmMw: the weight function must be symmetric')
end

% computation of the reduced matrix in vectors (a,b)

if mod(N,2) == 0
    even = 1;
    Nc = N / 2;
else
    even = 0;
    Nc = fix(N / 2) +1;
end


absd = ab(:,2);
absq = sqrt(absd);

a = zeros(1,Nc);
b = a;

switch even
    case 1
        % N even
        a(1) = absd(2);
        b(1) = absq(2) * absq(3);
        
        k = [2:Nc-1];
        a(k) = absd(2*k-1) + absd(2*k);
        b(k) = absq(2*k) .* absq(2*k+1);
        a(Nc) = absd(N) + absd(N-1);
        start = 1;
        
        J = diag(a) + diag(b(1:Nc-1),1) + diag(b(1:Nc-1),-1);
        t = sort(eig(J));
        w = weights_3t(t',a,b);
        % w are the squares of the first components
        w = w' / 2;
    case 0
        % N odd
        a(1) = absd(2);
        b(1) = absq(2) * absq(3);
        
        k = [2:Nc-1];
        a(k) = absd(2*k-1) + absd(2*k);
        b(k) = absq(2*k) .* absq(2*k+1);
        a(Nc) = absd(N);
        start = 2;
        
        % the first node must be zero
        J = diag(a) + diag(b(1:Nc-1),1) + diag(b(1:Nc-1),-1);
        t = sort(eig(J));
        t(1) = 0;
        w = weights_3t(t',a,b);
        w = [w(1); w(2:end)' / 2];
    otherwise
        error('this is not possible')
end

xwp = sqrt(t);

xw(:,1) = [-xwp(end:-1:start,1); xwp];
xw(:,2) = ab(1,2) * ([w(end:-1:start); w]);









function tw=quad_trig_conversion(xw,omega)

%--------------------------------------------------------------------------
% Object: 
%--------------------------------------------------------------------------
% Inputs:
%--------------------------------------------------------------------------
% Outputs:
%--------------------------------------------------------------------------
% Authors G. Meurant and A. Sommariva
% June 2012
%--------------------------------------------------------------------------

tw(:,1)=2*asin(sin(omega/2)*xw(:,1));
tw(:,2)=xw(:,2);









function w=weights_3t(t,a,b)

%--------------------------------------------------------------------------
% Object: 
% Squares of the 1st components of eigenvectors from the 3-term
% recurrence relation of the orthogonal polynomials
%--------------------------------------------------------------------------
% Inputs:
% t: nodes
% a,b: coefficients of the 3-term recurrence
%--------------------------------------------------------------------------
% Outputs
% w: squares of the first components of the eigenvectors
%--------------------------------------------------------------------------
% Authors G. Meurant and A. Sommariva
% June 2012
%--------------------------------------------------------------------------

N = length(t);

P = zeros(N,N);
P(1,:) = ones(1,N);
P(2,:) = (t - a(1)) / b(1);

for k = 3:N
    k1 = k - 1;
    k2 = k - 2;
    P(k,:) = ((t - a(k1)) .* P(k1,:) - b(k2) * P(k2,:)) / b(k1);
end

P2 = P .* P;
w = 1 ./ sum(P2);


