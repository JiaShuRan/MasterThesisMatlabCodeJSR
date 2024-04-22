function ab=r_jacobi(N,a,b) 

% r_jacobi is used to compute the a_n-coefficients and the b_n-coefficients
% of the three term recurrence p_{n+1}(x) = (x-a_n)p_n(x) - b_np_{n-1}(x)
% Input: N the degree of the polynomial p_{n+1}(x)
%        a 0 (for Gauss-Legendre)
%        b 0 (for Gauss-Legendre)
% Output:ab is a N*2 matrix, where the first column is [a_1,a_2,...,a_N]^T
%        the second column is [b_1,b_2,...,b_N]^T.

% R_JACOBI Recurrence coefficients for monic Jacobi polynomials.
%
%    ab=R_JACOBI(n,a,b) generates the first n recurrence 
%    coefficients for monic Jacobi polynomials with parameters 
%    a and b. These are orthogonal on [-1,1] relative to the
%    weight function w(t)=(1-t)^a(1+t)^b. The n alpha-coefficients
%    are stored in the first column, the n beta-coefficients in
%    the second column, of the nx2 array ab. The call ab=
%    R_JACOBI(n,a) is the same as ab=R_JACOBI(n,a,a) and
%    ab=R_JACOBI(n) the same as ab=R_JACOBI(n,0,0).
%
%    Supplied by Dirk Laurie, 6-22-1998; edited by Walter
%    Gautschi, 4-4-2002.

nu=(b-a)/(a+b+2); % 0
mu=2^(a+b+1)*gamma(a+1)*gamma(b+1)/gamma(a+b+2); % 2
if N==1
 ab=[nu mu]; return
end

N=N-1;
n=1:N;
nab=2*n+a+b; % (N-1)*1 vector
nuadd=(b^2-a^2)*ones(1,N)./(nab.*(nab+2)); % (N-1)*1 vector 
A=[nu nuadd]; % N*1 vector
n=2:N; % (N-2)*1 vector
nab=nab(n); % (N-2)*1 vector
B1=4*(a+1)*(b+1)/((a+b+2)^2*(a+b+3)); % 1/3
B=4*(n+a).*(n+b).*n.*(n+a+b)./((nab.^2).*(nab+1).*(nab-1)); % (N-2)*1 vector
abadd=[mu; B1; B']; % N*1 vector
ab=[A' abadd];