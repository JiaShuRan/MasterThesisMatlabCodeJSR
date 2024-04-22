

function [AEinfMV,AE2MV,beta0MV,lambdaV,XW,XWR,JzMV,HzMV,Wg_norm2]=...
    demo_interval(lambda_index,a,sigma,XW,XWR)

%--------------------------------------------------------------------------
% OBJECT
%--------------------------------------------------------------------------
% Numerical experiment in "Hybrid hyperinterpolation over general regions".
% Region: unit-interval [-1,1].
%--------------------------------------------------------------------------
% Usage:
% >> demo_interval
%--------------------------------------------------------------------------
% Note:
% The routine uses 'binornd' that requires Statistics and Machine Learning 
% Toolbox.
%--------------------------------------------------------------------------
% Dates:
% Written on January 1, 2023: A. Sommariva.
% Modified on April 23, 2023: A. Sommariva.
%--------------------------------------------------------------------------
% COPYRIGHT
%--------------------------------------------------------------------------
% Copyright (C) 2023- 
%
% Authors:
% Alvise Sommariva
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Degrees of precision in numerical experiments: can be a vector.
%-------------------------------------------------------------------------
LV=250;      % Hyperinterpolant tot degree.
NV=LV+1;     % Points of the rule used for hyperintepolation.
NR=300;      % Points used by the reference rule (testing L2 errors).

%--------------------------------------------------------------------------
% Noise and choice of lasso, hybrid, hard thresh. parameter.
%--------------------------------------------------------------------------

% defining hybrid parameter (knowning hyp. coeffs, it is the k-th in 
% descending absolute value!)
if nargin < 1, lambda_index=1; end

noise=1;            % 0: no noise, 1: noise (see parameters below)

% In the reference paper of Lasso (Table 1): a=0, sigma=0.2.
if noise
    if nargin <2,a=0; end      % defining impulse noise (in experiment 2)
    if nargin <3,sigma=0.4;end % defining gaussian noise (in experiment 2)
else
    a=0; sigma=0; % no noise.
end

% * Function to approximate (see approx. line 312):
funct_example=3; % set 1 or 2 for polynomials of degree L or L/2 (approx.)

% No table or stats.
display_stats=0;

%--------------------------------------------------------------------------
% Special settings.
%--------------------------------------------------------------------------

% Plot domain and nodes: do_plot=1 (yes), do_plot=0 (no).
do_plot=0;

% Number of tests for reconstructing functions of a type on this domain.
ntests=100;







% ........................ Main code below ................................


AEinfMV=[]; AE2MV=[]; beta0MV=[]; % vectors used for statistics
JzMV=[]; HzMV=[];

for k=1:length(NV)
    N=NV(k); % Quadrature points.
    L=LV(k); % Hyperinterpolant degree.

    % Test points
    if nargin < 5, XWR=quad_gausslegendre(NR); end
    if isempty(XWR), XWR=quad_gausslegendre(NR); end
    XR=XWR(:,1); WR=XWR(:,2);

    % define quadrature rule for hyperinterpolation at degree N.
    if nargin < 4, XW=quad_gausslegendre(N); end
    if isempty(XW), XW=quad_gausslegendre(N); end
    X=XW(:,1); W=XW(:,2);

    % Vandermonde matrix at nodes.
    % compute hyperinterpolant coefficients
    V=vandermonde_gausslegendre(L,X);

    % .. testing AE_L2err hyperinterpolation error for each "f" at "deg" ..

    poly_coeffs=[];
    lambdaV=[];

    for j=1:ntests

        % ... define function to approximate ...
        g=define_function(funct_example,L);

        % ... evaluate function to approximate ...
        gX=feval(g,X);

        % ... Add noise (if present) ...

        % add impulse noise
        pert_impulse=0;
        if a > 0
            pert_impulse=a*(1-2*rand(length(gX),1))*binornd(1,0.5);
            while norm(pert_impulse) == 0
                pert_impulse=a*(1-2*rand(length(gX),1))*binornd(1,0.5);
            end
        end

        % add gaussian noise
        pert_gauss=0;
        if sigma > 0
            var=sigma^2;
            pert_gauss=sqrt(var)*randn(size(gX));
            while norm(pert_gauss) == 0
                pert_gauss=sqrt(var)*randn(size(gX));
            end
        end

        % add gaussian + impulse noise
        pert=pert_impulse+pert_gauss;

        % perturbed values
        gX_pert=gX+pert;

        % ... determine polynomial hyperinterpolant ...
        coeff0=(gX_pert.*W)'*V; coeff0=coeff0';

        degs=0:L;

        % test hyperinterpolant with or withour filters.

       lambdas=sort(abs(coeff0),'descend');
       lambdaL=lambdas(lambda_index);

        for ktest=1:6
            switch ktest
                case 1
                    hypermode='tikhonov';
                    parms.lambda=lambdaL;
                    parms.mu=[];
                    parms.b=ones(size(coeff0));
                    coeff=hyperfilter(hypermode,coeff0,degs,parms);

                case 2
                    hypermode='filtered';
                    parms.lambda=[];
                    parms.mu=[];
                    parms.b=[];
                    coeff=hyperfilter(hypermode,coeff0,degs,parms);
                case 3
                    hypermode='lasso';
                    parms.lambda=lambdaL;
                    parms.mu=ones(size(coeff));
                    parms.b=[];
                    coeff=hyperfilter(hypermode,coeff0,degs,parms);
                case 4
                    hypermode='hybrid';
                    parms.lambda=lambdaL;
                    parms.mu=ones(size(coeff0));
                    parms.b=ones(size(coeff0));
                    parms.w=W;
                    parms.pert=pert;
                    parms.hybrid=0; % establishes it is a pre-choosen parameter.
                    coeff=hyperfilter(hypermode,coeff0,degs,parms);
              case 5
                    hypermode='hard';
                    parms.lambda=lambdaL;
                    parms.mu=[];
                    parms.b=[];
                    coeff=hyperfilter(hypermode,coeff0,degs,parms);
                case 6
                    hypermode='hyperinterpolation';
                    parms.lambda=[];
                    parms.mu=[];
                    parms.b=[];
                    coeff=hyperfilter(hypermode,coeff0,degs,parms);
            end

            gXR=feval(g,XR);
            VR=vandermonde_gausslegendre(L,XR);
            pXR=VR*coeff;

            % errors
            AEinfV(ktest,j)=norm(gXR-pXR,inf); % absolute error (inf norm)
            AE2V(ktest,j)=sqrt(WR'*((gXR-pXR).^2)); % absolute error (2 norm)
            beta0V(ktest,j)=sum(abs(coeff) > 0);

            % evaluate J(coeff) and H(coeff), that are error relevant 
            % parameters, as observed in Thm 5.1.

            JzV(ktest,j)=evaluate_J(coeff,coeff0);
            HzV(ktest,j)=evaluate_H(coeff,W,pert,V);

        end

        lambdaV=[lambdaV lambdas];

    end

    % averages of the errors (vectors 5 x 1)
    AEinfM=mean(AEinfV,2);
    AE2M=mean(AE2V,2);
    beta0M=mean(beta0V,2);
    JzM=mean(JzV,2);
    HzM=mean(HzV,2);

    if display_stats
    fprintf('\n       ........ table at degree: %2.0f ........ \n \n ',N);
    HypType=categorical({'tikhonov'; 'filtered'; 'lasso'; 'hybrid'; ...
        'hard'; 'hyperint.'});
    T = table(HypType,AEinfM,AE2M,beta0M,JzM,HzM); disp(T)
    end
    
    AEinfMV=[AEinfMV AEinfM]; AE2MV=[AE2MV AE2M]; beta0MV=[beta0MV beta0M];
    JzMV=[JzMV JzM]; HzMV=[HzMV HzM];

end

Wg_norm2=(norm(sqrt(W).*gX,2))^2;




%  fprintf('\n \t Hyp. coeff.: %5.0f card rule 1: %5.0f  card rule 2: %5.0f ',...
%      length(coeff0),length(X),length(XR));





%--------------------------------------------------------------------------
% ATTACHED ROUTINES.
%--------------------------------------------------------------------------

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





function plot_error(f,coeff,n)

XR=-1:0.01:1; XR=XR';

fXR=feval(f,XR);
VR=vandermonde_gausslegendre(n,XR);
pXR=VR*coeff;
aerr=abs(fXR-pXR);
rerr=aerr./abs(fXR);

semilogy(XR,rerr);





function V=vandermonde_gausslegendre(deg,x)

%--------------------------------------------------------------------------
% Object:
% Orthonormal Gauss-Legendre Vandermonde matrix.
% This basis is numerically almost orthonormal even at high degrees.
%--------------------------------------------------------------------------
% Input:
% deg: degree of the polynomial basis;
% x  : M points in which the basis is evaluated.
%--------------------------------------------------------------------------
% Output:
% V  : Vandermonde matrix of degree "deg" at "x". It is a "(deg+1) x M"
%      matrix.
%--------------------------------------------------------------------------

V=ones(size(x));

if deg >= 1
    V=[V x];
end

for n=1:deg-1
    VO=V(:,end); VOO=V(:,end-1);
    VL=((2*n+1)*x.*VO-n*VOO)/(n+1);
    V=[V VL];
end

N=0:deg;
s=1./sqrt(2./(2*N+1));
S=repmat(s,size(V,1),1);
V=V.*S;





function g=define_function(funct_example,L)

%--------------------------------------------------------------------------
% Object:
% Define function to be tested.
%--------------------------------------------------------------------------
% Input:
% funct_example: parameter defining the function to be tested
%                1: random polynomial of form (c0+c1*x)^L.
%                2: random polynomial of form (c0+c1*x)^n with
%                   n=max(floor(L/2)-1,0).
%                3: exp(-x^2)
%--------------------------------------------------------------------------
% Output:
% g: function to test.
%--------------------------------------------------------------------------


switch funct_example

    case 1 % test exactness hyperinterpolation
        nexp=L;
        c0=rand(1); c1=rand(1);
        g=@(x) (c0+c1*x).^nexp;

    case 2 % test exactness filt. hyperinterpolation
        nexp=max(floor(L/2)-1,0);
        c0=rand(1); c1=rand(1); c2=rand(1);
        g=@(x) (c0+c1*x).^nexp;

    case 3 % function of that type

        %  Test function "exp(-x.^2)" has been choosen on the Lasso
        %  paper.

        funct_example_sub=2;

        switch funct_example_sub
            case 1
                g=@(x) exp(-2*x); fstring='exp(-2*x)';
            case 2
                g=@(x) exp(-x.^2); fstring='exp(-x.^2)';
            case 3
                g=@(x) cos(-x.*2); fstring='cos(-x.*2)';
            case 4
                g=@(x) sin(x)+x.^2; fstring='sin(x)+x.^2';
        end

end





function Jz=evaluate_J(z,alpha)

%--------------------------------------------------------------------------
% Object:
% Evaluate function J(z)=sum( (z(l))^2-2*z(l)*alpha(l) )
%--------------------------------------------------------------------------
% Input:
% z    : vector of dimension d x 1
% alpha: vector of dimension d x 1
%--------------------------------------------------------------------------
% Output:
% Jz: value of J(z)=sum( (z(l))^2-2*z(l)*alpha(l) )
%--------------------------------------------------------------------------
% Reference:
% Quantity relevant in Thm. 5.1 of the paper
% "Hybrid hyperinterpolation over general regions"
% Congpei An · Alvise Sommariva · Jia-Shu Ran
%--------------------------------------------------------------------------

% Jz=sum(z.^2-2*z.*alpha);
Jz=z'*z -2*z'*alpha;




function Hz=evaluate_H(z,w,err,V)

%--------------------------------------------------------------------------
% Object:
% Evaluate function H(z)=2*sum_l( z(l) * sum_j( w(j)*err(j)*V(l,j) ) )
%--------------------------------------------------------------------------
% Input:
% z    : vector of dimension d x 1
% w    : vector of dimension N x 1
% err  : vector of dimension N x 1
% V    : matrix of dimension d x N
%--------------------------------------------------------------------------
% Output:
% Hz: value of H(z)=2*sum_l( z(l) * sum_j( w(j)*err(j)*V(l,j) ) )
%--------------------------------------------------------------------------
% Reference:
% Quantity relevant in Thm. 5.1 of the paper
% "Hybrid hyperinterpolation over general regions"
% Congpei An · Alvise Sommariva · Jia-Shu Ran
%--------------------------------------------------------------------------

inner_term=V'*(w.*err);
outer_term=z'*inner_term;
Hz=2*outer_term;





