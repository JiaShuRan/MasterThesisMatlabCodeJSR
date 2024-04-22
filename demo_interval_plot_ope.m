function demo_interval_plot_ope
% Codes based on Alvise Sommariva (University of Padova)
% Date: 26 Sep, 2023
% demo_interval_plot_new is used to plot original function
% and its recovery function via different hyperinterpolants
% on the interval [-1,1]

%--------------------------------------------------------------------------
% Degrees of precision in numerical experiments: can be a vector.
%-------------------------------------------------------------------------
%LV=250;      % Hyperinterpolant total degree.
L=250;

NV=251;     % Points of the rule used for hyperintepolation. 
NR=300;      % Points used by the reference rule (testing L2 errors).


%--------------------------------------------------------------------------
% Noise and choice of lasso, hybrid, hard thresh. parameter.
%--------------------------------------------------------------------------



% No table or stats.
display_stats=1;

% * Function to approximate (see approx. line 312):
funct_example=3; % set 1 or 2 for polynomials of degree L or L/2 (approx.)


% ........................ Main code below ................................

% vectors used for statistics
AEinfMV=[]; AE2MV=[]; beta0MV=[]; 

% Test points
XWR=quad_gausslegendre(NR); XR=XWR(:,1); WR=XWR(:,2);

% define quadrature rule for hyperinterpolation at degree NV.
XW=quad_gausslegendre(NV); X=XW(:,1); W=XW(:,2);


poly_coeffs=[];
lambdaV=[];

% ... define function to approximate ...
g=define_function(funct_example,L);

% ... evaluate function to approximate ...
gX=feval(g,X);


sigma = 0.15;

% add gaussian noise
pert = sigma*randn(size(gX));


% perturbed values
gX_pert=gX+pert;


Candiate = [4,10];

LV_initial = 10;
LV_final = 250;
LV_step = 10;

% Vandermonde matrix at nodes.
V1=vandermonde_gausslegendre(50,X);

% ... determine polynomial hyperinterpolant ...
coeff1=(gX_pert.*W)'*V1; coeff1=coeff1';
noise1=(pert.*W)'*V1; noise1 = noise1';
lambdak=sort(abs(coeff1),'descend');
for k=1:length(coeff1)
    coeff1_new = coeff1.*(abs(coeff1) > lambdak(k));
    Jz(k) = 2*coeff1_new'*noise1 - coeff1_new'*coeff1_new;
    sgn_coeff = sgnfun(coeff1_new)';
    noise1_new = noise1.*(abs(coeff1) > lambdak(k));
    Hz(k) = sum(abs(sgn_coeff))*lambdak(k)^2 - 2*lambdak(k)*sgn_coeff'*noise1_new;
end

%Hz_great_0 = Hz.*(Hz >=0);

for p=1:length(Hz)
    if Hz(p) >= 0
        Hz_great_0(p) = Hz(p);
    else
        Hz_great_0(p) = NaN;
    end
end


Val_f = sum((gX.^2).*W);
Val_f = ones(1,length(coeff1))*Val_f;

AErr_lasso = sqrt(Jz+Val_f+Hz);
AErr_hard = sqrt(Jz+Val_f);






% Vandermonde matrix at nodes.
V=vandermonde_gausslegendre(LV,X);

% ... determine polynomial hyperinterpolant ...
coeff0=(gX_pert.*W)'*V; coeff0=coeff0';
noise0=(pert.*W)'*V; noise0 = noise0';

lambdas=sort(abs(coeff0),'descend');

lambdaL=lambdas(4);

degs=0:LV;

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
    VR=vandermonde_gausslegendre(LV,XR);
    pXR(:,ktest)=VR*coeff;

    % errors
    AEinfV(jtest,ktest)=norm(gXR-pXR(:,ktest),inf); % absolute error (inf norm)
    AE2V(jtest,ktest)=sqrt(WR'*((gXR-pXR(:,ktest)).^2)); % absolute error (2 norm)
    beta0V(jtest,ktest)=sum(abs(coeff) > 0);

end




%% Plot

%% Step 3

figure(1)
subplot(1,2,1)

plot(1:length(noise1),Hz,'-.o','linewidth',1,'MarkerSize',10,'color','k'), box on, set(gca,'fontsize',16),
set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),
hold on

plot(1:length(noise1),Hz_great_0,'o','linewidth',1,'MarkerSize',10,'color','red',"MarkerFaceColor","red"), box on, set(gca,'fontsize',16),
set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),
hold on

xlabel({'\textbf{Sparsity} $(c)$'},'interpreter','latex','fontsize',35);ylabel({'\textbf{Values}'},'interpreter','latex','fontsize',35);
legend({'\textbf{Negative}','\textbf{Nonnegative}'},'interpreter','latex','fontsize',30);
%title({'$J(z)$ \textbf{for hard thresholding and Lasso hyperinterpolations}'},'interpreter','latex','fontsize',35);

subplot(1,2,2)
semilogy(1:length(noise1),sqrt(Jz+Val_f),'bd','linewidth',1,'markersize',18), box on, set(gca,'fontsize',16),
set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),
hold on

semilogy(1:length(noise1),sqrt(Hz+Jz+Val_f),'rpentagram','linewidth',1,'markersize',15,"MarkerFaceColor","red"), box on, set(gca,'fontsize',16),
set(gca, 'XMinorGrid', 'on'), set(gca, 'YMinorGrid', 'on'),
hold on

xlabel({'\textbf{Sparsity} (d)'},'interpreter','latex','fontsize',35);ylabel({'$L_2$ \textbf{Errors}'},'interpreter','latex','fontsize',35);
%legend({'\textbf{Hard thresholding hyper.}','\textbf{Lasso hyper.}'},'interpreter','latex','fontsize',30);
%title({'$L_2$ \textbf{Error estimate for hard thresholding and Lasso hyperinterpolations}'},'interpreter','latex','fontsize',35);




end


%% ATTACHED ROUTINES.

function sgn_coeff = sgnfun(coeff0)

for k=1:length(coeff0)
    if coeff0(k) > 0
        sgn_coeff(k) = 1;
    elseif coeff0(k) < 0
        sgn_coeff(k) = -1;
    else
        sgn_coeff(k) = 0;
    end
end

end


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

end



function plot_error(f,coeff,n)

XR=-1:0.01:1; XR=XR';

fXR=feval(f,XR);
VR=vandermonde_gausslegendre(n,XR);
pXR=VR*coeff;
aerr=abs(fXR-pXR);
rerr=aerr./abs(fXR);

semilogy(XR,rerr);

end




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

end



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
                %g=@(x) 4*(x.^5).*sin(10*x); fstring = '4*(x.^5).*sin(10*x)';
            case 3
                g=@(x) cos(-x.*2); fstring='cos(-x.*2)';
            case 4
                g=@(x) sin(x)+x.^2; fstring='sin(x)+x.^2';
        end

    case 4
        g=@(x) 1./(1+25*x.^2);
end



end