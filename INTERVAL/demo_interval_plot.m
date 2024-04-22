function demo_interval_plot
% Codes based on Alvise Sommariva (University of Padova)
% Date: 14 August, 2023
% demo_interval_plot is used to plot original function
% and its recovery function via different hyperinterpolants
% on the interval [-1,1]

%--------------------------------------------------------------------------
% Degrees of precision in numerical experiments: can be a vector.
%-------------------------------------------------------------------------
LV=250;      % Hyperinterpolant total degree.
NV=251;     % Points of the rule used for hyperintepolation. 
NR=300;     % Points used by the reference rule (testing L2 errors).


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

% Vandermonde matrix at nodes.
V=vandermonde_gausslegendre(LV,X);

poly_coeffs=[];
lambdaV=[];

% ... define function to approximate ...
g=define_function(funct_example,LV);

% ... evaluate function to approximate ...
gX=feval(g,X);


sigma = 0.1;

% add gaussian noise
pert = sigma*randn(size(gX));


% perturbed values
gX_pert=gX+pert;

% ... determine polynomial hyperinterpolant ...
coeff0=(gX_pert.*W)'*V; coeff0=coeff0';


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
            AEinfV(ktest)=norm(gXR-pXR(:,ktest),inf); % absolute error (inf norm)
            AE2V(ktest)=sqrt(WR'*((gXR-pXR(:,ktest)).^2)); % absolute error (2 norm)
            beta0V(ktest)=sum(abs(coeff) > 0);
            coeff_all(:,ktest)=coeff;
         

end


%% Plot
figure(1)
Color = parula(7);
fontsize_baseline = 35;
fontsize_baselinet = 25;
fontsize_baselinea = 15;

%Primal and noisy function
subplot(2,4,1)
%axes('position',[0.075 0.55 0.2 0.4]), 
plot(X,gX,'linewidth',2,'color','k'), box on,...
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('目标函数','fontsize', fontsize_baselinet),...
     grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-1,1,0,1.5]),
subplot(2,4,5)
%axes('position',[0.075 0.05 0.2 0.4]), 
plot(X,gX_pert,'linewidth',2,'color','k'),box on,...
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('含噪声的目标函数','fontsize', fontsize_baselinet),...
     grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-1,1,0,1.5]),

% hyper. and its error
subplot(2,4,2)
%axes('position',[0.3 0.55 0.2 0.4]),
plot(XR,pXR(:,6),'linewidth',2,'color','k'), box on,...
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('超插值', 'fontsize', fontsize_baselinet),  grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-1,1,0,1.5]),
%axes('position',[0.3 0.05 0.2 0.4]), 
subplot(2,4,6)
plot(XR,abs(gXR-pXR(:,6)),'linewidth',2,'color','k'),set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('Absolute error', 'interpreter','latex','fontsize', fontsize_baseline),...
    title('绝对值误差','fontsize', fontsize_baselinet),box on, grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-1,1,0,1]),

% hard hyper. and its error
%axes('position',[0.525 0.55 0.2 0.4]),
subplot(2,4,3)
plot(XR,pXR(:,5),'linewidth',2,'color','k'), box on,...
    set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('硬阈值超插值','fontsize', fontsize_baselinet),  grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-1,1,0,1.5]),
%axes('position',[0.525 0.05 0.2 0.4]), 
subplot(2,4,7)
plot(XR,abs(gXR-pXR(:,5)),'linewidth',2,'color','k'),set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('Absolute error', 'interpreter','latex','fontsize', fontsize_baseline),...
    title('绝对值误差','fontsize', fontsize_baselinet),box on, grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-1,1,0,1]),

%hybrid hyper. and its error
%axes('position',[0.75 0.55 0.2 0.40]),
subplot(2,4,4)
plot(XR,pXR(:,4),'linewidth',2,'color','k'),box on,set(gca, 'fontsize', fontsize_baselinea), ...
     xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
     title('混合超插值','fontsize', fontsize_baselinet),...
     grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-1,1,0,1.5]),
%axes('position',[0.75 0.05 0.2 0.40]), 
subplot(2,4,8)
plot(XR,abs(gXR-pXR(:,4)),'linewidth',2,'color','k'),set(gca, 'fontsize', fontsize_baselinea), xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('Absolute error','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('绝对值误差','fontsize', fontsize_baselinet),box on,grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-1,1,0,1]),


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