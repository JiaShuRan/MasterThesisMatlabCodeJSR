function demo_disk_plot
% Author: Jia-Shu Ran
% Date: 15 August, 2023
% demo_interval_plot is used to plot original function
% and its recovery function via different hyperinterpolants
% on the unit-disk

LV=16;        % Hyperinterpolant tot degree.
NV=80;        % Degree of precision of the rule.
NR=100;       % Reference degree of cubature rule.

%--------------------------------------------------------------------------
% Noise and choice of lasso, hybrid, hard thresh. parameter.
%--------------------------------------------------------------------------

% * Function to approximate:
% 1. degree L poly., 2. degree floor(L/2)-1 poly. 3. test functions 
% (see line 300 approx).
funct_example=3;  

AEinfMV=[]; AE2MV=[]; beta0MV=[]; % vectors used for statistics

% Test points
XYWR=cub_disk_productrule(NR); XR=XYWR(:,1); YR=XYWR(:,2); WR=XYWR(:,3);

% define quadrature rule for hyperinterpolation at degree N.
XYW=cub_disk_productrule(NV); X=XYW(:,1); Y=XYW(:,2); W=XYW(:,3);

% Vandermonde matrix at nodes.
% compute hyperinterpolant coefficients
[V,degs]=vandermonde_logan_shepp(LV,[X Y]);

% ... define function to approximate ...
g=define_function(funct_example,LV);

% ... evaluate function to approximate ...
gXY=feval(g,X,Y);

% add impulse noise 
a =0; %pert=a*(1-2*rand(length(gXY),1))*binornd(1,0.5);

sigma = 0.1;

% add gaussian noise
pert = sigma*randn(size(gXY)) + a*(1-2*rand(length(gXY),1))*binornd(1,0.5);

% perturbed values
gXY_pert=gXY+pert;

% ... determine polynomial hyperinterpolant ...
coeff0=(gXY_pert.*W)'*V; coeff0=coeff0';

degs=0:LV;

% test hyperinterpolant with or withour filters.

lambdas=sort(abs(coeff0),'descend');
lambdaL=lambdas(12);

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

    gXYR=feval(g,XR,YR);
    [VR,degs]=vandermonde_logan_shepp(LV,[XR YR]);


    pXYR(:,ktest)=VR*coeff;

    % errors
    AEinfV(ktest)=norm(gXYR-pXYR(:,ktest),inf); % absolute error (inf norm)
    AE2V(ktest)=sqrt(WR'*((gXYR-pXYR(:,ktest)).^2)); % absolute error (2 norm)
    beta0V(ktest)=sum(abs(coeff) > 0);

    [pXYR_re, abs_err] = evaluate_disk(g,coeff,LV);

    pXYR_re_plot(ktest).matrix = pXYR_re;
    abs_err_plot(ktest).matrix = abs_err;

end

%[noisy_function,real_function ]= noise_fun(g,a);
[noisy_function,real_function ]= noise_fun(g,sigma,a);
t = -1:0.025:1;
[XX, YY] = meshgrid(t);
%% Plotting
fontsize_baselinea = 10;
fontsize_baseline = 15;
fontsize_baselinet = 25;

% Primal and noisy function
%subplot(2,4,1)
axes('position',[0.075 0.55 0.2 0.4]),
mesh(XX,YY,real_function,'edgecolor','flat'), set(gca, 'fontsize', fontsize_baselinea), box on,...
    xlabel('$x_1$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$x_2$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('目标函数','fontsize', fontsize_baselinet),...
    grid on, axis([-1,1,-1,1,-0.5,1.5]),
axes('position',[0.075 0.05 0.2 0.4]),
%subplot(2,4,5)
mesh(XX,YY,noisy_function,'edgecolor','flat'), set(gca, 'fontsize', fontsize_baselinea),box on,...
    xlabel('$x_1$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$x_2$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('含噪声的目标函数','fontsize', fontsize_baselinet),...
    grid on,...
    axis([-1,1,-1,1,-0.5,1.5]),

% hyper. and its error
axes('position',[0.3 0.55 0.2 0.4]),
mesh(XX,YY,pXYR_re_plot(6).matrix,'edgecolor','flat'),set(gca, 'fontsize', fontsize_baselinea),box on,...
    xlabel('$x_1$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$x_2$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('超插值','fontsize', fontsize_baselinet),  grid on,axis([-1,1,-1,1,-.5,1.5]),

axes('position',[0.3 0.05 0.2 0.4]),
mesh(XX,YY,abs_err_plot(6).matrix,'edgecolor','flat'),set(gca, 'fontsize', fontsize_baselinea), xlabel('$x_1$','interpreter','latex', 'fontsize', fontsize_baseline), ylabel('$x_2$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('Absolute error', 'interpreter','latex','fontsize', fontsize_baseline),...
    title('绝对值误差','fontsize', fontsize_baselinet),box on,  grid on,...
    set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),axis([-1,1,0,1.5]),
axis([-1,1,-1,1,0,0.5]),

% Hard hyepr. and its error
axes('position',[0.525 0.55 0.2 0.4]),
mesh(XX,YY,pXYR_re_plot(5).matrix,'edgecolor','flat'), set(gca, 'fontsize', fontsize_baselinea),box on,...
    xlabel('$x_1$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$x_2$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('硬阈值超插值','fontsize', fontsize_baselinet),...
    grid on,  set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),
axis([-1,1,-1,1,-.5,1.5]),

axes('position',[0.525 0.05 0.2 0.4]),
mesh(XX,YY,abs_err_plot(5).matrix,'edgecolor','flat'),set(gca, 'fontsize', fontsize_baselinea),xlabel('$x_1$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$x_2$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('Absolute error','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('绝对值误差', 'fontsize', fontsize_baselinet),box on,  grid on,...
    set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),
axis([-1,1,-1,1,0,0.5]),

% Hybrid hyper. and its error
axes('position',[0.75 0.55 0.2 0.4]),
mesh(XX,YY,pXYR_re_plot(4).matrix,'edgecolor','flat'),set(gca, 'fontsize', fontsize_baselinea),box on,...
    xlabel('$x_1$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$x_2$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('$f(x)$','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('混合超插值', 'fontsize', fontsize_baselinet),...
    grid on, set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),
axis([-1,1,-1,1,-.5,1.5]),

axes('position',[0.75 0.05 0.2 0.4]),
mesh(XX,YY,abs_err_plot(4).matrix,'edgecolor','flat'),set(gca, 'fontsize', fontsize_baselinea),xlabel('$x_1$','interpreter','latex', 'fontsize', fontsize_baseline),ylabel('$x_2$','interpreter','latex', 'fontsize', fontsize_baseline),...%ylabel('Absolute error','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('绝对值误差','fontsize', fontsize_baselinet),box on, grid on,...
    set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),  axis([-1,1,-1,1,0,0.5])

end

%% Functions used in this programm

function XYW=cub_disk_productrule(N)

% PRODUCT TYPE RULE ON THE UNIT DISK: ade at least "N". 

n=ceil(N/2);
gaussian_type=1;

switch gaussian_type
    case 1
        % Gaussian rule in [0,1]: ade=2*n -> n+1 points.
        ab=r_jacobi(n+1,0,0); xw=gauss(n+1,ab);
        r0=xw(:,1); wr0=xw(:,2); r=(r0+1)/2; wr=wr0/2;
        %     case 2
        %         xw=lobatto_jacobi(n,0,0); % ADE: 2*n+1 using n+2 points.
        %         r0=xw(:,1); wr0=xw(:,2); r=(r0+1)/2; wr=wr0/2;
        %     case 3
        %         xw=radau_jacobi(n,1,0,0); % ADE: 2*n using n+1 points.
        %         r0=xw(:,1); wr0=xw(:,2); r=(r0+1)/2; wr=wr0/2;
end

% Trapezoidal rule.
m=0:2*n; th=2*pi*m/(2*n+1); th=th';
wth=2*pi*ones(size(th))/(2*n+1);

% Nodes in polar coordinates.
[R,TH]=meshgrid(r,th);
Rv=R(:); THv=TH(:);
X=Rv.*cos(THv); Y=Rv.*sin(THv);

[W1,W2]=meshgrid(wr.*r,wth);
W=W1.*W2; W=W(:);

XYW=[X Y W];

end





function plot_error_disk(f,coeff,n)

t=-1:0.01:1;
[XR,YR]=meshgrid(t);

fXYR=feval(f,XR,YR);
VR=vandermonde_logan_shepp(n,[XR(:) YR(:)]);
pXYR0=VR*coeff;

pXYR=reshape(pXYR0,size(XR,1),size(XR,2));

val=(ones(size(XR))-XR.^2-YR.^2) >= 0;
err=abs(fXYR-pXYR).*val;

plot3(XR,YR,err);
hold on;
AZ=20; EL=45;
view(AZ,EL);
hold off;

end





function g=define_function(funct_example,L)

% function to test

switch funct_example

    case 1 % test exactness hyperinterpolation
        nexp=L;
        c0=rand(1); c1=rand(1); c2=rand(1);
        g=@(x,y) (c0+c1*x+c2*y).^nexp;

    case 2 % test exactness filt. hyperinterpolation
        nexp=max(floor(L/2)-1,0);
        c0=rand(1); c1=rand(1); c2=rand(1);
        g=@(x,y) (c0+c1*x+c2*y).^nexp;

    case 3 % function of that type

        %  Test function.

        funct_example_sub=1;

        switch funct_example_sub
            case 1
                g=@(x,y) (1-x.^2-y.^2).*exp(x.*cos(y));
                fstring='(1-x.^2-y.^2).*exp(x.*cos(y))';
            case 2
                g=@(x,y) exp((x.^6).*cos(y));
                fstring='exp((x.^6).*cos(y))';
        end
end

end

function [pXYR_re, abs_err] = evaluate_disk(f,coeff,n)

t = -1:0.025:1;
[XR, YR] = meshgrid(t);

fXYR = feval(f, XR, YR);
VR = vandermonde_logan_shepp(n, [XR(:) YR(:)]);

pXYR0 = VR*coeff;

pXYR = reshape(pXYR0,size(XR,1),size(XR,2));

val = (ones(size(XR))- XR.^2 - YR.^2) >=0;

val = val./val;

abs_err = abs(fXYR-pXYR).*val;
pXYR_re = pXYR.*val;
end

%function [noisy_function,real_function ]= noise_fun(f,a)
function [noisy_function,real_function ]= noise_fun(f,sigma,a)
t = -1:0.025:1;
[XR, YR] = meshgrid(t);

fXYR = feval(f, XR, YR);

val = (ones(size(XR))- XR.^2 - YR.^2) >=0;

val = val./val;

%noisy_function = (fXYR + a*(1-2*rand(length(fXYR),1))*binornd(1,0.5) ).*val;
noisy_function = (fXYR + sigma*randn(size(fXYR)) + a*(1-2*rand(length(fXYR),1))*binornd(1,0.5)).*val;
real_function = fXYR.*val;
end

