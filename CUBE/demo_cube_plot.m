function demo_cube_plot
% Codes based on Alvise Sommariva (University of Padova)
% Date: 18 August, 2023
% demo_cube_plot is used to plot original function
% and its recovery function via different hyperinterpolants
% in the cube

LV=10;      % Hyperinterpolant tot degree.
%NV=2*LV;    % Degree of the rule.
NV=40;      % Quadrature exactness
NR=50;      % Reference rule for computing L2 errors.

% * Function to approximate:
% 1. degree L poly., 2. degree floor(L/2)-1 poly. 3. test functions
% (see line 274 approx).
funct_example=3;

% The degree d of the polynomial space is d = (LV+3)(LV+2)(LV+1)/6.
% The quadrature points N = (NV+2)^3/4.



% ........ Numerical approximation, varying the degree in "nV" ............

AEinfMV=[]; AE2MV=[]; beta0MV=[]; % vectors used for statistics

% Define quadrature rule for hyperinterpolation at degree N.
XYZW=cub_cube(NV); X=XYZW(:,1); Y=XYZW(:,2); Z=XYZW(:,3); W=XYZW(:,4);

% Test points
XYZWR=cub_cube(NR); XR=XYZWR(:,1); YR=XYZWR(:,2); ZR=XYZWR(:,3); WR=XYZWR(:,4);

% Vandermonde matrix at nodes and polynomial degrees.
[V,dbox,duples]=dCHEBVAND0(LV,[X Y Z]);
degs=sum(duples,2);

% define function (see attached file at the bottom)
g=choose_function(funct_example,LV);

% ... evaluate function to approximate ...
gXYZ=feval(g,X,Y,Z);


sigma = 0.2; var=sigma^2; pert=sqrt(var)*randn(size(gXYZ));

% add Gaussian noise

gXYZ_pert=gXYZ+pert;

% ... determine polynomial hyperinterpolant ...
coeff0=(gXYZ_pert.*W)'*V; coeff0=coeff0';

% test hyperinterpolant with or withour filters.

lambdas=sort(abs(coeff0),'descend');
lambdaL=lambdas(8);

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
                    coeff=coeff0;
            end

            gXYZR=feval(g,XR,YR,ZR);
            [VR]=dCHEBVAND0(LV,[XR YR ZR]);
            pXYZR(:,ktest)=VR*coeff;

            % errors
            AEinfV(ktest)=norm(gXYZR-pXYZR(:,ktest),inf); % absolute error (inf norm)
            AE2V(ktest)=sqrt(WR'*((gXYZR-pXYZR(:,ktest)).^2)); % absolute error (2 norm)
            beta0V(ktest)=sum(abs(coeff) > 0);

            [pXYZR_3D, abs_err_3D] = evaluate_cube(g,coeff,LV);

            pXYZR_3D_plot(ktest).matrix = pXYZR_3D;
            abs_err_3D_plot(ktest).matrix = abs_err_3D;
end

[XR_3, YR_3, ZR_3] = meshgrid(-1:0.1:1);


gXYZR_3D = feval(g, XR_3, YR_3, ZR_3);
gXYZR_pert_3D = gXYZR_3D + sqrt(var)*randn(size(gXYZR_3D));

%% Plot

fontsize_baselinea = 10;
fontsize_baseline = 10;
fontsize_baselinet = 25;
marksize = 80;

colormap(jet)


xslice = [-.25,.5,1]; 
yslice = [0,1]; 
zslice = [-1,0];

%Primal and noisy function
axes('position',[0.03,0.55,0.2,0.4]) 
slice(XR_3,YR_3,ZR_3,gXYZR_3D,xslice,yslice,zslice),set(gca, 'fontsize', fontsize_baselinea),...
title('目标函数', 'fontsize', fontsize_baselinet),...
grid on, set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),  view(-36,15), colorbar('eastoutside')

axes('position',[0.03,0.05,0.2,0.4])
slice(XR_3,YR_3,ZR_3,gXYZR_pert_3D,xslice,yslice,zslice),set(gca, 'fontsize', fontsize_baselinea),...
title('含噪声的目标函数','fontsize', fontsize_baselinet),...
grid on, set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),  view(-36,15), colorbar('eastoutside'),caxis([0,0.7])

% Hyper and its error
axes('position',[0.27,0.55,0.2,0.4])
slice(XR_3,YR_3,ZR_3,pXYZR_3D_plot(6).matrix,xslice,yslice,zslice),set(gca, 'fontsize', fontsize_baselinea),...
title('超插值','fontsize', fontsize_baselinet),...
grid on, set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'), view(-36,15), colorbar('eastoutside'), caxis([0,0.7])

axes('position',[0.27,0.05,0.2,0.4])
slice(XR_3,YR_3,ZR_3,abs_err_3D_plot(6).matrix,xslice,yslice,zslice),set(gca, 'fontsize', fontsize_baselinea),...
title('绝对值误差','fontsize', fontsize_baselinet),...
grid on, set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),  view(-36,15), colorbar('eastoutside'), caxis([0,0.25])

% Hard hyper. and its error
axes('position',[0.51,0.55,0.2,0.4])
slice(XR_3,YR_3,ZR_3,pXYZR_3D_plot(5).matrix,xslice,yslice,zslice),set(gca, 'fontsize', fontsize_baselinea),...
title('硬阈值超插值','fontsize', fontsize_baselinet),...
grid on, set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'), view(-36,15),colorbar('eastoutside'),caxis([0,.7])
 
axes('position',[0.51,0.05,0.2,0.4])
slice(XR_3,YR_3,ZR_3,abs_err_3D_plot(5).matrix,xslice,yslice,zslice),set(gca, 'fontsize', fontsize_baselinea),... 
title('绝对值误差','fontsize', fontsize_baselinet),...
grid on, set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),  view(-36,15),colorbar('eastoutside'), caxis([0,0.25])

% Hybrid hyper. and its error
axes('position',[0.75,0.55,0.2,0.4])
slice(XR_3,YR_3,ZR_3,pXYZR_3D_plot(4).matrix,xslice,yslice,zslice),set(gca, 'fontsize', fontsize_baselinea),...
title('混合超插值','fontsize', fontsize_baselinet),...
grid on,set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),  view(-36,15), colorbar('eastoutside'),caxis([0,.7])

axes('position',[0.75,0.05,0.2,0.4])
slice(XR_3,YR_3,ZR_3,abs_err_3D_plot(4).matrix,xslice,yslice,zslice),set(gca, 'fontsize', fontsize_baselinea),... 
title('绝对值误差','fontsize', fontsize_baselinet),...
grid on, set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),  view(-36,15),colorbar('eastoutside'),caxis([0,0.25])


 end

%% Function used in this programm

function g=choose_function(funct_example,L)

switch funct_example

    case 1 % test exactness hyperinterpolation
        nexp=L;
        c0=rand(1); c1=rand(1); c2=rand(1); c3=rand(1);
        g=@(x,y,z) (c0+c1*x+c2*y+c3*z).^nexp;

    case 2 % test exactness filt. hyperinterpolation
        nexp=max(floor(L/2)-1,0);
        c0=rand(1); c1=rand(1); c2=rand(1); c3=rand(1);
        g=@(x,y,z) (c0+c1*x+c2*y+c3*z).^nexp;

    case 3 % function of that type

        funct_example_sub=0;

        switch funct_example_sub
            case 0 % Test function used in the accompanying paper.
                g=@(x,y,z) exp(-1./(x.^2+y.^2+z.^2));
                % fstring='exp(-1./(x.^2+y.^2+z.^2))';
            case 1
                g=@(x,y,z) (1-x.^2-y.^2-z.^2).*exp(x.*cos(y));
                % fstring='(1-x.^2-y.^2-z.^2).*exp(x.*cos(y))';
            case 2
                g=@(x,y,z) exp((x.^6).*cos(y+2*z));
                % fstring='exp((x.^6).*cos(y+2*z))';
        end

end        
end

function [pXYZR_3D, abs_err_3D] = evaluate_cube(f,coeff,n)

    TT = -1:0.1:1;
     
    [XR_3, YR_3, ZR_3] = meshgrid(TT);

    fXYZR_3D = feval(f, XR_3, YR_3, ZR_3);

    [V,dbox,duples]=dCHEBVAND0(n,[XR_3(:) YR_3(:) ZR_3(:)]);

    pXYZR0 = V*coeff;

    pXYZR_3D = reshape(pXYZR0,size(XR_3,1),size(YR_3,1),size(ZR_3,1));

    abs_err_3D = abs(fXYZR_3D-pXYZR_3D);
end