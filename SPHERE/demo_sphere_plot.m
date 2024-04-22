function demo_sphere_plot
% Codes based on Alvise Sommariva (University of Padova)
% Date: 17 August, 2023
% demo_interval_plot is used to plot original function
% and its recovery function via different hyperinterpolants
% on the unit-sphere

LV=15;     % Hyperinterpolant maximum degree.
NV=30; %NV=2*LV;   % Degree of the rule used in hyperinterpolation.
NR=50;     % Degree of precision of the reference rule (estimate L2 error).


% * Function to approximate:
% 1. degree L poly., 2. degree floor(L/2)-1 poly. 3. test functions
% (see line 228 approx).
funct_example=3;

% ........ Numerical approximation, varying the degree in "N" ............

AEinfMV=[]; AE2MV=[]; beta0MV=[]; % vectors used for statistics

% Define quadrature rule for hyperinterpolation at degree N.
XYZW=cub_sphere(NV); X=XYZW(:,1); Y=XYZW(:,2); Z=XYZW(:,3); W=XYZW(:,4);

% Define quadrature rule for L2 error at degree NR.
XYZWR=cub_sphere(NR); XR=XYZWR(:,1); YR=XYZWR(:,2); ZR=XYZWR(:,3);  WR=XYZWR(:,4);

% Vandermonde matrix at nodes and polynomial degrees.
[V,degs]=vandermonde_sphharm(LV,[X Y Z]);

% define function (see attached file at the bottom)
g=choose_function(funct_example,LV);

% ... evaluate function to approximate ...
gXYZ=feval(g,X,Y,Z);

a = 0.02; sigma = 0.02;

pert_impulse=a*(1-2*rand(length(gXYZ),1))*binornd(1,0.5);

var=sigma^2; pert_gauss=sqrt(var)*randn(size(gXYZ));

% add gaussian + impulse noise
pert=pert_impulse+pert_gauss;

% perturbed values
gXYZ_pert=gXYZ+pert;

% ... determine polynomial hyperinterpolant ...
coeff0=(gXYZ_pert.*W)'*V; coeff0=coeff0';

% ... determine polynomial hyperinterpolant (noise) ...
noise0=(pert.*W)'*V; noise0=noise0';


lambdas=sort(abs(coeff0),'descend');
lambdaL=lambdas(9);

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
            [VR,degs]=vandermonde_sphharm(LV,[XR YR,ZR]);
            pXYZR(:,ktest)=VR*coeff;

            % errors
            AEinfV(ktest)=norm(gXYZR-pXYZR(:,ktest),inf); % absolute error (inf norm)
            AE2V(ktest)=sqrt(WR'*((gXYZR-pXYZR(:,ktest)).^2)); % absolute error (2 norm)
            beta0V(ktest)=sum(abs(coeff) > 0);
end

%% Plot

[Fmax, imax] = max(gXYZR);
[Fmin, imin] = min(gXYZR);
scale = -0.35;
scale1 = -scale;
FS = 1 + (scale/(Fmax-Fmin))*(gXYZR-Fmin);
FS1 = 1 + (scale1/(Fmax-Fmin))*(gXYZR-Fmin);
fontsize_baselinet = 35;

tri = convhull([XR YR ZR]);

% exact function and its noisy function
axes('position',[0.025,0.55,0.2,0.47])
fg = trisurf(tri,XR.*FS, YR.*FS, ZR.*FS, gXYZR,'facecolor','interp');
set(fg,'EdgeColor', 'none');
title('目标函数','fontsize',fontsize_baselinet,'position',[0,0,1.25])
colormap(jet(255));
view(-25,25), axis vis3d, axis equal tight, colorbar('south'), %caxis([1.32,1.47])
axis off

tri1 = convhull([X Y Z]);

axes('position',[0.025 0.05 0.2 0.47])
[Fmax, imax] = max(gXYZ_pert);
[Fmin, imin] = min(gXYZ_pert);

FS = 1 + (scale/(Fmax-Fmin))*(gXYZ_pert-Fmin);
fg = trisurf(tri1,X.*FS, Y.*FS, Z.*FS, gXYZ_pert,'facecolor','interp');
set(fg,'EdgeColor', 'none');
title('含噪声的目标函数','fontsize',fontsize_baselinet,'position',[0,0,1.45])
colormap(jet(255));
view(-25,25), axis vis3d, axis equal tight, colorbar('south'),%caxis([1.32,1.47]),
axis off

% Hyper. and its error

axes('position',[0.275 0.55 0.2 0.47])
[Fmax, imax] = max(pXYZR(:,5));
[Fmin, imin] = min(pXYZR(:,5));

FS = 1 + (scale/(Fmax-Fmin))*(pXYZR(:,6)-Fmin);
fg = trisurf(tri,XR.*FS, YR.*FS, ZR.*FS, pXYZR(:,6),'facecolor','interp');
set(fg,'EdgeColor', 'none');
title('超插值','fontsize',fontsize_baselinet,'position',[0,0,1.25])
colormap(jet(255));
view(-25,25), axis vis3d, axis equal tight, colorbar('south'),caxis([min(gXYZR),max(gXYZR)]),
axis off

axes('position',[0.275,0.05,0.2,0.47])
[Fmax, imax] = max(abs(pXYZR(:,6)-gXYZR));
[Fmin, imin] = min(abs(pXYZR(:,6)-gXYZR));

FS1 = 1 + (scale1/(Fmax-Fmin))*(abs(pXYZR(:,6)-gXYZR)-Fmin);
fg = trisurf(tri,XR.*FS1, YR.*FS1, ZR.*FS1, abs(pXYZR(:,6)-gXYZR),'facecolor','interp');
set(fg,'EdgeColor', 'none');
title('绝对值误差','fontsize',fontsize_baselinet,'position',[0,0,1.45])
colormap(jet(255));
view(-25,25), axis vis3d, axis equal tight, colorbar('south'),caxis([min(abs(pXYZR(:,4)-gXYZR)),max(abs(pXYZR(:,4)-gXYZR))]),
axis off

% Hard hyper. and its error
%subplot(4,4,[9,10,13,14])
axes('position',[0.525 0.55 0.2 0.47])
[Fmax, imax] = max(pXYZR(:,5));
[Fmin, imin] = min(pXYZR(:,5));

FS = 1 + (scale/(Fmax-Fmin))*(pXYZR(:,5)-Fmin);
fg = trisurf(tri,XR.*FS, YR.*FS, ZR.*FS, pXYZR(:,5),'facecolor','interp');
set(fg,'EdgeColor', 'none');
title('硬阈值超插值','fontsize',fontsize_baselinet,'position',[0,0,1.25])
colormap(jet(255));
view(25,25), axis vis3d, axis equal tight, colorbar('south'),caxis([min(gXYZR),max(gXYZR)]),
axis off

axes('position',[0.525,0.05,0.2,0.47])
[Fmax, imax] = max(abs(pXYZR(:,5)-gXYZR));
[Fmin, imin] = min(abs(pXYZR(:,5)-gXYZR));

FS1 = 1 + (scale1/(Fmax-Fmin))*(abs(pXYZR(:,5)-gXYZR)-Fmin);
fg = trisurf(tri,XR.*FS1, YR.*FS1, ZR.*FS1, abs(pXYZR(:,5)-gXYZR),'facecolor','interp');
set(fg,'EdgeColor', 'none');
title('绝对值误差','fontsize',fontsize_baselinet,'position',[0,0,1.45])
colormap(jet(255));
view(25,25), axis vis3d, axis equal tight, colorbar('south'),caxis([min(abs(pXYZR(:,4)-gXYZR)),max(abs(pXYZR(:,4)-gXYZR))]),
axis off

% hybrid hyper. and its error
axes('position',[0.775 0.55 0.2 0.47])    
[Fmax, imax] = max(pXYZR(:,4));
[Fmin, imin] = min(pXYZR(:,4));

FS = 1 + (scale/(Fmax-Fmin))*(pXYZR(:,4)-Fmin);
fg = trisurf(tri,XR.*FS, YR.*FS, ZR.*FS, pXYZR(:,4),'facecolor','interp');
set(fg,'EdgeColor', 'none');
title('混合超插值','fontsize',fontsize_baselinet,'position',[0,0,1.25]),
colormap(jet(255));
view(25,25), axis vis3d, axis equal tight, colorbar('south'),caxis([min(gXYZR),max(gXYZR)]),
axis off

axes('position',[0.775,0.05,0.2,0.47])    
[Fmax, imax] = max(abs(pXYZR(:,4)-gXYZR));
[Fmin, imin] = min(abs(pXYZR(:,4)-gXYZR));

FS1 = 1 + (scale1/(Fmax-Fmin))*(abs(pXYZR(:,4)-gXYZR)-Fmin);
fg = trisurf(tri,XR.*FS1, YR.*FS1, ZR.*FS1, abs(pXYZR(:,4)-gXYZR),'facecolor','interp');
set(fg,'EdgeColor', 'none');
title('绝对值误差','fontsize',fontsize_baselinet,'position',[0,0,1.45]),
colormap(jet(255));
view(25,25), axis vis3d, axis equal tight, colorbar('south'),caxis([min(abs(pXYZR(:,4)-gXYZR)),max(abs(pXYZR(:,4)-gXYZR))]),
axis off

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

        funct_example_sub=3;

        switch funct_example_sub
            case 0
                g=@(x,y,z) exp(-1./(x.^2+y.^2+z.^2));
                % fstring='exp(-1./(x.^2+y.^2+z.^2))';
            case 1
                g=@(x,y,z) (1-x.^2-y.^2-z.^2).*exp(x.*cos(y));
                % fstring='(1-x.^2-y.^2-z.^2).*exp(x.*cos(y))';
            case 2
                g=@(x,y,z) exp((x.^6).*cos(y+2*z));
                % fstring='exp((x.^6).*cos(y+2*z))';
            case 3 % Test function used in the accompanying paper.
                delta2=9*gamma(5/2)/(2*gamma(3));
                %delta2=12*gamma(7/2)/(2*gamma(4));
                phi_tilde=@(r) (((1-r).*(1-r > 0)).^6).*(35*r.^2+18*r+3);
                r1=@(x,y,z) ((x-1).^2+y.^2+z.^2)/delta2; %z1=[1 0 0]
                r2=@(x,y,z) ((x+1).^2+y.^2+z.^2)/delta2; %z2=[-1 0 0]
                r3=@(x,y,z) (x.^2+(y-1).^2+z.^2)/delta2; %z3=[0 1 0]
                r4=@(x,y,z) (x.^2+(y+1).^2+z.^2)/delta2; %z4=[0 -1 0]
                r5=@(x,y,z) (x.^2+y.^2+(z-1).^2)/delta2; %z5=[0 0 1]
                r6=@(x,y,z) (x.^2+y.^2+(z+1).^2)/delta2; %z6=[0 0 -1]
                g=@(x,y,z) (1/3)*(phi_tilde(r1(x,y,z))+...
                    phi_tilde(r2(x,y,z))+phi_tilde(r3(x,y,z))+...
                    phi_tilde(r4(x,y,z))+phi_tilde(r5(x,y,z))+...
                    phi_tilde(r6(x,y,z)));

        end


end

end
