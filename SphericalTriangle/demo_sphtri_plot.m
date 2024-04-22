function demo_sphtri_plot

% This function is used to plot the denoising effect of hyper. and its 
% variants over the spherical triangle with vertices A=[1,0,0], B=[0,1,0]
% and C=[0,0,1].
% Date: 20 Sep, 2023
% Codes based on Alvise Sommariva (University of Padova)

clear; clf;

domain_example=0;

%--------------------------------------------------------------------------
% Degrees of precision in numerical experiments: can be a vector.
%--------------------------------------------------------------------------
%nV=1:10;
LV=8;
NV=11;
NR=30;
%--------------------------------------------------------------------------
% Noise and lasso parameter.
%--------------------------------------------------------------------------

noise=1;          % 0: no noise, 1: noise
a=0.1;              % defining impulse noise (in experiment 2)
sigma=0.2; %0.02;   % defining gaussian noise (in experiment 2)
%lambda=10^(-2);   % defining lasso parameter
pos=0;            % extraction type.
domain_structure.domain='spherical-triangle';

% ....... Special settings .......

% Approximation type parameter "pts_type".
%     case 1, pts_type='Hyperinterpolation full set';
%     case 2, pts_type='Hyperinterpolation compressed set';
%
% Note: the case "2" should be used for mild values of "n", say at most 15.
pts_type=1;

% Plot domain and nodes: do_plot=1 (yes), do_plot=0 (no).
do_plot=1;

% testing functions:
% 1. polynomial of degree n,
% 2. polynomial of degree floor(n/2)-1
% 3. gaussian like exponential
%funct_example=7;

funct_example=7;

% ....... Apply settings to define domain, pointsets, and functions .......

% Domain
vertices=define_domain(domain_example);

% Test points
%[XYZWR,dbox]=define_cub_rule(domain_structure,30);
P1=vertices(1,:); P2=vertices(2,:); P3=vertices(3,:);
XYZWR = cub_sphtri(NR,P1',P2',P3',pos);

XR=XYZWR(:,1:end-1); WR=XYZWR(:,end);

% ........ Numerical approximation, varying the degree in "nV" ............
fprintf('\n \t ');
AEinfMV=[]; AE2MV=[]; beta0MV=[]; % vectors used for statistics

n = NV;
dimpoly=(n+1)^2;

% ... extract hyperinterpolation set (notice that ade=2*n) ...
if pts_type == 2 % compressed set
    % fprintf('\n \t * Compressed set')
    XYZW = cub_sphtri(2*n,P1',P2',P3',pos);
    [pts,weights,momerr,dbox] =...
        dCATCH(2*n,XYZW(:,1:3),XYZW(:,4));
else % full set
    % fprintf('\n \t * Full set')
    XYZW = cub_sphtri(2*n,P1',P2',P3',pos);
    pts=XYZW(:,1:end-1); weights=XYZW(:,end);
end

% .. testing AE_L2err hyperinterpolation error for each "f" at "deg" ..
g = define_function(funct_example);

% ... evaluate function to approximate ...
g_pts=feval(g,pts(:,1),pts(:,2),pts(:,3));

% ... Add noise (if present) ...
[g_pts_pert,pert] = add_noise(noise,a,sigma,g_pts);

% ... determine polynomial hyperinterpolant ...
[coeff0,R,jvec,dbox,degs] = dHYPERFIT2(LV,pts,weights,g_pts_pert,...
    [],[],domain_structure,dimpoly);

if iscell(jvec), degs=degs(jvec{1}); else, degs=degs(jvec); end

lambdas=sort(abs(coeff0),'descend');
lambdaL=lambdas(6);

%LV=NV;

% test hyperinterpolant with or withour filters.
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
            parms.w=weights;
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

    % ... evaluate hyperinterpolant at initial pointset ...
    p_XR(:,ktest)=dPOLYVAL2(LV,coeff,XR,R,jvec,dbox,domain_structure,dimpoly);

    % ... estimating hyperinterpolant error ...
    g_XR=feval(g,XR(:,1),XR(:,2),XR(:,3));

    AEinfV(ktest)=norm(g_XR-p_XR(:,ktest),inf); % absolute error (inf norm)
    AE2V(ktest)=sqrt(WR'*((g_XR-p_XR(:,ktest)).^2)); % absolute error (2 norm)
    beta0V(ktest)=sum(abs(coeff) > 0);

end


%% Plot denoising effect
figure(1)

fontsize_baselinet = 35;

% C. compute triangulation (valid only for sph.poly. lying in some hemisphere)

% 1. rotate to North Pole

CC=mean(vertices); CC=CC/norm(CC);

% ................ rotation matrix centroid to north pole .................

[az,el,r] = cart2sph(CC(1),CC(2),CC(3));
phi=az; theta=pi/2-el;
cp=cos(phi); sp=sin(phi); ct=cos(theta); st=sin(theta);
R1=[ct 0 -st; 0 1 0; st 0 ct]; R2=[cp sp 0; -sp cp 0; 0 0 1];
rotmat=R1*R2; inv_rotmat=rotmat';

% ........................ rotate vertices to north pole...................

vertices_NP=(rotmat*vertices')';

% ................... stereographic map from south pole ...................

% ....... vertices .......

XX_SP=vertices_NP(:,1); YY_SP=vertices_NP(:,2); ZZ_SP=vertices_NP(:,3);
rat=1./(1+ZZ_SP);

XX_SPm=rat.*XX_SP; YY_SPm=rat.*YY_SP; % vertices on the plane

% ....... points .......

rat2=1./(1+XR(:,3));
X_SP=rat2.*XR(:,1); % points on the plane
Y_SP=rat2.*XR(:,2);

% ...... triangulation on the plane ....

tri = delaunay(X_SP,Y_SP);

% exact function and its noisy function
axes('position',[0.025,0.55,0.2,0.47])
fg = trisurf(tri,XR(:,1), XR(:,2), XR(:,3), g_XR,'facecolor','interp');
set(fg,'EdgeColor', 'none');
title('目标函数','fontsize',fontsize_baselinet,'position',[0,0,1])
colormap(jet(255));
view(125,8), axis vis3d, axis equal tight, colorbar('south'), %caxis([min(g_XR)-0.02,max(g_XR)+0.02])
axis off

% ....... points for noisy function .......

rat21=1./(1+pts(:,3));
X_SP1=rat21.*pts(:,1); % points on the plane
Y_SP1=rat21.*pts(:,2);

% ...... triangulation on the plane ....

tri1 = delaunay(X_SP1,Y_SP1);

axes('position',[0.025 0.05 0.2 0.47])
fg = trisurf(tri1, pts(:,1), pts(:,2), pts(:,3), g_pts_pert,'facecolor','interp');
set(fg,'EdgeColor', 'none');
title('含噪声的目标函数','fontsize',fontsize_baselinet,'position',[0,0,1])
colormap(jet(255));
view(125,8), axis vis3d, axis equal tight, colorbar('south'),%caxis([min(g_pts_pert)-0.02,max(g_pts_pert)+0.02]) %caxis([-3.0,3.0]),
axis off

% Hyper. and its error

axes('position',[0.275 0.55 0.2 0.47])

fg = trisurf(tri, XR(:,1), XR(:,2), XR(:,3), p_XR(:,6),'facecolor','interp');
set(fg,'EdgeColor', 'none');
title('超插值','fontsize',fontsize_baselinet,'position',[0,0,1])
colormap(jet(255));
view(125,8), axis vis3d, axis equal tight, colorbar('south'),caxis([min(p_XR(:,5)),max(p_XR(:,5))]),
axis off

axes('position',[0.275,0.05,0.2,0.47])
fg = trisurf(tri, XR(:,1), XR(:,2), XR(:,3), abs(p_XR(:,6)-g_XR),'facecolor','interp');
set(fg,'EdgeColor', 'none');
title('绝对值误差','fontsize',fontsize_baselinet,'position',[0,0,1])
colormap(jet(255));
view(125,8), axis vis3d, axis equal tight, colorbar('south'),caxis([min(abs(p_XR(:,4)-g_XR)),max(abs(p_XR(:,4)-g_XR))])
axis off

% Hard hyper. and its error
axes('position',[0.525 0.55 0.2 0.47])
fg = trisurf(tri, XR(:,1), XR(:,2), XR(:,3), p_XR(:,5),'facecolor','interp');
set(fg,'EdgeColor', 'none');
title('硬阈值超插值','fontsize',fontsize_baselinet,'position',[0,0,1])
colormap(jet(255));
view(125,8), axis vis3d, axis equal tight, colorbar('south'),caxis([min(p_XR(:,5)),max(p_XR(:,5))]),
axis off

axes('position',[0.525,0.05,0.2,0.47])
fg = trisurf(tri, XR(:,1), XR(:,2), XR(:,3), abs(p_XR(:,5)-g_XR),'facecolor','interp');
set(fg,'EdgeColor', 'none');
title('绝对值误差','fontsize',fontsize_baselinet,'position',[0,0,1])
colormap(jet(255));
view(125,8), axis vis3d, axis equal tight, colorbar('south'),caxis([min(abs(p_XR(:,4)-g_XR)),max(abs(p_XR(:,4)-g_XR))])
axis off

% Hybrid hyper. and its error
%subplot(4,4,[11,12,15,16])
axes('position',[0.775 0.55 0.2 0.47])    
fg = trisurf(tri, XR(:,1), XR(:,2), XR(:,3), p_XR(:,4),'facecolor','interp');
set(fg,'EdgeColor', 'none');
title('混合超插值','fontsize',fontsize_baselinet,'position',[0,0,1])
colormap(jet(255));
view(125,8), axis vis3d, axis equal tight, colorbar('south'),caxis([min(p_XR(:,5)),max(p_XR(:,5))]),
axis off

axes('position',[0.775,0.05,0.2,0.47])    
fg = trisurf(tri, XR(:,1), XR(:,2), XR(:,3), abs(p_XR(:,4)-g_XR),'facecolor','interp');
set(fg,'EdgeColor', 'none');
title('绝对值误差','fontsize',fontsize_baselinet,'position',[0,0,1])
colormap(jet(255));
view(125,8), axis vis3d, axis equal tight, colorbar('south'),caxis([min(abs(p_XR(:,4)-g_XR)),max(abs(p_XR(:,4)-g_XR))])
axis off


% [Fmax, imax] = max(g_XR);
% [Fmin, imin] = min(g_XR);
% scale = 0.5;
% FS = 1 + (scale/(Fmax-Fmin))*(g_XR-Fmin);
% 
% 
% tri = convhull([XR(:,1) XR(:,2) XR(:,3)]);
% 
% % exact function and its noisy function
% axes('position',[0.025,0.55,0.2,0.47])
% fg = trisurf(tri,XR(:,1).*FS, XR(:,2).*FS, XR(:,3).*FS, g_XR,'facecolor','interp');
% set(fg,'EdgeColor', 'none');
% title('\textbf{Exact function}','interpreter','latex','fontsize',fontsize_baselinet,'position',[0,0,1.1])
% colormap(jet(255));
% view(125,18), axis vis3d, axis equal tight, colorbar('south'), caxis([min(g_XR)-0.02,max(g_XR)+0.02])
% axis off
% 
% tri1 = convhull([pts(:,1) pts(:,2) pts(:,3)]);
% 
% axes('position',[0.025 0.05 0.2 0.47])
% [Fmax, imax] = max(g_pts_pert);
% [Fmin, imin] = min(g_pts_pert);
% 
% FS = 1 + (scale/(Fmax-Fmin))*(g_pts_pert-Fmin);
% fg = trisurf(tri1, pts(:,1).*FS, pts(:,2).*FS, pts(:,3).*FS, g_pts_pert,'facecolor','interp');
% set(fg,'EdgeColor', 'none');
% title('\textbf{Noisy function}','interpreter','latex','fontsize',fontsize_baselinet,'position',[0,0,1.2])
% colormap(jet(255));
% view(125,38), axis vis3d, axis equal tight, colorbar('south'),caxis([min(g_pts_pert)-0.02,max(g_pts_pert)+0.02]) %caxis([-3.0,3.0]),
% axis off
% 
% % Hard thresholding hyper. and its error
% 
% axes('position',[0.275 0.55 0.2 0.47])
% [Fmax, imax] = max(p_XR(:,5));
% [Fmin, imin] = min(p_XR(:,5));
% 
% FS = 1 + (scale/(Fmax-Fmin))*(p_XR(:,5)-Fmin);
% fg = trisurf(tri, XR(:,1).*FS, XR(:,2).*FS, XR(:,3).*FS, p_XR(:,5),'facecolor','interp');
% set(fg,'EdgeColor', 'none');
% title('\textbf{Hard thresholding hyper.}','interpreter','latex','fontsize',fontsize_baselinet,'position',[0,0,1.1])
% colormap(jet(255));
% view(125,18), axis vis3d, axis equal tight, colorbar('south'),caxis([min(p_XR(:,5))-0.02,max(p_XR(:,5))+0.02])%caxis([1.32,1.47])
% axis off
% 
% axes('position',[0.275,0.05,0.2,0.47])
% [Fmax, imax] = max(abs(p_XR(:,5)-g_XR));
% [Fmin, imin] = min(abs(p_XR(:,5)-g_XR));
% 
% FS = 1 + (scale/(Fmax-Fmin))*(abs(p_XR(:,5)-g_XR)-Fmin);
% fg = trisurf(tri, XR(:,1).*FS, XR(:,2).*FS, XR(:,3).*FS, abs(p_XR(:,5)-g_XR),'facecolor','interp');
% set(fg,'EdgeColor', 'none');
% title('\textbf{Error}','interpreter','latex','fontsize',fontsize_baselinet,'position',[0,0,1.3])
% colormap(jet(255));
% view(125,38), axis vis3d, axis equal tight, colorbar('south'),caxis([min(abs(p_XR(:,5)-g_XR))-0.02,max(abs(p_XR(:,5)-g_XR))+0.02])
% axis off
% 
% % Lasso hyper. and its error
% %subplot(4,4,[9,10,13,14])
% axes('position',[0.525 0.55 0.2 0.47])
% [Fmax, imax] = max(p_XR(:,3));
% [Fmin, imin] = min(p_XR(:,3));
% 
% FS = 1 + (scale/(Fmax-Fmin))*(p_XR(:,3)-Fmin);
% fg = trisurf(tri, XR(:,1).*FS, XR(:,2).*FS, XR(:,3).*FS, p_XR(:,3),'facecolor','interp');
% set(fg,'EdgeColor', 'none');
% title('\textbf{Lasso hyper.}','interpreter','latex','fontsize',fontsize_baselinet,'position',[0,0,1.1])
% colormap(jet(255));
% view(125,18), axis vis3d, axis equal tight, colorbar('south'),caxis([min(p_XR(:,3))-0.02,max(p_XR(:,3))+0.02])%caxis([1.32,1.47])
% axis off
% 
% axes('position',[0.525,0.05,0.2,0.47])
% [Fmax, imax] = max(abs(p_XR(:,3)-g_XR));
% [Fmin, imin] = min(abs(p_XR(:,3)-g_XR));
% 
% FS = 1 + (scale/(Fmax-Fmin))*(abs(p_XR(:,3)-g_XR)-Fmin);
% fg = trisurf(tri, XR(:,1).*FS, XR(:,2).*FS, XR(:,3).*FS, abs(p_XR(:,3)-g_XR),'facecolor','interp');
% set(fg,'EdgeColor', 'none');
% title('\textbf{Error}','interpreter','latex','fontsize',fontsize_baselinet,'position',[0,0,1.3])
% colormap(jet(255));
% view(125,8), axis vis3d, axis equal tight, colorbar('south'),caxis([min(abs(p_XR(:,3)-g_XR))-0.02,max(abs(p_XR(:,3)-g_XR))+0.02])
% axis off
% 
% % Classical hyper. and its error
% %subplot(4,4,[11,12,15,16])
% axes('position',[0.775 0.55 0.2 0.47])    
% [Fmax, imax] = max(p_XR(:,6));
% [Fmin, imin] = min(p_XR(:,6));
% 
% FS = 1 + (scale/(Fmax-Fmin))*(p_XR(:,6)-Fmin);
% fg = trisurf(tri, XR(:,1).*FS, XR(:,2).*FS, XR(:,3).*FS, p_XR(:,6),'facecolor','interp');
% set(fg,'EdgeColor', 'none');
% title('\textbf{Hyperinterpolation}','interpreter','latex','fontsize',fontsize_baselinet,'position',[0,0,1.1])
% colormap(jet(255));
% view(125,18), axis vis3d, axis equal tight, colorbar('south'),caxis([min(p_XR(:,6))-0.02,max(p_XR(:,6))+0.02])%caxis([1.32,1.47])
% axis off
% 
% axes('position',[0.775,0.05,0.2,0.47])    
% [Fmax, imax] = max(abs(p_XR(:,6)-g_XR));
% [Fmin, imin] = min(abs(p_XR(:,6)-g_XR));
% 
% FS = 1 + (scale/(Fmax-Fmin))*(abs(p_XR(:,6)-g_XR)-Fmin);
% fg = trisurf(tri, XR(:,1).*FS, XR(:,2).*FS, XR(:,3).*FS, abs(p_XR(:,6)-g_XR),'facecolor','interp');
% set(fg,'EdgeColor', 'none');
% title('\textbf{Error}','interpreter','latex','fontsize',fontsize_baselinet,'position',[0,0,1.35])
% colormap(jet(255));
% view(125,38), axis vis3d, axis equal tight, colorbar('south'),caxis([min(abs(p_XR(:,6)-g_XR))-0.02,max(abs(p_XR(:,6)-g_XR))+0.02])
% axis off


end



%% Function used in this program

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

function vertices=define_domain(example)


switch example

    case 0
        % large domain
        vertices=[ 1 0 0;
            0 1 0;
            0 0 1];

    case 1
        a=0.5; b=0.1; % MEDIUM-LARGE SIZE (Australia)
        vertices=[0 0 1;
            0 sqrt(a) sqrt(1-a);
            sqrt(b) sqrt(a/2) sqrt(1-a/2-b)];

end

end


function g = define_function(funct_example)

switch funct_example

    case 1 % test exactness hyperinterpolation
        nexp=n;
        c0=rand(1); c1=rand(1); c2=rand(1); c3=rand(1);
        g=@(x,y,z) (c0+c1*x+c2*y+c3*z).^nexp;

    case 2 % test exactness filt. hyperinterpolation
        nexp=floor(n/2);
        c0=rand(1); c1=rand(1); c2=rand(1); c3=rand(1);
        g=@(x,y,z) (c0+c1*x+c2*y+c3*z).^nexp;

    case 3 % exponential type
        c0=rand(1); c1=rand(1); c2=rand(1); c3=rand(1);
        g=@(x,y,z) exp(-(c0*x.^2+c1*y.^2+c2*z.^2+c3));

    case 4
        g=@(x,y,z) 1+x+y.^2+x.^2.*y+x.^4+y.^5+x.^2.*y.^2.*z.^2;

    case 5
        g=@(x,y,z) cos(10*(x+y+z));

    case 6
        g=@(x,y,z) exp(-(x.^2+y.^2+(z-1).^2));

    case 7
         g=@(x,y,z) exp(-((x-1/sqrt(3)).^2 + (y-1/sqrt(3)).^2+ (z-1/sqrt(3)).^2 ) );
    case 8 % Test function used in the accompanying paper.
        g=@(x,y,z) 0.75*exp(-((9*x-2).^2)/4- ((9*y-2).^2)/4- ((9*z-2).^2)/4)+...
                   0.75*exp(-((9*x+1).^2)/49- ((9*y+1).^2)/10- ((9*z+1).^2)/10)+...
                   0.5*exp(-((9*x-7).^2)/4-((9*y-3).^2)/4-((9*z-5).^2)/4)-...
                   0.2*exp(-((9*x-4).^2) - ((9*y-7).^2)-((9*z-5).^2));

end
end

function [g_pts_pert,pert] = add_noise(noise,a,sigma,g_pts);
if noise
    % add impulse noise
    pert_impulse=0;
    if a > 0
        pert_impulse=a*(1-2*rand(length(g_pts),1))*binornd(1,0.5);
        while norm(pert_impulse) == 0
            pert_impulse=a*(1-2*rand(length(g_pts),1))*binornd(1,0.5);
        end
    end

    % add gaussian noise
    pert_gauss=0;
    if sigma > 0
        var=sigma^2;
        pert_gauss=sqrt(var)*randn(size(g_pts));
        while norm(pert_gauss) == 0
            pert_gauss=sqrt(var)*randn(size(g_pts));
        end
    end

    % add gaussian + impulse noise
    pert=pert_impulse+pert_gauss;

    % perturbed values
    g_pts_pert=g_pts+pert;
else
    g_pts_pert=g_pts;
end
end
