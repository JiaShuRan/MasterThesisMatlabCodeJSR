
function plot_s2(domain_type,x,pts1,pts2,title_str,R)

%--------------------------------------------------------------------------
% Object:
% Plotting pointsets on regions of the sphere.
%--------------------------------------------------------------------------
% Input:
% domain_type: 'sphere','spherical-rectangle','spherical-triangle',
%       'spherical-caps'.
%
% x     : * spherical rectangle representation / spherical representation:
%         It is a 2 x 2 matrix, containing all (theta_1,theta2) such that
%         theta_1 in (x(1,1),x(1,2)) and theta_2 in (x(2,1),x(2,2)), and
%         the points in cartesian coordinates are
%                    x_1=sin(theta_1)*cos(theta_2),
%                    x_2=sin(theta_1)*sin(theta_2),
%                    x_3=cos(theta_1)
%         * spherical triangle representation:
%         It is a 3 x 3 matrix, where the k-th row contains the cartesian
%         coordinates of the k-th vertex of the domain.
%
% pts1: spherical/cartesian coordinates of first set of points on the unit
%      sphere possibly belonging to the spherical rectangle/triangle
%      defined by "x".
%      The k-th row of this matrix determines the k-th point.
%         The points in cartesian coordinates are
%                    x_1=sin(theta_1)*cos(theta_2),
%                    x_2=sin(theta_1)*sin(theta_2),
%                    x_3=cos(theta_1)
%
% pts2: spherical/cartesian coordinates of first set of points on the unit
%      sphere possibly belonging to the spherical rectangle/triangle
%      defined by "x".
%      The k-th row of this matrix determines the k-th point.
%      The points in cartesian coordinates are
%                    x_1=sin(theta_1)*cos(theta_2),
%                    x_2=sin(theta_1)*sin(theta_2),
%                    x_3=cos(theta_1)
%
% title_str: string with the title of the domain
%
% R: sphere radius
%--------------------------------------------------------------------------
% Dates:
% Written on 20/11/2020: A. Sommariva and M. Vianello;
%
% Modified on:
% 28/11/2020: A. Sommariva.
%--------------------------------------------------------------------------


% ......................... troubleshooting ...............................

if nargin < 2, x=[]; end
if nargin < 3, pts1=[]; end
if nargin < 4, pts2=[]; end
if nargin < 5, title_str=''; end
if nargin < 6, R=1; end

if isempty(x), x=[0 pi; 0 2*pi]; end
if isempty(pts1), pts1=[]; end
if isempty(pts2), pts2=[]; end
if isempty(title_str), title_str=''; end
if isempty(R), R=1; end


% ......................... plot sphere ...................................

% ... plot transparent sphere ...
axis tight;
pts_size1=3; pts_size2=6;
color1='black'; color2='red';

% ......................... plot pointset 1 ...............................
if isempty(pts1) == 0, plot_set(pts1,color1,'*',pts_size1); end

hold on;

% ......................... plot pointset 2 ...............................
if isempty(pts2) == 0, plot_set(pts2,color2,'o',pts_size2); end

switch domain_type
    case 'spherical-caps'
        fprintf('\n \t Plot cap \n \n');
        theta_intv=x(1,:);
        funx = @(theta,phi) sin(theta).*cos(phi);
        funy = @(theta,phi) sin(theta).*sin(phi);
        funz = @(theta,phi) cos(theta);
        fsurf(funx,funy,funz,[theta_intv(1) theta_intv(2) -pi pi])
        alpha 0.1
        
    case 'spherical-triangle'
        
        fprintf('\n \t Plot cap (sph-tri) \n \n');
        
        % .... plot cap ....
        
        theta(1)=acos(x(1,3)/norm(x(1,:)));
        theta(2)=acos(x(2,3)/norm(x(2,:)));
        theta(3)=acos(x(3,3)/norm(x(3,:)));
        theta_intv=[0 1.4*max(theta)];
        
        funx = @(theta,phi) sin(theta).*cos(phi);
        funy = @(theta,phi) sin(theta).*sin(phi);
        funz = @(theta,phi) cos(theta);
        fsurf(funx,funy,funz,[theta_intv(1) theta_intv(2) -pi pi]);
        alpha 0.1 % the lower, the more transparent the region
        
        % .... plot sph. triangle boundary ....
         A=x(1,:); B=x(2,:); C=x(3,:);
         plot_side(A,B); plot_side(B,C); plot_side(C,A);
        
         % .... hide ticks, labels, axis and background ....
         
         set(gca,'xtick',[]) % no ticks / labels
         set(gca,'xticklabel',[])
         set(gca,'ytick',[])
         set(gca,'yticklabel',[])
         set(gca,'ztick',[])
         set(gca,'zticklabel',[])
         
         set(gca,'color','none') % no background
         
         set(gca,'Visible','off') % no axis
         
    otherwise
        fprintf('\n \t Plot sphere \n \n');
        [X,Y,Z] = sphere(40);
        %plot3(R*X,R*Y,R*Z,'MarkerFaceColor','red');
        surf(R*X,R*Y,R*Z,'MarkerFaceColor','red');
        xlabel('X','FontSize',20,'FontWeight','bold')
        ylabel('Y','FontSize',20,'FontWeight','bold')
        zlabel('Z','FontSize',20,'FontWeight','bold')
        title('Spherical triangle','FontSize',25,'FontWeight','bold')
        grid on
        colorbar;
        view(115,20)
end

hold on;

% ... plot axis lines ...
Lfact=0.6;

% plot_arrow([0 0 0],[Lfact*R 0 0],'x');
% plot_arrow([0 0 0],[0 Lfact*R 0],'y');
% plot_arrow([0 0 0],[0 0 Lfact*R],'z');



% ......................... adding title ..................................
%title(title_str);
hold off;









function plot_set(pts,color_str,marker_type,marker_size)

%--------------------------------------------------------------------------
% Object:
% Plotting pointsets on regions of the sphere.
%--------------------------------------------------------------------------
% Input:
% pts: spherical or cartesian coordinates of a pointset on the
%       unit sphere.
%      The k-th row of this matrix determines the k-th point.
%      The points in cartesian coordinates are
%                    x_1=sin(theta_1)*cos(theta_2),
%                    x_2=sin(theta_1)*sin(theta_2),
%                    x_3=cos(theta_1)
%
% color_str: color of the pointset, e.g. color_str='magenta';
% marker_type: marker type of the pointset, e.g. marker_type='.';
% marker_size: marker size of the pointset, e.g. marker_size=10.
%--------------------------------------------------------------------------
% Dates:
% Written on 20/11/2020: A. Sommariva and M. Vianello;
%
% Modified on:
% 28/11/2020: A. Sommariva.
%--------------------------------------------------------------------------

% ......................... troubleshooting ...............................

if nargin < 1, pts=[]; end
if nargin < 2, color_str='magenta'; end
if nargin < 3, marker_type='o'; end
if nargin < 4, marker_size=10; end

if isempty(marker_type), marker_type='*'; end
if isempty(marker_size), marker_size=10; end


% ......................... plotting points ...............................

if size(pts,2) == 2 % spherical coordinates
    theta_1=pts(:,1); theta_2=pts(:,2);
    s1=sin(theta_1); c1=cos(theta_1); s2=sin(theta_2); c2=cos(theta_2);
    xx=s1.*c2; yy=s1.*s2; zz=c1;
else % cartesian coordinates
    xx=pts(:,1); yy=pts(:,2); zz=pts(:,3);
end

plot3(xx,yy,zz,marker_type,'Color',color_str,'MarkerSize',marker_size,...
     'MarkerEdgeColor',color_str,'MarkerFaceColor',color_str);









function plot_arrow(p1,p2,str)

%--------------------------------------------------------------------------
% Object:
% Plot a line joining the row vectors "p1" and "p2" with an arrow.
%--------------------------------------------------------------------------

dp = p2-p1;                         
quiver3(p1(1),p1(2),p1(3),dp(1),dp(2),dp(3),0,'LineWidth',2)
text(p1(1),p1(2),str)








function plot_side(A,B)

%--------------------------------------------------------------------------
% Object:
% Plot "geodesic line" between A and B.
%--------------------------------------------------------------------------
% Input:
% A,B: row vectors of two points of the sphere to be joined.
%--------------------------------------------------------------------------

R=norm(A);
t=linspace(0,1,100);
C=[];
for k=1:length(t)
    CC=A+t(k)*(B-A); CC=CC/norm(CC); C=[C; R*CC];
end
plot3(C(:,1),C(:,2),C(:,3),'k-','LineWidth',2);






function fill_sphtri(P1,P2,P3)

%--------------------------------------------------------------------------
% Object:
% Fill regions of the sphere.
%--------------------------------------------------------------------------
% Input:
% P1,P2,P3: row vectors.
%--------------------------------------------------------------------------

R=norm(P1);
N=10000;
vals=rand(N,3);
vals=vals./sum(vals,2);
P=[];
for k=1:N
    val=vals(k,:);
    PL=val(1)*P1+val(2)*P2+val(3)*P3;
    PL=R*PL/norm(PL);
    P=[P; PL];
end

% color_str='white';
plot3(P(:,1),P(:,2),P(:,3),'y*','Color',color_str,'MarkerSize',3,...
    'MarkerEdgeColor',color_str,'MarkerFaceColor',color_str);




