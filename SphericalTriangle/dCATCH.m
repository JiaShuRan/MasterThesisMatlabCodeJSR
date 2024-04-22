function [nodes,w,momerr,dbox] = dCATCH(deg,X,u,LHDM_options,verbose,...
    dim_poly)

%--------------------------------------------------------------------------
% Object:
% This routine implements the Caratheodory-Tchakaloff d-variate discrete
% measure compression.
% Examples of its application are probability measures (designs) or
% quadrature formulas.
% Moments are invariant (close to machine precision) up to degree "deg"
% and adapt to the (numerical) dimension of the polynomial space on the
% point set X.
% The routine works satisfactorily for low/moderate degrees, depending on 
% the dimension.
%--------------------------------------------------------------------------
% Input:
% deg: polynomial exactness degree;
% X: d-column array of point coordinates;
% * u: 1-column array of nonnegative weights, or nonnegative scalar in
%    case of equal weights;
% * LHDM_options: structure containing the values of optimization 
%       parameters (see "LHDM.m" for details); 
% * verbose: 0: no information from routine;
%            1: relevant information from routine;
% * dim_poly: dimension of polynomial space (if known in advance, useful
%   for instance on the sphere or its portions).
% Note: the variables with an asterisk "*" are not mandatory and can be 
% also set as empty matrix.
%--------------------------------------------------------------------------
% Output:
% nodes: d-column array of extracted mass points coordinates; the variable
%    "nodes" has rows that are also rows of "X", i.e. nodes represent a 
%    subset of "X";
% w: 1-column array of corresponding new positive weights;
% momerr: moment reconstruction error;
% * dbox: variable that defines a hyperrectangle with 
%    sides parallel to the axis, containing the domain (or pointset X in 
%    the discrete case). 
%    It is a matrix with dimension "2 x d", where "d" is the dimension of
%    the space in which it is embedded the domain. 
%    For instance, for a 2-sphere, it is "d=3", for a 2 dimensional  
%    polygon it is "d=2".
%    As example, the set "[-1,1] x [0,1]" is described as 
%                          "dbox=[-1 0; 1 1]".
%--------------------------------------------------------------------------
% Data:
% Written on 26/07/2020 by M. Dessole, F. Marcuzzi, M. Vianello.
% Last update by:
% 04/01/2020: A. Sommariva.
%--------------------------------------------------------------------------

% .........................  Function Body ................................

% ..... troubleshooting .....

if nargin < 6, dim_poly=[]; end
if nargin < 5, verbose=[]; end
if nargin < 4, LHDM_options=[]; end
if nargin < 3, u=[]; end

if isempty(verbose), verbose=0; end
if isempty(LHDM_options)
    dim = size(X,2);
    LHDM_options = struct( 'k', ceil(nchoosek(2*deg+dim,dim)/(deg*dim)),...
        'init', false, 'thres', 0.2222, 'thres_w', 0.8);
end
if isempty(u), u=1; end


% ..... Main code below .....

% Vandermonde-like matrix of a u-orthogonal polynomial basis on X
[U,~,~,~,dbox]=dORTHVAND2(deg,X,u,[],[],[],dim_poly);

if verbose
    fprintf('Vandermonde matrix size = %d x %d \n', size(U,1), size(U,2));
end

if size(U,1)<=size(U,2)
    if verbose
        fprintf('Vandermonde matrix not underdetermined: no compression');
        fprintf('\n');
    end
    % no compression expected
    nodes=X;
    
    if isscalar(u), u=u*ones(size(Q,1),1); end % weights saved as scalar
    w=u;
    momerr = 0;
    
else
    
    % further orthogonalization to reduce the conditioning
    [Q,~]=qr(U,0);
    
    % new moments
    if isscalar(u), u=u*ones(size(Q,1),1); end % weights saved as scalar
    orthmom=Q'*u;
    [nodes, w, momerr, ~]= NNLS(X, u, U, Q, orthmom, LHDM_options, verbose);
    
end










function [nodes, w, momerr, e]= NNLS(X, u, U, Q, orthmom, options, verbose)

%--------------------------------------------------------------------------
% Object:
% Caratheodory-Tchakaloff points and weights via accelerated NNLS
%--------------------------------------------------------------------------

tic;
if  isfield(options,'lsqnonneg')
    if options.lsqnonneg
        if verbose
            fprintf('Matlab lsqnonneg \n');
        end
        [weights,~,~,~,output] = lsqnonneg(Q',orthmom);
        iter = output.iterations;
        cardP = [];
    else
        [weights,~,~,iter]=LHDM(Q',orthmom,options,verbose);
    end
else
    [weights,~,~,iter]=LHDM(Q',orthmom,options,verbose);
end
e = toc;

% indexes of nonvanishing weights and compression
ind=find(abs(weights)>0);
nodes=X(ind,:);
w=weights(ind);

% moment reconstruction error
momerr=norm(U(ind,:)'*w-U'*u);

% displaying results
if verbose
    fprintf('NNLS number of outer iterations = %d \n', iter);
    fprintf('NNLS elapsed time = %.6f s  \n', e);
    fprintf('initial design cardinality = %4.0f \n',size(X,1));
    fprintf('concentrated support cardinality = %4.0f \n',length(w));
    fprintf('compression ratio = %4.0f \n',size(X,1)/length(w));
    fprintf('moment reconstruction error = %4.2e \n \n',momerr);
end



