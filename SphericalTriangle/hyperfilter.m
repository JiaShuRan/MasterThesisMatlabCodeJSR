
function [coeff1,lambda]=hyperfilter(hypermode,coeff,degs,parms)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% This routine starting from hyperinterpolation coefficients computes those
% of between
% a) filtered hyperinterpolation (hypermode='filtered')
% b) lasso hyperinterpolation (hypermode='lasso')
% c) Tikhonov hyperinterpolation (hypermode='tikhonov')
% d) elastic hyperinterpolation (hypermode='elastic')
% e) hybrid hyperinterpolation (hypermode='hybrid')
%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
% hypermode: choose between 
%         'filtered', 'lasso', 'tikhonov', 'hard','hybrid'.
%  
% coeff: classic hyperinterpolation coefficients (approximation of Fourier
%        coefficients).
%
% degs:  if the coefficient "coeff(i)" refers to the polynomial "phi(i)" of
%        the orthonormal basis, then "degs(i)" is the total degree of 
%        "phi(i)".
% 
% parms: struct containing the following parameters:
%        parms.lambda: define "lambda" parameters, it may be
%                      a scalar (in the "lasso", "hybrid", "hard",case);
%                      a vector (in the "tikhonov"case);
%        parms.mu: "lasso" mu vector;
%        parms.b : vector in "tikhonov" and "elastic" case;
%        parms.pert: perturbation on data;
%        parms.w: cubature weights;
%--------------------------------------------------------------------------
% OUTPUT:
%--------------------------------------------------------------------------
% coeff1: computation of the "variant" hyperinterpolation coefficients.
%
% lambda: used lambda parameters for hyperinterpolation variants.
%--------------------------------------------------------------------------
% INFO:
%--------------------------------------------------------------------------
% Written by A. Sommariva.
% Last update: January 16, 2023.
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



% ................. Troubleshooting ................. 

if nargin < 4
    parms.lambda=10^(-1.5);
    parms.mu=ones(size(coeff));
    parms.b=ones(size(coeff));
    parms.w=1;
    parms.pert=0;
    parms.hybrid=0;
end

if isempty(parms)
    parms.lambda=10^(-1.5);
    parms.mu=ones(size(coeff));
    parms.b=ones(size(coeff));
    parms.w=1;
    parms.hybrid=0;
end

if length(parms.mu) == 1
    parms.mu=parms.mu*ones(size(coeff));
end

if length(parms.b) == 1
    parms.b=parms.b*ones(size(coeff));
end

if size(degs,2) > size(degs,1)
    degs=degs';
end

filt_type=1; % option: filtered hyperinterpolation function 
             % * see function "filter_hyps" attached to this file




% ................. Modify hyperinterpolation coefficients  ...............

switch hypermode

    case 'filtered'

        h=filter_hyps(filt_type);
        L=max(degs);
        coeff1=h(degs/L).*coeff;

        lambda=[];

    case 'lasso'

        lambda=parms.lambda; % regularization parameter
        mu=parms.mu;         % penalties parameters
        S=@(a,k) max(0,a-k)+min(0,a+k);
        L=length(degs);
        for l=1:L
            al=coeff(l);
            kl=lambda*mu(l);
            coeff1(l,1)=S(al,kl);
        end

    case 'tikhonov'

        % Regularization terms
        lambda=parms.lambda; % regularization parameters
        b=parms.b;

        % Apply Regularization
        term_in=1./(1+lambda*b.^2);
        coeff1=term_in.*coeff;
        

    case 'hard'

        % Regularization terms
        lambda=parms.lambda; % regularization parameters
        term=abs(coeff) > lambda;
        coeff1=term.*coeff;


    case 'hybrid'

        % "lambda_choice" parameter will follow this presets:
        %
        %      0: prechoosen lambda_1
        %      1: a priori parameter choice based on the weights and
        %         perturbation
        %      2: a priori parameter choice based on the weights and
        %         perturbation, and mu parameter
        %      3: Morozov (to be implemented)

       lambda=parms.lambda; % regularization parameters
       lambda1=lambda(1); 


        % filter type component
        h=filter_hyps(filt_type);
        L=max(degs);
        coeffF=h(degs/L);

        % lasso type component
        mu=parms.mu; % penalties parameters
        S=@(a,k) max(0,a-k)+min(0,a+k);
        L=length(degs);
        for l=1:L
            al=coeff(l);
            kl=lambda1*mu(l);
            coeffL(l,1)=S(al,kl);
        end

        coeff1=coeffF.*coeffL;

        lambda=lambda1;

    otherwise

        coeff1=coeff;
end




function h=filter_hyps(filt_type)

switch filt_type
    case 1
        h=@(x) (x >= 0 & x <= 1/2) + ((sin(pi*x)).^2).*(x > 0.5 & x < 1);
    otherwise
        h=@(x) (x >= 0 & x <= 1/2) + ((sin(pi*x)).^2).*(x > 0.5 & x < 1);
end

% function filt_coeffs=hyperfilter(hypermode,hyper_coeffs,degs,parms)
% 
% if nargin < 4 
%     parms.lambda=10^(-1.5);
%     parms.mu=ones(size(hyper_coeffs));
% end
% 
% if isempty(parms)
%     parms.lambda=10^(-1.5);
%     parms.mu=ones(size(hyper_coeffs));
% end
% 
% if length(parms.mu) == 1, parms.mu=parms.mu*ones(size(hyper_coeffs)); end
% 
% if size(degs,2) > size(degs,1), degs=degs'; end
% 
% switch hypermode
%     case 'filtered'
%         h=@(x) (x >= 0 & x <= 1/2) + ((sin(pi*x)).^2).*(x > 0.5 & x < 1);
%         L=max(degs);
%         filt_coeffs=h(degs/L).*hyper_coeffs;
%     case 'lasso'
%         lambda=parms.lambda; % regularization parameter
%         mu=parms.mu; % penalties parameters
%         S=@(a,k) max(0,a-k)+min(0,a+k);
%         L=length(degs);
%         for l=1:L
%             al=hyper_coeffs(l);
%             kl=lambda*mu(l);
%             filt_coeffs(l,1)=S(al,kl);
%         end
%     otherwise
%         filt_coeffs=hyper_coeffs;
% end
