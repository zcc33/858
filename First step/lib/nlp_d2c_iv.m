% function [d2c_x,dc_x] = nlp_d2c_iv(x, lam, Prob, varargin)
%
% TOMLAB gateway routine for computation of the 2nd part of the Lagrangian
% function, d2c(x) for problems with strict integer variables (not continuously 
% differentiable)
%
% nlp_d2c_iv computes
%
%           lam' * d2c(x)
%
% in
%
%   L(x,lam) =   f(x) - lam' * c(x)
% d2L(x,lam) = d2f(x) - lam' * d2c(x)
%
%
% nlp_d2c_iv calls dc = Prob.MIP.FUNCS.dc
%        as d2c_x = feval(dc,x,Prob,varargin)
% unless empty, in which case Prob.FUNCS.dc is called as
%        feval(Prob.FUNCS.dc,x,Prob,varargin);
% If both Prob.MIP.FUNCS.dc and Prob.FUNCS.dc are empty, nlp_d2c_iv calls 
% nlp_dc_iv as
%           d2c_x = nlp_dc_iv(x,Prob,varargin)
% If the problem has continuously differentiable variables nlp_d2c_iv also calls 
% the routine nlp_d2c as
%        as d2c_x = nlp_d2c(x,lam,Prob,varargin)
%
% Prob.FUNCS.d2c returns lam' * d2c(x), i.e. if there are m constraints,
% d2c is the weighted sum of m matrices of size n by n. Each weight is
% the Lagrange parameter lam(i). This means that
%
% d2c = sum_i=1:m  lam(i) * d2c(i) / dx^2
%
% d2c(i) / dx^2 is the Hessian matrix (second order matrix) for constraint i
%
% The global counter variable n_d2c is incremented

% Bjorn Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2014 by Tomlab Optimization Inc., $Release: 8.1.0$
% Written Nov 19, 2014.   Last modified December 2, 2014.

function [d2c_x,dc_x] = nlp_d2c_iv(x, lam, Prob, varargin)

global NLP_xdc NLP_dc
global n_d2c NARG
x=x(:);

m = Prob.mNonLin;
n = min(length(x),Prob.N);

tempnarg = NARG;
narg5 = Prob.MIP.FUNCS.narg(5);
dc = Prob.MIP.FUNCS.dc;

if ~isempty(NLP_xdc) && all(NLP_xdc == x)
   dc_x = NLP_dc;
else
   if ~isempty(dc)
      if narg5 > 2
         dc_x = feval(dc,x,Prob,varargin);
      else
         dc_x = feval(dc,x,Prob);
      end
   elseif ~isempty(Prob.FUNCS.dc)
      dc_x = feval(Prob.FUNCS.dc,x,Prob,varargin);
   else
      dc_x = nlp_dc_iv(x,Prob,varargin);
   end
end

% Use the information prepared by iniIntVars to call the function with integer 
% values for the strictly integer variables.
SIV_c = Prob.MIP.SIV_c;
anySIV_c = any(SIV_c);
if any(anySIV_c)
   x(anySIV_c) = round(x(anySIV_c));
end
if any(~anySIV_c)
   % Restore original function before calling nlp_d2c
   Prob.FUNCS.c = Prob.MIP.FUNCS.c;
   Prob.FUNCS.dc = dc;
   Prob.FUNCS.d2c = Prob.MIP.FUNCS.d2c;
   NARG(6) = Prob.MIP.FUNCS.narg(6);
   d2c_x = nlp_d2c(x,lam,Prob,varargin);
else
   n_d2c = n_d2c + 1;
   NARG(6) = Prob.MIP.FUNCS.narg(6);
end

if isempty(dc)
   dc = 'nlp_dc_iv';
   narg5 = 3;
end

% The implementation below is based on the procedure described in 'A Comparative 
% Study of SQP-Type Algorithms for Nonlinear and Nonconvex Mixed-Integer 
% Optimization' written by Oliver Exler, Thomas Lehmann and Klaus Schittkowski.
% The publication is available at 
% http://www.klaus-schittkowski.de/minlp_comp_study.pdf
%
% Additional modifications to the original procedure are:
% (1) Only the constraints second order derivative is estimated.
% (2) No best neighbourhood point is returned.
% (3) The check on bounds requires the variable to be at least two integer steps
%     from the bound. This must be done to account for the additional step taken
%     in the gradient estimation procedure
z_p = x;
z_m = x;
x_L = Prob.x_L;
x_U = Prob.x_U;

ix = find(anySIV_c);
for i = ix(:)'
    c_i = find(SIV_c(:,i)); % constraints in which x(i) is nondifferentiable
    % neighborhood point z^{+1}
    z_p(i) = x(i) + 1;
    % neighborhood point z^{-1}
    z_m(i) = x(i) - 1;
    d2c_x_i = zeros(m,n);
    if x_L(i) < z_m(i) && z_p(i) < x_U(i)
       if narg5 > 2
          dc_p = feval(dc,z_p,Prob,varargin);
          dc_m = feval(dc,z_m,Prob,varargin);
       else
          dc_p = feval(dc,z_p,Prob);
          dc_m = feval(dc,z_m,Prob);
       end
       d2c_x_i(c_i,i) = (dc_p(c_i,i) - dc_m(c_i,i))/2;
    elseif z_m(i) == x_L(i) && z_p(i) <= x_U(i)
       if narg5 > 2
          dc_p = feval(dc,z_p,Prob,varargin);
       else
          dc_p = feval(dc,z_p,Prob);
       end
       d2c_x_i(c_i,i) = (dc_p(c_i,i) - dc_x(c_i,i));
    elseif z_p(i) == x_U(i) && x_L(i) <= z_m(i)
       if narg5 > 2
          dc_m = feval(dc,z_m,Prob,varargin);
       else
          dc_m = feval(dc,z_m,Prob);
       end
       d2c_x_i(c_i,i) = (dc_x(c_i,i) - dc_m(c_i,i))/2;
    else
       warning(['nlp_G_iv','Cannot estimate 2nd part of Hessian for strict' ...
         ' integer variable  x(' num2str(i) '). Variable might be binary']);
    end
    
    % Multiply by lam and add to dc_x
    d2c_x(:,i) = d2c_x(:,i) + d2c_x_i'*lam(:);
    
    % Make symmetric. Symmetry cannot be assumed since some variables might not
    % be continuously differentiable, but done regardless in good faith.
    d2c_x(i,:) = d2c_x(:,i);
    
    % restore z_p and z_m to x
    z_p(i) = x(i);
    z_m(i) = x(i);
end

NLP_xdc = x;
NLP_dc = dc_x;
NARG = tempnarg;


% MODIFICATION LOG
%
% 141119   bjo   Created based on nlp_H_iv
% 141126   bjo   First version completed.
% 141202   bjo   Added dc_x as second output argument
% 150128   bjo   Added missing update of NLP_xdc and NLP_dc