% function [H_x,g_x] = nlp_H_iv(x, Prob, varargin)
%
% TOMLAB gateway routine for the computation of the Hessian matrix H(x) for 
% problems with strict integer variables (not continuously differentiable)
%
% H_x(i,j) is d2f/dx_i,dx_j
% Currently the symmetry is not taken into account
%
% nlp_H_iv calls Prob.MIP.FUNCS.g as
%           g_x = feval(g,x,Prob,varargin)
% unless empty, in which case Prob.FUNCS.g is called as
%           g_x = feval(Prob.FUNCS.g,x,Prob,varargin);
% If both Prob.MIP.FUNCS.g and Prob.FUNCS.g are empty, nlp_H_iv calls nlp_g_iv
% as
%           g_x = nlp_g_iv(x,Prob,varargin)
% If the problem has continuously differentiable variables nlp_H_iv also calls 
% the routine nlp_H as
%           H_x = nlp_H(x,Prob,varargin)
%
% The global counter variable n_H is incremented

% Bjorn Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2014 by Tomlab Optimization Inc., $Release: 8.1.0$
% Written Nov 19, 2014.   Last modified December 2, 2014.

function [H_x,g_x] = nlp_H_iv(x, Prob, varargin)

% Communication nlp_g/H
global NLP_g NLP_xg
global n_H NARG

narg2 = Prob.MIP.FUNCS.narg(2);
tempnarg = NARG;

g = Prob.MIP.FUNCS.g;
n = length(x);
if ~isempty(NLP_xg) && all(NLP_xg == x)
   g_x = NLP_g;
else
   if ~isempty(g)
      if narg2 > 2
         g_x = feval(g,x,Prob,varargin);
      else
         g_x = feval(g,x,Prob);
      end
   elseif ~isempty(Prob.FUNCS.g)
      g_x = feval(Prob.FUNCS.g,x,Prob,varargin);
   else
      g_x = nlp_g_iv(x,Prob,varargin);
   end
end

% Use the information prepared by iniIntVars
SIV_f = Prob.MIP.SIV_f;
if any(SIV_f)
   x(SIV_f) = round(x(SIV_f));
end

if any(~SIV_f)
   % Restore original functions before calling nlp_H
   Prob.FUNCS.f = Prob.MIP.FUNCS.f;
   Prob.FUNCS.g = g;
   Prob.FUNCS.H = Prob.MIP.FUNCS.H;
   % TODO: Find out how to avoid evaluation on SIV_ix
   % Current issue:
   % nlp_H->FD-Hess->fdng->nlp_f evaluates with noninteger values, despite setting
   % NLP_g(~SIV_f) = NaN; % Simplest way
   Prob.g_k = zeros(Prob.N,1);
   Prob.g_k(~SIV_f) = NaN;
   NARG(3) = Prob.MIP.FUNCS.narg(3);
   H_x = nlp_H(x,Prob,varargin);
else
   H_x = nan(n,n);
   n_H = n_H + 1;
   NARG = Prob.MIP.FUNCS.narg;
end
if ~any(SIV_f)
   NARG(3) = tempnarg;
   return
end
if isempty(g)
   g = 'nlp_g_iv';
   narg2 = 3;
end

% The implementation below is based on the procedure described in 'A Comparative 
% Study of SQP-Type Algorithms for Nonlinear and Nonconvex Mixed-Integer 
% Optimization' written by Oliver Exler, Thomas Lehmann and Klaus Schittkowski.
% The publication is available at 
% http://www.klaus-schittkowski.de/minlp_comp_study.pdf
%
% Additional modifications to the original procedure are:
% (1) The constraint derivatives are not estimated.
% (2) The best neighborhood point selection is not performed
% (3) The check on bounds requires the variable to be at least two integer steps
%     from the bound. This must be done to account for the additional step taken
%     in the gradient estimation procedure
z_p = x;
z_m = x;
x_L = Prob.x_L;
x_U = Prob.x_U;

ix = find(SIV_f);
for i = ix(:)'
    % neighborhood point z^{+1}
    z_p(i) = x(i) + 1;
    % neighborhood point z^{-1}
    z_m(i) = x(i) - 1;
    if x_L(i) < z_m(i) && z_p(i) < x_U(i)
       if narg2 > 2
          g_p = feval(g,z_p,Prob,varargin);
          g_m = feval(g,z_m,Prob,varargin);
       else
          g_p = feval(g,z_p,Prob);
          g_m = feval(g,z_m,Prob);
       end
       H_x(:,i) = (g_p - g_m)/2;
       H_x(i,:) = H_x(:,i);
    elseif z_m(i) == x_L(i) && z_p(i) <= x_U(i)
       if narg2 > 2
          g_p = feval(g,z_p,Prob,varargin);
       else
          g_p = feval(g,z_p,Prob);
       end
       H_x(:,i) = (g_p - g_x);
       H_x(i,:) = H_x(:,i);
    elseif z_p(i) == x_U(i) && x_L(i) <= z_m(i)
       if narg2 > 2
          g_m = feval(g,z_m,Prob,varargin);
       else
          g_m = feval(g,z_m,Prob);
       end
       H_x(:,i) = (g_x - g_m)/2;
       H_x(i,:) = H_x(:,i);
    else
       warning(['nlp_G_iv','Cannot estimate Hessian for strict integer' ... 
                'variable  x(' num2str(i) '). Variable might be binary']);
    end
    % restore z_p and z_m to x
    z_p(i) = x(i);
    z_m(i) = x(i);
end

NLP_xg = x;
NLP_g = g_x;
NARG = tempnarg;

% MODIFICATION LOG
%
% 141119   bjo   Created based on nlp_H
% 141126   bjo   First version completed.
% 141202   bjo   Added g_x as second output argument
% 140128   bjo   Added missing update of NLP_xg and NLP_g