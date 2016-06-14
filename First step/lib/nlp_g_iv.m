% function [g_x,f_x,f_bn,x_bn] = nlp_g_iv(x, Prob, varargin)
%
% TOMLAB gateway routine for the computation of the gradient vector g(x) for 
% problems with strict integer variables (not continuously differentiable)
%
% Additional outputs are x_bn, the point in the direct bounded neighborhood of x
% with the lowest objective function f_bn(x_bn).
%
% nlp_g_iv calls f = Prob.MIP.FUNCS.f as
%           f_x = feval(f,x,Prob,varargin)
% unless empty, in which case Prob.FUNCS.f is called as
%           f_x = feval(Prob.FUNCS.f,x,Prob,varargin);
% If both Prob.MIP.FUNCS.f and Prob.FUNCS.f are empty, nlp_g_iv calls nlp_f_iv
% as
%           f_x = nlp_f_iv(x,Prob,varargin)
% If the problem has continuously differentiable variables nlp_g_iv also calls 
% the routine nlp_g as
%        as g_x = nlp_g(x,Prob,varargin)
%
% The global counter variable n_g is incremented

% Bjorn Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2014 by Tomlab Optimization Inc., $Release: 8.1.0$
% Written Nov 19, 2014.   Last modified December 2, 2014.

function [g_x,f_x,f_bn,x_bn] = nlp_g_iv(x, Prob, varargin)

% % Communication nlp_f/g and nlp_c/dc
global NLP_x NLP_f
% Communication nlp_g/H
global NLP_g NLP_xg
global n_g NARG

x = x(:);
if ~isempty(NLP_xg) && all(x == NLP_xg)
   if ~isempty(NLP_f) && ~isempty(NLP_x) && all(x == NLP_x)
      f_x = NLP_f;
   else
      f_x = nlp_f_iv(x,Prob,varargin);
   end
   g_x = NLP_g;
   f_bn = []; x_bn = [];
   return
end

n = length(x);

tempnarg = NARG;
narg1 = Prob.MIP.FUNCS.narg(1);

f = Prob.MIP.FUNCS.f;
if ~isempty(NLP_x) && all(NLP_x == x)
   f_x = NLP_f;
else
   if ~isempty(f)
      if narg1 > 2
         f_x = feval(f,x,Prob,varargin);
      else
         f_x = feval(f,x,Prob);
      end
   elseif ~isempty(Prob.FUNCS.f)
      f_x = feval(Prob.FUNCS.f,x,Prob,varargin);
   else
      f_x = nlp_f_iv(x,Prob,varargin);
   end
end

% Use the information prepared by iniIntVars to call the function with integer 
% values for the strictly integer variables.
SIV_f = Prob.MIP.SIV_f;
if any(SIV_f)
   x(SIV_f) = round(x(SIV_f));
end
if any(~SIV_f)
   % Restore original functions before calling nlp_g
   Prob.FUNCS.f = f;
   Prob.FUNCS.g = Prob.MIP.FUNCS.g;
   NARG(2) = Prob.MIP.FUNCS.narg(2);
   Prob.g_k = zeros(n,1);
   Prob.g_k(~SIV_f) = NaN;
   g_x = nlp_g(x, Prob, varargin);
else
   g_x = nan(n,1);
   n_g = n_g + 1;
end
if ~any(SIV_f)
   f_bn = []; x_bn = [];
   NARG = tempnarg;
   return
end

% The implementation below is based on the procedure described in 'A Comparative 
% Study of SQP-Type Algorithms for Nonlinear and Nonconvex Mixed-Integer 
% Optimization' written by Oliver Exler, Thomas Lehmann and Klaus Schittkowski.
% The publication is available at 
% http://www.klaus-schittkowski.de/minlp_comp_study.pdf
%
% Additional modifications to the original procedure are:
% (1) Only the derivative of the objective function is estimated.
% (2) The best neighborhood point selection is based only on the objective.
z_p = x;
z_m = x;
x_L = Prob.x_L;
x_U = Prob.x_U;
f_bn = Inf;
x_bn = x;

ix = find(SIV_f);
for i = ix(:)'
    % neighborhood point z^{+1}
    z_p(i) = x(i) + 1;
    % neighborhood point z^{-1}
    z_m(i) = x(i) - 1;
    if x_L(i) < x(i) && x(i) < x_U(i)
       if narg1 > 2
          f_p = feval(f,z_p,Prob,varargin);
          f_m = feval(f,z_m,Prob,varargin);
       else
          f_p = feval(f,z_p,Prob);
          f_m = feval(f,z_m,Prob);
       end
       if f_p < f_bn
          f_bn = f_p;
          x_bn = z_p;
       end
       if f_m < f_bn
          f_bn = f_m;
          x_bn = z_m;
       end
       g_x(i) = (f_p-f_m)/2;
    elseif x(i) == x_L(i)
       if narg1 > 2
          f_p = feval(f,z_p,Prob,varargin);
       else
          f_p = feval(f,z_p,Prob);
       end
       if f_p < f_bn
          f_bn = f_p;
          x_bn = z_p;
       end
       g_x(i) = (f_p-f_x);
    elseif x(i) == x_U(i)
       if narg1 > 2
          f_m = feval(f,z_p,Prob,varargin);
       else
          f_m = feval(f,z_p,Prob);
       end
       if f_m < f_bn
          f_bn = f_m;
          x_bn = z_m;
       end
       g_x(i) = (f_x-f_m);
    end
    % restore z_p and z_m to x
    z_p(i) = x(i);
    z_m(i) = x(i);
end

NLP_xg = x;
NLP_g  = g_x;
NLP_x = x;
NLP_f = f_x;

NARG = tempnarg;

% MODIFICATION LOG
%
% 141119   bjo   Created based on nlp_g
% 141126   bjo   First version completed.
% 141202   bjo   Added f_x as second output argument
% 150128   bjo   Add missing update of NLP_f