% function [f_x,g_x,c_x,dc_x,f_bn,c_bn,x_bn] = nlp_fgcdc_iv(x, Prob, varargin)
%
% TOMLAB gateway routine for the computation of the gradient vector g(x) for 
% problems with strict integer variables (not continuously differentiable)
%
% Additional outputs are x_bn, the point in the direct bounded neighborhood of x
% with the lowest objective function f_bn(x_bn).
%
% nlp_fgcdc_iv calls f = Prob.MIP.FUNCS.f as
%           f_x = feval(f,x,Prob,varargin)
% unless empty, in which case Prob.FUNCS.f is called as
%           f_x = feval(Prob.FUNCS.f,x,Prob,varargin);
% If both Prob.MIP.FUNCS.f and Prob.FUNCS.f are empty, nlp_fgcdc_iv calls
% nlp_f_iv as
%           f_x = nlp_f_iv(x,Prob,varargin)
% If the problem has continuously differentiable variables nlp_fgcdc_iv also calls 
% the routine nlp_g as
%        as g_x = nlp_g(x,Prob,varargin)
%
% The global counter variable n_g is incremented

% Bjorn Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2014 by Tomlab Optimization Inc., $Release: 8.1.0$
% Written Dec 1, 2014.   Last modified December 1, 2014.

function [f_x,g_x,c_x,dc_x,f_bn,c_bn,x_bn] = nlp_fgcdc_iv(x, Prob, varargin)

% Communication nlp_f/g
global NLP_x NLP_f
% Communication nlp_g/H
global NLP_g NLP_xg
% Communication nlp_c/dc
global NLP_xc NLP_xdc NLP_c NLP_dc
% Counter variables
global n_g n_dc
% Vector with number of arguments
global NARG

x = x(:);
f_x = []; g_x = []; c_x = []; dc_x = []; f_bn = []; c_bn = []; x_bn = [];

n = min(length(x),Prob.N);
m = Prob.mNonLin;

tempnarg = NARG;
narg1 = Prob.MIP.FUNCS.narg(1);
narg4 = Prob.MIP.FUNCS.narg(4);

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

c = Prob.MIP.FUNCS.c;
if ~isempty(NLP_xc) && all(NLP_xc == x)
   c_x = NLP_c;
else
   if ~isempty(c)
      if narg4 > 2
         c_x = feval(c,x,Prob,varargin);
      else
         c_x = feval(c,x,Prob);
      end
   elseif ~isempty(Prob.FUNCS.c)
      c_x = feval(Prob.FUNCS.c,x,Prob,varargin);
   else
      c_x = nlp_c_iv(x,Prob,varargin);
   end
end
if isempty(NLP_dc)
   NLP_dc = zeros(m,n);
end
if ~isempty(NARG)
   NARG(1) = Prob.MIP.FUNCS.narg(1);
   NARG(4) = Prob.MIP.FUNCS.narg(4);
end

if ~isempty(NLP_xg) && all(x == NLP_xg)
   g_x = NLP_g;
end

if ~isempty(NLP_xdc) && all(x == NLP_xdc)
   dc_x = NLP_dc;
end

if ~isempty(g_x) && ~isempty(dc_x)
   % No estimation required
   return;
end

ConIntVars = Prob.MIP.ConIntVars;
SIV = ~ConIntVars;

% Use the information prepared by iniIntVars to call the objective and 
% constraint functions with integer values for the strictly integer variables.
if any(SIV)
   x(SIV) = round(x(SIV));
end
if any(ConIntVars)
   % Restore original functions before calling nlp_g and nlp_dc
   Prob.FUNCS.f = f;
   Prob.FUNCS.g = Prob.MIP.FUNCS.g;
   Prob.FUNCS.c = c;
   Prob.FUNCS.dc = Prob.MIP.FUNCS.dc;
   NLP_dc(:,ConIntVars) = NaN;
   NARG(2) = Prob.MIP.FUNCS.narg(2);
   NARG(5) = Prob.MIP.FUNCS.narg(5);
   Prob.g_k = zeros(n,1);
   Prob.dc_k = zeros(m,n);
   Prob.g_k(ConIntVars) = NaN;
   Prob.dc_k(ConIntVars) = NaN;
   g_x = nlp_g(x, Prob, varargin);
   dc_x = nlp_dc(x, Prob, varargin);
else
   g_x = nan(n,1);
   dc_x = nan(m,n);
   n_g = n_g + 1;
   n_dc = n_dc + 1;
   NARG(2) = Prob.MIP.FUNCS.narg(2);
   NARG(5) = Prob.MIP.FUNCS.narg(5);
end
if ~any(SIV)
   f_bn = []; x_bn = [];
   c_bn = []; x_bn = [];
   NARG = tempnarg;
   return
end

% The implementation below is based on the procedure described in 'A Comparative 
% Study of SQP-Type Algorithms for Nonlinear and Nonconvex Mixed-Integer 
% Optimization' written by Oliver Exler, Thomas Lehmann and Klaus Schittkowski.
% The publication is available at 
% http://www.klaus-schittkowski.de/minlp_comp_study.pdf
z_p = x;
z_m = x;
x_L = Prob.x_L;
x_U = Prob.x_U;
f_bn = Inf;
x_bn = x;
c_bn = Inf;
x_bn = x;
c_L = Prob.c_L;
c_U = Prob.c_U;

% TODO: Efficiency improvements: Set these initially in the Prob-struct
cEQ = c_L == c_U;
cIE = c_L ~= c_U;
% Safeguard
if ~any(cEQ)
   cEQ = [];
end
if ~any(cIE)
   cIE = [];
end

eps = Prob.optParam.DiffInt;

ix = find(SIV);
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
       if narg4 > 2
          c_p = feval(c,z_p,Prob,varargin);
          c_m = feval(c,z_m,Prob,varargin);
       else
          c_p = feval(c,z_p,Prob);
          c_m = feval(c,z_m,Prob);
       end
       c_in_inf = c_p(cEQ);
       c_eq_inf = min(c_p(cIE) - c_L(cIE),min(0,c_U(cIE) - c_p(cIE)));
       if f_p < f_bn && ...
          norm([c_in_inf;c_eq_inf],Inf) <= eps
          f_bn = f_p;
          c_bn = c_p;
          x_bn = z_p;
       end
       c_in_inf = c_m(cEQ);
       c_eq_inf = min(c_m(cIE) - c_L(cIE),min(0,c_U(cIE) - c_m(cIE)));
       if f_m < f_bn && ...
          norm([c_in_inf;c_eq_inf],Inf) <= eps
          f_bn = f_m;
          c_bn = c_m;
          x_bn = z_m;
       end
       g_x(i) = (f_p - f_m)/2;
       dc_x(:,i) = (c_p - c_m)/2;
    elseif x(i) == x_L(i)
       if narg1 > 2
          f_p = feval(f,z_p,Prob,varargin);
       else
          f_p = feval(f,z_p,Prob);
       end
       if narg4 > 2
          c_p = feval(c,z_p,Prob,varargin);
       else
          c_p = feval(c,z_p,Prob);
       end
       c_in_inf = c_p(cEQ);
       c_eq_inf = min(c_p(cIE) - c_L(cIE),min(0,c_U(cIE) - c_p(cIE)));
       if f_p < f_bn && ...
          norm([c_in_inf;c_eq_inf],Inf) <= eps
          f_bn = f_p;
          c_bn = c_p;
          x_bn = z_p;
       end
       g_x(i) = f_p - f_x;
       dc_x(:,i) = c_p - c_x;
    elseif x(i) == x_U(i)
       if narg1 > 2
          f_m = feval(f,z_p,Prob,varargin);
       else
          f_m = feval(f,z_p,Prob);
       end
       if narg4 > 2
          c_m = feval(c,z_m,Prob,varargin);
       else
          c_m = feval(c,z_m,Prob);
       end
       c_in_inf = c_m(cEQ);
       c_eq_inf = min(c_m(cIE) - c_L(cIE),min(0,c_U(cIE) - c_m(cIE)));
       if f_m < f_bn && ...
          norm([c_in_inf;c_eq_inf],Inf) <= eps
          f_bn = f_m;
          c_bn = c_m;
          x_bn = z_m;
       end
       g_x(i) = f_x - f_m;
       dc_x(:,i) = c_x - c_m;
    end
    % restore z_p and z_m to x
    z_p(i) = x(i);
    z_m(i) = x(i);
end

NLP_x = x;
NLP_f = f_x;
NLP_xg = x;
NLP_g  = g_x;
NLP_xc = x;
NLP_c = c_x;
NLP_xdc = x;
NLP_dc = dc_x;

NARG = tempnarg;

% MODIFICATION LOG
%
% 141201   bjo   Created based on nlp_g_iv and nlp_dc_iv
% 150128   bjo   Added missing updates of global NLP_-variables.
% 150209   bjo   Added safeguard against cases without equality or nonequality
%                constraints
%                Added missing indexing when calculating distances to constraint
%                bounds