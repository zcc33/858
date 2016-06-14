% function [dc_x,c_x,c_bn,x_bn] = nlp_dc_iv(x, Prob, varargin)
%
% TOMLAB gateway routine for computation of the constraint Jacobian matrix dc(x)
% for problems with strict integer variables (not continuously differentiable)
%
% nlp_dc_iv calls c = Prob.MIP.FUNCS.c as
%           c_x = feval(c,x,Prob,varargin)
% unless empty, in which case Prob.FUNCS.c is called as
%           c_x = feval(Prob.FUNCS.c,x,Prob,varargin);
% If both Prob.MIP.FUNCS.c and Prob.FUNCS.c are empty, nlp_dc_iv calls nlp_c_iv
% as
%           dc_x = nlp_c_iv(x,Prob,varargin)
% If the problem has continuously differentiable variables nlp_dc_iv also calls 
% the routine nlp_dc as
%        as dc_x = nlp_dc(x,Prob,varargin)
%
% dc is a mL x n matrix, where mL = number of nonlinear constraints
%                               n = number of variables
%
% Additional outputs are x_bn, the point in the direct bounded neighborhood of x
% with the lowest Inf-norm of the distances between the constraint functions 
% c_bn(x_bn) and their noninfinite bounds.
%
% The global counter variable n_dc is incremented
%
% If Prob.ConsDiff > 0, dc is estimated numerically
% If Prob.CheckNaN > 0, NaN elements in dc are estimated numerically
%
% Numerical method is dependent on Prob.ConsDiff and Prob.CheckNaN

% Bjorn Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2014 by Tomlab Optimization Inc., $Release: 8.1.0$
% Written Nov 19, 2014.   Last modified December 2, 2014.

function [dc_x,c_x,c_bn,x_bn] = nlp_dc_iv(x, Prob, varargin)

global NLP_xc NLP_xdc NLP_c NLP_dc  % Communication nlp_c/dc
global n_dc NARG

x=x(:);
if ~isempty(NLP_xdc) && all(x == NLP_xdc)
   if ~isempty(NLP_xc) && all(x == NLP_xc)
      c_x = NLP_c;
   else
      c_x = nlp_c_iv(x,Prob,varargin);
   end
   dc_x = NLP_dc;
   return
end

m = Prob.mNonLin;
n = min(length(x),Prob.N);

narg4 = Prob.MIP.FUNCS.narg(4);
tempnarg = NARG;

% Use the information prepared by iniIntVars to call the function with integer 
% values for the strictly integer variables.
SIV_c = Prob.MIP.SIV_c;
anySIV_c = any(SIV_c);
if any(anySIV_c)
   x(anySIV_c) = round(x(anySIV_c));
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
    
if any(~anySIV_c)
   % Restore original function before calling nlp_dc
   NARG(5) = Prob.MIP.FUNCS.narg(5);
   NLP_dc(:,~anySIV_c) = NaN;
   Prob.FUNCS.c = c;
   Prob.FUNCS.dc = Prob.MIP.FUNCS.dc;
   Prob.dc_k = zeros(m,n);
   Prob.dc_k(~SIV_c) = NaN;
   dc_x = nlp_dc(x, Prob, varargin);
else
   dc_x = nan(m,n);
   n_dc = n_dc + 1;
   NARG(5) = Prob.MIP.FUNCS.narg(5);
end
if ~any(SIV_c)
   NARG = tempnarg;
   c_bn = []; x_bn = [];
   return
end

% The implementation below is based on the procedure described in 'A Comparative 
% Study of SQP-Type Algorithms for Nonlinear and Nonconvex Mixed-Integer 
% Optimization' written by Oliver Exler, Thomas Lehmann and Klaus Schittkowski.
% The publication is available at 
% http://www.klaus-schittkowski.de/minlp_comp_study.pdf
%
% Additional modifications to the original procedure are:
% (1) Only the constraint derivative is estimated.
% (2) The best neighborhood point selection is based on the minimum value of 
%     the Inf-norm of the vector of distances between the constraint functions
%     and their noninfinite lower and upper bounds
z_p = x;
z_m = x;
x_L = Prob.x_L;
x_U = Prob.x_U;
c_bn = Inf;
c_bn_in = Inf*ones(m,1);
x_bn = x;
c_L = Prob.c_L;
c_U = Prob.c_U;
nIL = ~isinf(c_L);
nIU = ~isinf(c_U);
cLnIL = c_L(nIL);
cUnIU = c_U(nIU);

ix = find(any(SIV_c));
for i = ix(:)'
    c_i = find(SIV_c(:,i));
    % neighborhood point z^{+1}
    z_p(i) = x(i) + 1;
    % neighborhood point z^{-1}
    z_m(i) = x(i) - 1;
    if x_L(i) < x(i) && x(i) < x_U(i)
       if narg4 > 2
          c_p = feval(c,z_p,Prob,varargin);
          c_m = feval(c,z_m,Prob,varargin);
       else
          c_p = feval(c,z_p,Prob);
          c_m = feval(c,z_m,Prob);
       end
       c_p_in = min(norm(c_p(nIL)-cLnIL,Inf),norm(cUnIU-c_p(nIU),Inf));
       if c_p_in < c_bn_in
          c_bn = c_p;
          c_bn_in = c_p_in;
          x_bn = z_p;
       end
       c_m_in = min(norm(c_m(nIL)-cLnIL,Inf),norm(cLnIL-c_m(nIL),Inf));
       if c_m_in < c_bn_in
          c_bn = c_m;
          c_bn_in = c_m_in;
          x_bn = z_m;
       end
       dc_x(c_i,i) = (c_p(c_i)-c_m(c_i))/2;
    elseif x(i) == x_L(i)
       if narg4 > 2
          c_p = feval(c,z_p,Prob,varargin);
       else
          c_p = feval(c,z_p,Prob);
       end
       c_p_in = min(norm(c_p(nIL)-cLnIL,Inf),norm(cUnIU-c_p(nIU),Inf));
       if c_p_in < c_bn_in
          c_bn = c_p;
          c_bn_in = c_p_in;
          x_bn = z_p;
       end
       dc_x(c_i,i) = (c_p(c_i)-c_x(c_i));
    elseif x(i) == x_U(i)
       if narg4 > 2
          c_m = feval(c,z_m,Prob,varargin);
       else
          c_m = feval(c,z_m,Prob);
       end
       c_m_in = min(norm(c_m(nIL)-cLnIL,Inf),norm(cLnIL-c_m(nIL),Inf));
       if c_m_in < c_bn_in
          c_bn = c_m;
          c_bn_in = c_m_in;
          x_bn = z_m;
       end
       dc_x(c_i,i) = (c_x(c_i)-c_m(c_i));
    end
    % restore z_p and z_m to x
    z_p(i) = x(i);
    z_m(i) = x(i);
end
NLP_xdc = x;
NLP_dc  = dc_x;
NLP_xc = x;
NLP_c = c_x;
NARG = tempnarg;

% MODIFICATION LOG
%
% 141119   bjo   Created based on nlp_dc
% 141126   bjo   First version completed.
% 141202   bjo   Added c_x as second output argument
% 150128   bjo   Add missing update of NLP_xc and NLP_c