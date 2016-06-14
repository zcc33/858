% nlp2_fgH.m
%
% function [f, g, H] = nlp2_fgH(x, Prob, varargin)
%
% nlp2_fgH returns the function f(x), the gradient g(x), and the
% Hessian H(x) for a nonlinear problem.
%
% nlp2_fgH calls the user function that returns the function, and possibly
% the user function that returns the gradient vector and the Hessian matrix.
%
% H(x) may be sparse
%
% if nargout > 1, nlp_fgH also computes the gradient g at x, g(x).
%
% if nargout > 2, nlp_fgH also computes the Hessian H at x, H(x).
%
% nlp2_fgH is used when implementing the OPTIM TB 2.0 compatibility interface

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2014 by Tomlab Optimization Inc., Sweden. $Release: 8.0.0$
% Written July 6, 1999.    Last modified March 2, 2014.

function [f, g, H] = nlp2_fgH(x, Prob, varargin)

nargin;
global n_f n_g n_H NARG

fgH=Prob.FUNCSX.f;

n_f = n_f + 1;

if NARG(1) > 2
   f = feval(fgH, x, Prob, varargin{:});
elseif NARG(1) == 1
   f = feval(fgH, x);
else
   f = feval(fgH, x, Prob);
end

if nargout > 1
   gFunc=Prob.FUNCSX.g;
   if ~isempty(gFunc)
      n_g = n_g + 1;
      if NARG(2) > 2
         g = feval(gFunc, x, Prob, varargin{:});
      elseif NARG(2) == 1
         g = feval(gFunc, x);
      else
         g = feval(gFunc, x, Prob);
      end
   else
      g=[];
   end
   if nargout > 2
      HFunc=Prob.FUNCSX.H;
      if ~isempty(HFunc)
         n_H = n_H + 1;
         if NARG(3) > 2
            H = feval(HFunc, x, Prob, varargin{:});
         elseif NARG(3) == 1
            H = feval(HFunc, x);
         else
            H = feval(HFunc, x, Prob);
         end
      else
         H=[];
      end
   end
end

% MODIFICATION LOG:
%
% 011204 hkh  H=[] was misplaced
% 030111 hkh  Revision for v4.0
% 060814 med  FUNCS used for callbacks instead
% 140302 hkh   Assume Prob used if NARG value == 0
