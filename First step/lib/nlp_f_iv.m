% function [f,Result] = nlp_f_iv(x, Prob, varargin)
%
% TOMLAB gateway routine for the computation of function values f(x) for 
% problems with strict integer variables (not continuously differentiable)
%
% nlp_f_iv calls the routine
%       as [f,Result] = nlp_f(x,Prob,varargin)
%
% The global counter variable n_f is incremented by nlp_f

% Bjorn Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2014 by Tomlab Optimization Inc., $Release: 8.1.0$
% Written Nov 19, 2014.   Last modified November 26, 2014.

function [f,Result] = nlp_f_iv(x, Prob, varargin)

global n_f BUILDP NARG NARGO PartSep
global mad_f
% Communication nlp_f/g and nlp_c/dc
global NLP_x NLP_f NLP_pSepIndex

if ~isempty(NARG)
   fnarg = NARG(1);
   NARG(1) = Prob.MIP.FUNCS.narg(1);
end
% Restore original function before making any calls
Prob.FUNCS.f = Prob.MIP.FUNCS.f;

% Use the information prepared by iniIntVars to call the function with integer 
% values for the strictly integer variables.
SIV_f = Prob.MIP.SIV_f;
x(SIV_f) = round(x(SIV_f));
NLP_x = x;

Result = [];
if Prob.GATEF > 0
   % [f,Result]=feval(Prob.GATE.f, x, Prob, varargin{:} );
   f = feval(Prob.FUNCS.f, x(:), Prob);
   return
elseif Prob.simType > 0
   f = sim_fc(x(:), Prob, varargin{:} );
   return
end
x = x(:);

if ~isempty(NLP_x)
   if length(x)~=length(NLP_x)
      NLP_x=[];
   elseif PartSep
      if NLP_pSepIndex == Prob.PartSep.index
         if all(x == NLP_x) 
            f=NLP_f;
            return
         end
      end
   end
end

n_f  = n_f+1;
Func = Prob.FUNCS.f;
N    = min(length(x),Prob.N);

if isempty(Func)
   NLP_x=[];
   NLP_f=[];
   f=NaN;
   return
else
   if isempty(NARG)
      p = xnargin(Func);
   else
      p = NARG(1);
   end
   if isempty(NARGO)
      NARGO(1) = nargout(Func);
   end
   
   if Prob.ADObj == 1
      if NARGO(1) > 1
         if p > 2
            [mad_f,Result]=feval(Func, fmad(x(1:N),speye(N)), Prob,varargin{:});
         elseif p == 1
            [mad_f,Result]=feval(Func, fmad(x(1:N),speye(N)));
         else
            [mad_f,Result]=feval(Func, fmad(x(1:N),speye(N)), Prob);
         end
      else
         if p > 2
            mad_f=feval(Func, fmad(x(1:N),speye(N)), Prob, varargin{:});
         elseif p == 1
            mad_f=feval(Func, fmad(x(1:N),speye(N)));
         else
            mad_f=feval(Func, fmad(x(1:N),speye(N)), Prob);
         end
      end
      f=getvalue(mad_f);
   else
      BUILDP0=BUILDP;
      BUILDP=[];
      if NARGO(1) > 1
         if p > 2
            [f,Result]=feval(Func, x(1:N), Prob, varargin{:} );
         elseif p == 1
            [f,Result]=feval(Func, x(1:N));
         else
            [f,Result]=feval(Func, x(1:N), Prob);
         end
      else
         if p > 2
            f=feval(Func, x(1:N), Prob, varargin{:} );
         elseif p == 1
            f=feval(Func, x(1:N));
         else
            f=feval(Func, x(1:N), Prob);
         end
      end
      BUILDP=BUILDP0;
   end
end
if isnan(f), f=Inf; end

NLP_x=x;
NLP_f=f;

if PartSep
   % Number of partially separable functions.
   NLP_pSepIndex = Prob.PartSep.index;    % Function computed
end

if BUILDP > 0
   pbuild(x,f);
end

if ~isempty(NARG)
   NARG(1) = fnarg;
end

% MODIFICATION LOG
%
% 141119   bjo   Created based on nlp_f
% 141126   bjo   First version completed.
% 141223   bjo   Merged with all code from nlp_f to avoid recursion when calling
%                subsolvers from MIPNLP-solvers.