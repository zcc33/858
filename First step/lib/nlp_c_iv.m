% function c = nlp_c_iv(x, Prob, varargin)
%
% TOMLAB gateway routine for computation of the constraint vector c(x) for 
% problems with strict integer variables (not continuously differentiable)
%
% nlp_c_iv calls the routine nlp_c
%        as c = nlp_c(x,Prob,varargin)
%
% The global counter variable n_c is incremented by nlp_c

% Bjorn Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2014 by Tomlab Optimization Inc., $Release: 8.1.0$
% Written Nov 19, 2014.   Last modified November 26, 2014.

function c = nlp_c_iv(x, Prob, varargin)

global NARG

if ~isempty(NARG)
   cnarg = NARG(4);
   NARG(4) = Prob.MIP.FUNCS.narg(4);
end

% Call the constraints function with only integer values for the strictly 
% integer variables.
SIV_c = any(Prob.MIP.SIV_c);
x(SIV_c) = round(x(SIV_c));
% Restore original function before calling nlp_c
Prob.FUNCS.c = Prob.MIP.FUNCS.c;
c = nlp_c(x,Prob,varargin);

if ~isempty(NARG)
   NARG(4) = cnarg;
end

% MODIFICATION LOG
%
% 141119   bjo   Created based on nlp_dc
% 141126   bjo   First version completed.