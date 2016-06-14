% Defines objective function for mixed-integer nonlinear programming 
% (MINLP) problems from 'A Collection of 186 Test Problems for Nonlinear 
% Mixed-Integer Programming in Fortran' by Klaus Schittkowski.
%
% function f = mi_test_probs_f(x, Prob)
function f = mi_test_probs_f(x, Prob)

x=x(:);

% Call MEX
[f, c] = mi_test_probs(1, x);

% MODIFICATION LOG
%
% 130414 bjo  Wrote file