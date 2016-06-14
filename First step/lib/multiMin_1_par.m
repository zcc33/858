% function F0 = multiMin_1_par(X0,Prob)
%
% Compute loop 1 in multiMin with parfor using Prob.Threads threads

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2001-2013 by Tomlab Optimization Inc., $Release: 7.10.0$
% Written July 13, 2006.    Last modified Sept 16, 2013.

function F0 = multiMin_1_par(X0,Prob)
M  = size(X0,2);
F0 = zeros(1,M);
parfor P = 1:M
    F0(1,P) = nlp_f(X0(:,P), Prob);         % Compute f(x_0) for each x_0
end

% MODIFICATION LOG:
%
% 110723  hkh  Written
% 130916  hkh  Utilize parfor
