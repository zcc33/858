% function Prob = bincont2lin(Prob, idx_prod, idx_bin, idx_cont)
%
% Adds constraints when modeling with binary variables which are multiplied
% by integer or continuous variables. This is the most efficient way to
% get rid off quadratic objectives or constraints.
%
% prod = bin * cont. The problem should be built with the extra variables
% prod in place of the bin*cont products. The indices of the unique product
% variables are needed to convert the problem properly.
%
% Four inequalities are added to the problem:
%
% prod <= cont - min(0,x_L) * (1 - bin)
% prod >= cont - max(0,x_U) * (1 - bin)
% prod <= x_U * bin
% prod >= x_L * bin
%
% By these constraints, prod will always equal bin * cont.
%
% INPUT PARAMETERS
%
% Prob         Problem structure to be converted
% idx_prod     Indices for product variables
% idx_bin      Indices for binary variables
% idx_cont     Indices for continuous/integer variables
%
% OUTPUT PARAMETERS
%
% Prob         Problem structure with added constraints

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2015 by Tomlab Optimization Inc., $Release: 6.1.0$
% Written Oct 4, 2005.   Last modified Feb 19, 2015.

function Prob = bincont2lin(Prob, idx_prod, idx_bin, idx_cont)

if nargin < 4
    error('bincont2lin requires 4 inputs');
end

len = length(idx_prod);

if len ~= length(idx_bin)
    error('idx_prod and idx_bin do not have the same length');
end

if len ~= length(idx_cont)
    error('idx_prod and idx_cont do not have the same length');
end

if ~all(isfinite(Prob.x_U(idx_cont)))
    error('An upper bound is required on the integer/continuous variables');
end

if ~all(isfinite(Prob.x_L(idx_cont)))
    error('A lower bound is required on the integer/continuous variables');
end

idx_prod = full(idx_prod(:));
idx_bin  = full(idx_bin(:));
idx_cont = full(idx_cont(:));

% Tighten bounds on prod
Prob.x_L(idx_prod) = max(Prob.x_L(idx_prod),min(0,Prob.x_L(idx_cont)));
Prob.x_U(idx_prod) = min(Prob.x_U(idx_prod),max(0,Prob.x_U(idx_cont)));

x_U = Prob.x_U(idx_cont);
x_L = Prob.x_L(idx_cont);

% prod <= cont - min(0,x_L) * (1 - bin)
A1 = sparse(repmat((1:len)',3,1),[idx_cont;idx_bin;idx_prod],[ones(len,1);min(0,x_L);-ones(len,1)],len,Prob.N);
bL1 = min(0,x_L);
bU1 = inf(len,1);

% prod >= cont - max(0,x_U) * (1 - bin)
A2 = sparse(repmat((1:len)',3,1),[idx_cont;idx_bin;idx_prod],[ones(len,1);max(0,x_U);-ones(len,1)],len,Prob.N);
bL2 = -inf(len,1);
bU2 = max(0,x_U);

% prod <= x_U * bin
A3 = sparse(repmat((1:len)',2,1),[idx_prod;idx_bin],[ones(len,1);-x_U],len,Prob.N);
bL3 = zeros(len,1);
bU3 = inf(len,1);

% prod >= x_L * bin
A4 = sparse(repmat((1:len)',2,1),[idx_prod;idx_bin],[ones(len,1);-x_L],len,Prob.N);
bL4 = -inf(len,1);
bU4 = zeros(len,1);

Prob.A = [Prob.A;A1;A2;A3;A4];
Prob.b_U = [Prob.b_U; bU1; bU2; bU3; bU4];
Prob.b_L = [Prob.b_L; bL1; bL2; bL3; bL4];

Prob.mLin = size(Prob.A,1);

% MODIFICATION LOG
%
% 051004 med  Created.
% 070220 med  Updated b_U vector
% 080522 med  Lower bounds now supported as well
% 150219 rut  Correct behavior when "cont" is negative
