% function rngInit = rngset(RandState)
%
% INPUT:
% RandState   Value RS determines how to obtain random numbers using the pseudo-random generator
%          RS >0  Call rng(RS), i.e. reinitialize random generator with value RS
%             0   Call rng('shuffle'), i.e. generate new RANDOM state for random generator
%             <0  Call rng('shuffle'), i.e. generate new RANDOM state for random generator
%             []  Call rng('default'), i.e. use default random generator in default status
%             isnan(RS), do not call rng, i.e. keep current random generator in current status 
%             isstruct(RS)  Call rng(RS), i.e. initialize the random generator with a user given state
%
%
% OUTPUT:
% rngInit     Current rng state, except when isnan(RandState)
%
%

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2014-2014 by Tomlab Optimization Inc., $Release: 8.1.0$
% Written Oct 6, 2014.   Last modified Oct 6, 2014.

function rngInit = rngset(RandState)

try
% Set initial state for pseudo random generator
if isstruct(RandState)
   rng('default');
   rng(RandState);
elseif isnan (RandState)
   % Do no initialization
   % rng('default');
elseif isempty(RandState)
   rng('default');
elseif RandState > 0
   rng('default');
   rng(RandState);
elseif RandState == 0
   rng('default');
   rng('shuffle');
   %rng(RandState);
else  % RandState < 0
   rng('default');
   rng('shuffle');
end
catch
    if ~exist('rng','file')
        rngInit = tlRng(RandState);
        return
    else
        rethrow(lasterror);
    end
end

% Even though RandState is Nan or [], we save the current random generator (rng) status value.
if isnan(RandState)
   rngInit = nan;
else
   rngInit = rng;
end

% MODIFICATION LOG
%
% 141006 hkh Written
