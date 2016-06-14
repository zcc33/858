% function rngInit = tlRng(RandState)
%
% INPUT:
% RandState   Value RS determines how to obtain random numbers using the pseudo-random generator
%          RS >0  Call rng(RS), i.e. reinitialize random generator with value RS
%             0   Call rng('shuffle'), i.e. generate new RANDOM state for random generator
%             <0  Call rng('shuffle'), i.e. generate new RANDOM state for random generator
%             []  Call rng('default'), i.e. use default random generator in default status
%             isnan(RS), do not call rng, i.e. keep current random generator in current status 
%             isstruct(RS)  Call rng(RS), i.e. initialize the random generator with a user given state

% Anders Goran, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2014-2014 by Tomlab Optimization Inc., $Release: 8.1.0$
% Written Oct 2, 2014.  Last modified Oct 16, 2014.

function rngInit = tlRng(RandState)

% R2011A introduced rng. Try to replicate proper behaviour through rand() on older
% Matlab versions. 

% No inputs - return state and exit.
if nargin < 1
   try
      rngInit = rng;
      return
   catch
      if ~exist('rng')
         rngInit = rand('state');
         return
      else
         rethrow(lasterror);
      end
   end
   return
end

try 
   % This is mostly like hkh's rngset but the isstruct case must come 
   % first: isnan on a struct throws an error.
   if isstruct(RandState)
      rng('default');
      rng(RandState);
      rngInit = rng;
   elseif isnan(RandState)
      % Do no initialization. Return NaN as rng status value.
      % rng('default');
      rngInit = NaN;
   elseif isempty(RandState)
      rng('default');
      rngInit = rng;
   elseif RandState > 0
      rng(RandState);
      rngInit = rng;
   elseif RandState == 0
      %rng(RandState);
      rng('default');
      rng('shuffle');
      rngInit = rng;
   else  % RandState < 0
      rng('default');
      rng('shuffle');
      rngInit = rng;
   end
   
   return
   
catch

   % Old behaviour if rng doesn't exist
   if ~exist('rng','file')
      if isnan(RandState)
         % Do no initialization
         rngInit = NaN;
      elseif length(RandState)>1 
         rand('state',RandState);
         rngInit = rand('state');
      elseif ~isempty(RandState)
         if RandState >= 0
           rand('state',RandState);
           rngInit = rand('state');
         else
           rngInit = rand('state');
         end
      else
         rand('state',sum(100*clock));
         rngInit = rand('state');
      end
      return
   else
      % Something else is wrong, give up at last.
      rethrow(lasterror);
   end
end

% MODIFICATION LOG
%
% 141002 ango Written
% 141016 ango Adapt to rngset 
