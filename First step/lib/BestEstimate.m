function E = BestEstimate(fmin, L, Lx, PCL, PCU, xIP)
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Best estimate, estimation of best solution obtainable from node in 
% mixed-integer branch-and-bound tree. Uses pseudo-costs.
%
% INPUT:
%
%   fmin       Lower bounds on integer solution
%              fmin must have length [maxNodes] and be indexed by values in L.
%   L          List of currently active nodes. 
%              L must have length [nNodes].
%   Lx         Integer variables values from latest relaxed LP solution
%              of the parent nodes to the nodes in L.
%              Lx must have size [nIntVars x nNodes].
%   PCL, PCU   Pseudo costs for integer variables. 
%              PCL and PCU must have length [nIntVars].
%   xIP        Integer values from one or more feasible integer solution, 
%              optional.
%              xIP must have size [nIntVars x nIntSolutions].
%
% OUTPUT:
%
%   E          Best estimate of integer solution for nodes in L. 
%              E is of the same size as L.
% 
if nargin < 5
   error('BestEstimate : Need at least 5 input arguments');
end

nI = size(Lx,1);
if nargin < 6
   xIP = [];
else
   if size(xIP,1) ~= size(Lx,1);
      error('BestEstimate : Illegal size of xIP, must have nI rows');
   end
end

nNodes = length(L);
E = zeros(nNodes,1);
for i = 1:nNodes
    f = Lx(:,i) - floor(Lx(:,i));
    if isempty(xIP)
       % simple optimistic estimate
       E(i) = fmin(L(i)) + sum(min(PCL.*f,PCU.*(1-f)));
    else
       % extended estimate using probability of round-off to integer solution, 
       % requires at least one already obtained feasible integer solution
      
       % calculate flip-percentage
       fp = sum(abs(Lx - xIP(:,1)) > 0.5)/size(xIP,1);
      
       nIP = size(xIP,2); % Number of solutions found
       wIP = 1-1/(nIP+1); % Weight given to obtained IP solutions
       % calculate q_j and E
       E(i) = fmin(L(i));
       for j = 1:nI
           q_j = 1-4*fp*min(f(j),1-f(j));
           xI = round(Lx(j));
           pxI = sum(xIP(j,:) == xI)/nIP;
           q_j = wIP*pxI + (1-wIP)*q_j;
           if f(j) <= 0.5
              E(i) = E(i) + f(j)*PCL(j)*q_j + (1-f(j))*PCU(j)*(1-q_j);
           else
              E(i) = E(i) + f(j)*PCL(j)*(1-q_j) + (1-f(j))*PCU(j)*q_j;
           end
       end
    end
end

% MODIFICATION LOG
%
% 130714  bjo  Written, separated from lp-nlp-bb minlp solver
% 130821  bjo  Fixed file description, extended input information