% Implementation of preprocessing and probing techniques for mixed integer 
% programming problems on the TOMLAB format, based on the article 
% 'Preprocessing and Probing Techniques for Mixed Integer Programming Problems', 
% by M.W.P. Savelsbergh, School of Industrial and Systems Engineering, Georgia 
% Institute of Technology, Atlanta, GA 30332-0205
%
% BJORN - TODO1: Perform update of changes instantly,
%                and update L_min,L_max too!
% Future updates: Add phase 2 and 3.
%
% function [Prob,IFflag] = mipPP(Prob,PriLev)

% Written Sep 27, 2013.  Last modified Oct 24, 2013.
function [Prob,INFEASIBLE] = mipPP(Prob,PriLev)

INFEASIBLE = 0;

if nargin < 2
   PriLev = 0;
end

if nargin < 1
   warning('mipPP called with no input.');
   return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identify binary variables and nonbinary variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bI = ((Prob.x_L == 0) & (Prob.x_U == 1));
if max(Prob.MIP.IntVars) == 1
   iv = find(Prob.MIP.IntVars);
   bI = intersect(find(bI),find(iv));
else
   if length(Prob.MIP.IntVars) > 1
      bI = intersect(find(bI),Prob.MIP.IntVars);
   elseif length(Prob.MIP.IntVars) == 1
      bI = intersect(find(bI),1:Prob.MIP.IntVars);
   else
      bI = [];
   end
   iv = Prob.MIP.IntVars;
end
cI = setdiff(1:Prob.N,bI)';

%%%%%%%%%%%%%%%%%%%%%%%
% Reformulate problem %
%%%%%%%%%%%%%%%%%%%%%%%
b_i_u = find(~isinf(Prob.b_U));
b_i_l = find(~isinf(Prob.b_L));
b =  [Prob.b_U(b_i_u); -Prob.b_L(b_i_l)];
AG = [Prob.A(b_i_u,:); -Prob.A(b_i_l,:)];
xl = Prob.x_L;
xu = Prob.x_U;

if isempty(AG)
   % Nothing to do
   return
end

%%%%%%%%%%%%%
% Main loop %
%%%%%%%%%%%%%
phase1_i = 0;   % iteration index of phase one
impcoef_i = []; % keep track of constraints with changed coefficients
CHANGES = 1;    % terminate loop when there have been no changes
while CHANGES
   
   % first run basic preprocessing and probing
   [AG,b,xl,xu,impcoef_i,CHANGES,INFEASIBLE] = basicPP(b,b_i_l,b_i_u,AG,bI, ...
                                                       cI,xl,xu,iv, ...
                                                       impcoef_i,PriLev);
   if INFEASIBLE
      return
   end

   % BJORN TODO : identify logical implications and reduce system by probing 
   %              binary variables.
   
   phase1_i = phase1_i + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reassemble problem into tomlab format. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constraints
mLinMax = length(b_i_l)+length(b_i_u)+length(impcoef_i);
A = zeros(mLinMax,Prob.N);
b_L = zeros(mLinMax,1);
b_U = zeros(mLinMax,1);

% loop
bil = 1;
biu = 1;
bils = length(b_i_u);
b_i_l = [b_i_l;Inf];
b_i_u = [b_i_u;Inf];
i = 1;
while bil <= length(b_i_l)-1 || biu <= length(b_i_u)-1
   if b_i_l(bil) < b_i_u(biu)
      A(i,:) = -AG(bils+bil,:);
      b_L(i) = -b(bils+bil);
      b_U(i) = Inf;
      bil = bil + 1;
      i = i + 1;
   elseif b_i_l(bil) > b_i_u(biu)
      A(i,:) = AG(biu,:);
      b_L(i) = -Inf;
      b_U(i) = b(biu);
      biu = biu + 1;
      i = i + 1;
   else % b_i_l(bil) == b_i_u(biu)
      % check if coefficients changed
      if ismember(bils+bil,impcoef_i) || ismember(biu,impcoef_i)
         A(i,:) = -AG(bils+bil,:);
         b_L(i) = -b(bils+bil);
         b_U(i) = Inf;
         bil = bil + 1;
         i = i + 1;
         A(i,:) = AG(biu,:);
         b_L(i) = Inf;
         b_U(i) = b(biu);
         biu = biu + 1;
         i = i + 1;
      else
         A(i,:) = AG(biu,:);
         % should be indentical to AG(biu,:) !!!
         b_L(i) = -b(bils+bil);
         b_U(i) = b(biu);
         bil = bil + 1;
         biu = biu + 1;
         i = i + 1;
      end
   end
end

Prob.A = A(1:i-1,:);
Prob.b_L = b_L(1:i-1);
Prob.b_U = b_U(1:i-1);
Prob.mLin = size(Prob.A,1);

% Variable bounds
Prob.x_L = xl;
Prob.x_U = xu;

function [L_max,L_min] = Lmaxmin(AG,xl,xu,bI,cI)
  
xlc = xl(cI);
xuc = xu(cI);
AG_p = AG(:,cI);
AG_n = AG_p;
AG_p(~(AG_p > 0)) = 0;
AG_n(~(AG_n < 0)) = 0;
L_max = AG_p*xuc + AG_n*xlc;
L_min = AG_p*xlc + AG_n*xuc;
L_max = sum(max(AG(:,bI),0),2) + L_max;
L_min = sum(min(AG(:,bI),0),2) + L_min;

% BJORN TODO : Add optional outputs with info about performed techniques

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basicPP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic preprocessing and probing techniques
%
% Preprocessing
% - Identification of feasibility
% - Identification of redundancy
% - Improving bounds
% Probing
% - Fixing binary variables
% - Improving coefficients (binary variables)
function [AG,b,xl,xu,impcoef_i,CHANGES,INFEASIBLE] = basicPP(b,b_i_l,b_i_u, ...
                                                             AG,bI,cI,xl,xu, ...
                                                             iv,impcoef_i, ...
                                                             PriLev)

CHANGES = 0; INFEASIBLE = 0;

% Calculate L_i_max and L_i_min
[L_max,L_min] = Lmaxmin(AG,xl,xu,bI,cI);

if isempty(L_max) || isempty(L_min)
   return
end

% Identification of infeasibility
b_if = L_min > b;

% Identification of redundancy
b_re = L_max <= b;

% Remove any redundant constraints already
if any(b_re)
   AG = AG(~b_re,:);
   b = b(~b_re);
   L_max = L_max(~b_re);
   L_min = L_min(~b_re);
end

ncI = length(cI);
nb  = length(b);

% Improvement of variable bounds
% TODO : improve this method for largescale sparse matrices
%         no need to do calculations for zero-coefficients

% For all constraint coefficients on nonbinary variables, split coefficients
% into positive and negative ones.
Ac_p = AG(:,cI);
Ac_n = Ac_p;
Ac_p(Ac_p < 0) = 0; % Positive linear coefficients on nonbinary variables
Ac_n(Ac_n > 0) = 0; % Negative linear coefficients on nonbinary variables
B = b*ones(1,ncI);  % Create a matrix of upper bound vectors for each variable
LLmin = L_min*ones(1,ncI); % Create a matrix of L_min vectors for each variable
nb_p = (B-(LLmin-Ac_p*diag(xl(cI))))./Ac_p; % 
nb_n = -((LLmin-Ac_n*diag(xu(cI)))-B)./Ac_n;
nb_p(isinf(nb_p)) = NaN;
nb_n(isinf(nb_n)) = NaN;
if nb > 1
   nb_p = min(nb_p)';
   nb_n = max(nb_n)';
else
   nb_p = nb_p';
   nb_n = nb_n';
end
nb_p(nb_p+eps >= xu(cI)) = NaN;
nb_n(nb_n-eps <= xl(cI)) = NaN;
icv = ismember(cI,iv);
nb_p(icv) = floor(nb_p(icv)+eps);
nb_n(icv) = ceil(nb_n(icv)-eps);
impbnd = [nb_n nb_p];

% TODO: Fix this, something is wrong, test with ecpMINLP on Problem 10 to see
%       an infinite loop occur due to the improvement not being correctly
%       applied

% Used by fixing of variables and improvement of coefficients
if ~isempty(bI)
   nbI = length(bI);
   LLmax = L_max*ones(1,nbI);
   LLmin = L_min*ones(1,nbI);
   B = b*ones(1,nbI);
   nI = AG(:,bI) < 0;
end

% Fixing of variables
fixvar = [];
uf = xl(bI) ~= xu(bI);
if any(uf)
   fixvar = -ones(nbI,1);
   fI = LLmin + abs(AG(:,bI)) > B;
   fixvar(any(fI,1)) = max(nI(fI),[],1);
   fixvar(~uf) = -1;
end

% Improvement of coefficients
% - test with mip_prob 22,23,35,(46),47
impcoef = [];
if ~isempty(bI)
   % Bjorn TODO: handle Inf in L_max
   impcoef = zeros(nb,nbI);
   impcoef_append = zeros(nb,1);
   delta = B - LLmax + abs(AG(:,bI));
   dI = delta > 0;
   pI = AG(:,bI) > 0;
   nI = AG(:,bI) < 0;
   impcoef(dI & nI) = delta(dI & nI);
   impcoef(dI & pI) = -delta(dI & pI);
   impcoef_append(any(dI & pI,2)) = min(-delta(dI & pI),[],2);
   impcoef = [impcoef impcoef_append];
end

if any(b_if)
  if PriLev > 0
     disp('mipPP: infeasible constraints detected.');
  end
  INFEASIBLE = 1;
  return
end

if PriLev > 0
  if any(b_re)
     disp(['mipPP: redundant constraints detected in phase 1, iteration ' num2str(phase1_i)]);
  end
  if any(any(~isnan(impbnd)))
     disp(['mipPP: improved bounds in phase 1, iteration ' num2str(phase1_i)]);
  end
  if any(fixvar > -1)
     disp(['mipPP: fixed binary variables in phase 1, iteration ' num2str(phase1_i)]);
  end
  if any(any(impcoef))
     disp(['mipPP: improved linear constraint coefficients in phase 1, iteration ' num2str(phase1_i)]);
  end
end

impcoef_i = union(impcoef_i,find(any(impcoef,2)));

if any(b_re)
  bil_i = b_re(length(b_i_u)+1:end);
  biu_i = b_re(1:length(b_i_u));
  b_i_l = b_i_l(~bil_i);
  b_i_u = b_i_u(~biu_i);
end

% Apply any changes

% Remove any redundant constraints
% if any(b_re)
%   AG = AG(~b_re,:);
%   b = b(~b_re);
%   CHANGES = 1;
% end

% Improved bound on variables
if any(any(~isnan(impbnd)))
  xl(cI) = max(impbnd(:,1),xl(cI));
  xu(cI) = min(impbnd(:,2),xu(cI));
  CHANGES = 1;
end

% Fixing of binary variables
fix_i = fixvar > -1;
if any(fix_i)
  fix_i = find(fix_i);
  xl(bI(fix_i)) = fixvar(fix_i);
  xu(bI(fix_i)) = fixvar(fix_i);
  CHANGES = 1;
end

% Improved constraint coefficients
if any(any(impcoef))
  AG(:,bI) = AG(:,bI) + impcoef(:,1:end-1);
  b = b + impcoef(:,end);
  CHANGES = 1;
end

% MODIFICATION LOG
%
% 130927 bjo  Written. Added identification of feasibility, redundancy and
%             improvements of bounds.
% 131002 bjo  Added fixing of variables and improvement of coefficients.
% 131022 bjo  Reformulation of problem back to TOMLAB format.
%             Optional print of any basic preprocessing performed added.
% 131101 bjo  Improved code efficiency.
% 141208 bjo  Removed old code. Fixed variable-bound updates for problems with
%             only one linear constraint
% 141211 bjo  Merged basicPP with applyPP.
% 150120 bjo  Break immediately when infeasible constraint detected.