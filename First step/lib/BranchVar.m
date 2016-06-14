% --------------------------------------------------------------------------
function [IntVar,xC,xI,UPC] = BranchVar(VarSel,x,Ridx,VarW,eps_I,...
                    PriLev,UPC,PCU,PCL,NODE)
% --------------------------------------------------------------------------
% BranchVar implements variable selection strategies in Branch and Bound algorithms
%
% INPUT:
%   VarSel    Variable selection strategy
%             = 1: Simple fractional strategy to find variable to branch on
%             = 2: Use Lagrange multipliers, distance to nearest integer value 
%             = 3: Pseudo-cost based
%   x         Current solution x point
%   Ridx      The integer variables that are currently real-valued in solution x
%   VarW      Weight for each variable in the variable selection phase.
%             A lower value gives higher priority.
%   eps_I     Tolerance how close to integer is OK
%   PriLev    Print level
%
% Inputs required only for pseudocost branching
%
%   UPC       List of nodes for which to update pseudocosts after next
%             solution of the relaxed outer-approximation based master 
%             problem
%   PCU, PCL  Pseudocosts for rounding each variable up or down.
%   NODE      Number of node in branch and bound tree
%
% OUTPUT:
%   IntVar    Flag if x is already integer valued
%   xC        Index for the variable selected
%   xI        The lower integer value in the branching: x(xC) <= xI and x(xC) >= xI+1
%   UPC       Updated input UPC.
%   PCU,PCL   Updated pseudocosts.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Oct 15, 2009.   Last modified Oct 16, 2009.

% Avoid returning empty xC
x_I = floor(x(Ridx));
% x_I = floor(x(Ridx)+eps_I*max(1,abs(x(Ridx))));
x_r = max(0,x(Ridx)-x_I);

switch VarSel
  case 1
    % VarSel = 1: Simple fractional strategy to find variable to branch on
    
    if isempty(VarW)
       % Variables with frac.part closest to 0.5
       %[best_frac xBest] = min(abs(x_r-0.5));
       r             = abs(x_r-0.5);
       [rBest iBest] = sort(r);
       irBest        = find(rBest < 0.5-eps_I);
    else
       r             = VarW(Ridx);
       r(x_r < eps_I*max(1,x_r))=Inf;
       [rBest iBest] = sort(r);
       irBest        = find(~isinf(rBest));
    end

    if isempty(irBest)
       IntVar = 1;
       xC     = 0;
       xI     = 0;
    else
       IntVar = 0;
       iCand  = iBest(irBest);
       Cand   = Ridx(iCand);

       if PriLev > 1
          fprintf('   %d nonintegers. ',length(irBest));
          if isempty(VarW)
             fprintf('Best fraction = %8.6f.',x_r(iCand(1)))
          else
             fprintf('Weight %8.3f. Fraction %8.6f.',r(iCand(1)),x_r(iCand(1)))
          end
          fprintf(' Base var = %4.0f. Var # = %4.0f',iCand(1),Cand(1))
          fprintf('\n')
       end
       xI = x_I(iCand(1));
       xC = Cand(1);
    end
  case 2
    % HKH This strategy is not OK yet
    % VarSel = 2: Use Lagrange multipliers, distance to nearest integer value 
    v           = (1-x_r).*VarW(Ridx);
    y           = -x_r.*VarW(Ridx);
    z           = min(v,y);
    if PriLev > 1
       xprint(v,'v:')
       xprint(y,'y:')
       xprint(z,'z:')
    end

    [xMin iMin] = min(z);
    xI          = x_I(iMin);
    xC          = Ridx(iMin);
    IntVar      = 0;
    if abs(x(Ridx)-xI) <= eps_I
       IntVar = 1;
    else
       IntVar = 0;
    end
    if PriLev > 1
       fprintf('xMin %f iMin %d, xI %d xC %d\n',xMin,iMin,xI,xC);
    end
  case 3
    % Pseudo-cost branching
    xI = floor(x(Ridx));
    xFrac = x(Ridx)-xI;
    D = PCL(Ridx).*xFrac + PCU(Ridx).*(1-xFrac);
    [ds di] = sort(D(:),1,'descend');
    xC = Ridx(di(1));
    % Save selected lower branching bound
    xI = xI(di(1));
    % Update pseudo costs after next CMP solution of these nodes
    UPC = [UPC NODE+1 NODE+2];
    if PriLev > 1
       fprintf('max(PCL*xf+PCU*(1-xf)) %f ipcMax %d, xI %d xC %d\n',ds(1),di(1),xI,xC);
    end
    if abs(xFrac(di(1))) <= eps_I
       IntVar = 1;
    else
       IntVar = 0;
    end
end

% MODIFICATION LOG
%
% 091015  hkh  Written
% 130717  bjo  Added pseudocost branching
