% nameprob.m
%
% Define different file names and corresponding description for the
% initialization files defining problems to be solved by TOMLAB.
% (The basic TOMLAB init (setup) files are normally named xxx_prob.m,
%  where e.g. xxx=con)
%
% function [F, N, D, probTypV] = nameprob (optType);
%
% INPUT:
%   optType Solver type number (See checkType for current list and numbers)
%     (UC, QP, CON, LS, LLS, Constrained LS (CLS), MIP, LP, Global UC (GLB), 
%     Global Con (GLC),MIQP, MINLP, SDP, BMI,CGO,ODE,EXP,MIQQ,GP)
%   mexType  Vector with allowed types of problems.
%     0 = Matlab, 1=AMPL, 2=CUTE, 3=NTS, 4=Helax many big files
%
% OUTPUT:
%   F        File names of all problem files
%   N        String matrix, used for the menu to choose problem files from
%   D        The default file for menu or driver routines, File Name F(D,:).
%   probTypV Vector of problem types, corresponding to:
%            UC(=1), QP(=2), CON(=3), LS(=4), LLS(=5), Constrained LS(=6),
%            MIP (=7), LP(=8), Global Unconstrained Optimization(=9),
%            Constrained Mixed-Integer Global Optimization(=10).
%            Mixed-integer quadratic programming (MIQP)  (=11)
%            Mixed-integer nonlinear programming (MINLP) (=12)
%            Semidefinite programming (SDP) Linear SDP with LMI  (=13)
%            Linear SDP with BMI constraints (BMI)  (=14)
%            Costly Global Optimization(=15).
%            ODEFit, ExpFit 
%            Mixed-Integer QP with Quadratic constraints,
%            Geometric Programming
%            One element for each problem file name.
%            See checkType for correct numbers

%            OC (=xx), Nonlinear time series(=xx).

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2014 by Tomlab Optimization Inc., $Release: 8.1.0$
% Written Feb 13, 1998.   Last modified Oct 2, 2014.

function [F, N, D, probTypV] = nameprob(optType,mexType)

if nargin < 2
   mexType=0:2;
   if nargin < 1
      optType=3;
   end
end

TomFile='TomlabProblem.mat';
load(TomFile,'tomProb');

D=1; % The default problem file is the first file by default

UC    = checkType('uc');
QP    = checkType('qp');
CON   = checkType('con');
LS    = checkType('ls');
LLS   = checkType('lls');
CLS   = checkType('cls');
MIP   = checkType('mip');
LP    = checkType('lp');
GLB   = checkType('glb');
GLC   = checkType('glc');
MIQP  = checkType('miqp');
MINLP = checkType('minlp');
SDP   = checkType('sdp');
BMI   = checkType('bmi');
CGO   = checkType('cgo');
ODE   = checkType('ode');
EXP   = checkType('exp');
NTS   = checkType('nts');
GP    = checkType('gp');
MCO   = checkType('mco');
LCP   = checkType('lcp');
MCP   = checkType('mcp');
OC    = checkType('oc');
MIQQ  = checkType('miqq');

switch optType
  case UC    % uc can handle uc, ls, and exp without constraints
    z=[UC, LS, ODE, EXP];
  case QP    % qp can handle qp and lp
    z=[QP, LP];
  case CON   % con handles many types, but not mi*, cgo, sdp/bmi, ode, nor
             % gp, mco, co, lcp, mcp, and nts
    z=[UC, QP, CON, LS, LLS, CLS, LP, GLB, GLC, ODE, EXP];
  case LS    % nlls (or ls) handles nlls, exp and nts
    z=[LS, ODE, EXP, NTS];
  case LLS   % lls 
    z=[LLS 100];
  case CLS   % cls handles cls, lls, ls, ode, exp and nts
    z=[CLS, LLS, LS, ODE, EXP, NTS];
  case MIP   % mip can handle only mip (+ LP)
    z=[MIP, LP];
  case LP    % lp can handle only lp
    z=[LP, 100];
  case GLB   % glb can handle uc, ls, ode, exp, and glb
    z=[GLB, UC, LS, ODE, EXP];
  case GLC   % glc can handle almost everything, not sdp, bmi, some minlp, 
             % nor gp, mco, co, lcp, mcp, and nts
    z=[UC, QP, CON, LS, LLS, CLS, MIP, LP, GLB, GLC, MIQP, MINLP, ODE, EXP];
  case MIQP  % miqp can handle mip (milp), qp and lp
    z=[MIQP, MIP, LP];
  case MINLP % minlp can handle almost everything, not sdp, bmi
    z=[MINLP, UC, QP, CON, LS, LLS, CLS, MIP, LP, GLB, GLC, MIQP, ODE, EXP];
  case SDP   % semidefinite programming can handle LP
    z=[SDP, LP];
  case BMI    % semidefinite programming with bmi can handle lp, sdp w lmi
    z=[BMI, LP, SDP];
  case CGO   % cgo can handle almost everything, not sdp, bmi, some minlp
    z=[CGO, UC, QP, CON, LS, LLS, CLS, MIP, LP, GLB, GLC, MIQP, MINLP, ODE, EXP];
  case ODE   % ode is not used for optType
    z=[ODE, 100];
  case EXP   % exp is not used for optType
    z=[EXP,100];
  case MIQQ  % miqq can handle mip (milp), qp and lp and miqp
    z=[MIQQ, MIQP, MIP, QP, LP];
  case GP    % gp can handle only gp
    z=[GP,100];
  case MCO   % mco can handle everything that nlp handles, but skip for now
    z=[MCO,100];
  case OC    % oc is not used for optType
    z=[OC,100];
  case LCP   % lcp
    z=[LCP,100];
  case MCP   % mcp 
    z=[MCP,100];
  case NTS   % nts is not used for optType
    z=[NTS,100];
end

n=length(z);
m=length(tomProb.probType);
n1=length(mexType);
if n1 == 1
   % Just to fool Matlab. Irregular behaviour for 1 element
   mexType=[mexType,100];
   n1=2;
end
m1=length(tomProb.mex);

ix=find(any( z(:)*ones(1,m)==(ones(n,1)*tomProb.probType)) & ...
        any(mexType(:)*ones(1,m1)==(ones(n1,1)*tomProb.mex)));

if ~isempty(ix)
   % Sort the problems of type SolvType first in list
   j=optType==tomProb.probType(ix);
   % Generate these problems backwards, so the latest defined is first
   jz=find(j);
   ix=ix([jz(length(jz):-1:1),find(~j)]);
   F=tomProb.File(ix,:);
   N=tomProb.Name(ix,:);
   probTypV=tomProb.probType(ix);
   D=1;
   for i=1:length(ix)
       % Check that the default file is really in the path
       % otherwise change the default file
       if exist(deblank(F(i,:)),'file')==2
          D=i;
          break;
       end
   end
   %if D~=1
   %   ix = [D,1:D-1,D+1:length(ix)];
   %   F  = F(ix,:);
   %   N  = N(ix,:);
   %   probTypV = probTypV(ix);
   %   D  = 1;
   %end
else
   F=[];
   N=[];
   probTypV=[];
end

% MODIFICATION LOG:
%
% 990727   hkh  Completely rewritten, using structures
% 990913   hkh  Modified to read tomProb from mat file
% 990913   hkh  Change the order of the default files, latest defined first.
% 990915   hkh  Safeguard against default file not in path
% 000820   hkh  Update for v3.0
% 000925   hkh  Change from solvType to optType
% 000927   hkh  mip was wrongly accessed for con, ls and cls.
% 001105   hkh  Insert lls, change exp and nts
% 020701   hkh  Add four new problem types
% 020702   hkh  Minor changes in what the types may handle
% 030117   hkh  Change miqq to bmi
% 030123   hkh  Correct the list of what different types can handle
% 050502   hkh  Move exp, add cgo, ode (nts, oc)
% 060602   hkh  Add mcp,lcp,mco,gp, make more general solution
% 141002   bjo  Replace all numerals representing problem types with 
%               variables obtained initially from checkType(strType)
