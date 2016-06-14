% SolverList returns a list of all solvers for a particular solvType
%
% function [SolvList,SolvTypeList] =  SolverList(solvType,LargeScale, Silent)
%
% INPUT:
% solvType   The TOMLAB solvType/optType number (See help checkType for the
%            numbers), or the optType name
% optType name  Description
%      uc       Unconstrained Optimization (UC), bounds.
%      qp       Quadratic Programming (QP)
%      con      Constrained Nonlinear Programming (NLP)
%      ls       Nonlinear Least Squares (NLLS), bounds.
%      lls      Linear Least Squares (LS)
%      cls      Constrained Nonlinear Least Squares
%      mip      Mixed-Integer Programming
%      lp       Linear Programming
%      glb      Global optimization (GO), bounds.
%      glc      Global optimization (GO), constraints.
%      miqp     Mixed-Integer Quadratic Programming (MIQP)
%      minlp    Mixed-Integer Nonlinear Programming (MINLP)
%      sdp      Semidefinite Programming (SDP), LMI
%      bmi      Linear SDP with BMI constraints
%      cgo      Costly Global optimization (CGO)
%      ode      Parameter estimation in ODE (ODEFit)
%      exp      Exponential sum fitting (ExpFit)
%      miqq     Mixed-Integer QP with Quadratic constraints
%      gp       Geometric Programming Problems
%      mco      Multi-Criteria Optimization
%      oc       Optimal control
%      lcp      Standard Linear Complementarity Problem (LCP)
%      mcp      Polyhedrally constrained variational
%               inequality Problem or
%               Mixed Complementarity Problem (MCP)
%      nts      Nonlinear Time Series
%
% LargeScale 1 = Large Scale problems, 0 = Small or medium sized
%
% Silent  If Silent == 1, totally silent
%         If Silent == 0 (default), an information text is displayed
%
% OUTPUT:
% SolvList     All solvers for a particular solvType
% SolvTypeList The solvType numbers corresponding to the elements in SolvList

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2014 by Tomlab Optimization Inc., $Release: 8.1.0$
% Written Nov 29, 1998.   Last modified Oct 8, 2014.
function [SolvList,SolvTypeList] =  SolverList(solvType,LargeScale, Silent)

if nargin < 3
   Silent = 0;
   if nargin < 2
      LargeScale = 0;    
      if nargin < 1
         solvType = [];
      end
   end
end

[TomV,os,TV]=tomlabVersion;

% How to update with new optType
%
% Add new 2-3 letter code to optType
% Add full name to optTypeDescr
% Add shorter description to optTypeDescrShort
% Set nToUse to the number of types to actually display
nToUse = 19;
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
MIQQ  = checkType('miqq');
% NTS   = checkType('nts');
GP    = checkType('gp');
MCO   = checkType('mco');
% LCP   = checkType('lcp');
MCP   = checkType('mcp');
OC    = checkType('oc');

optType(UC,:)     = 'uc   ';
optType(QP,:)     = 'qp   ';
optType(CON,:)    = 'con  ';
optType(LS,:)     = 'ls   ';
optType(LLS,:)    = 'lls  ';
optType(CLS,:)    = 'cls  ';
optType(MIP,:)    = 'mip  ';
optType(LP,:)     = 'lp   ';
optType(GLB,:)    = 'glb  ';
optType(GLC,:)    = 'glc  ';
optType(MIQP,:)   = 'miqp ';
optType(MINLP,:)  = 'minlp';
optType(SDP,:)    = 'sdp  ';
optType(BMI,:)    = 'bmi  ';
optType(CGO,:)    = 'cgo  ';
optType(ODE,:)    = 'ode  ';
optType(EXP,:)    = 'exp  ';
optType(MIQQ,:)   = 'miqq ';
optType(GP,:)     = 'gp   ';
optType(MCO,:)    = 'mco  ';
optType(MCP,:)    = 'mcp  ';
optType(OC,:)     = 'oc   ';
optType(end+1,:)  = '     ';

optTypeDescr(UC,:)    = 'Unconstrained Optimization (UC), bounds    ';
optTypeDescr(QP,:)    = 'Quadratic Programming (QP)                 ';
optTypeDescr(CON,:)   = 'Constrained Nonlinear Programming (NLP)    ';
optTypeDescr(LS,:)    = 'Nonlinear Least Squares (NLLS), bounds     ';
optTypeDescr(LLS,:)   = 'Linear Least Squares (LS)                  ';
optTypeDescr(CLS,:)   = 'Constrained Nonlinear Least Squares        ';
optTypeDescr(MIP,:)   = 'Mixed-Integer Linear Programming           ';
optTypeDescr(LP,:)    = 'Linear Programming                         ';
optTypeDescr(GLB,:)   = 'Global Optimization (GO), bounds           ';
optTypeDescr(GLC,:)   = 'Global Optimization (GO), constraints      ';
optTypeDescr(MIQP,:)  = 'Mixed-Integer Quadratic Programming (MIQP) ';
optTypeDescr(MINLP,:) = 'Mixed-Integer Nonlinear Programming (MINLP)';
optTypeDescr(SDP,:)   = 'Semidefinite Programming (SDP)             ';
optTypeDescr(BMI,:)   = 'Linear SDP with BMI constraints            ';
optTypeDescr(CGO,:)   = 'Costly Global Optimization (CGO)           ';
optTypeDescr(ODE,:)   = 'Parameter estimation in ODE (ODEFit)       ';
optTypeDescr(EXP,:)   = 'Exponential sum fitting (ExpFit)           ';
optTypeDescr(MIQQ,:)  = 'Mixed-Integer QP w Quadratic constr (MIQQ) ';
optTypeDescr(GP,:)    = 'Geometric Programming (GP)                 ';
optTypeDescr(MCO,:)   = 'Multi-Criteria Optimization (MCO)          ';
optTypeDescr(MCP,:)   = 'Var InEq or Mixed Complementarity   (mcp)  ';
optTypeDescr(OC,:)    = 'Optimal control (OC)                       ';
optTypeDescr(end+1,:) = '                                           ';

optTypeDescrShort(UC,:)      = 'UC                        ';
optTypeDescrShort(QP,:)      = 'QP                        ';
optTypeDescrShort(CON,:)     = 'NLP                       ';
optTypeDescrShort(LS,:)      = 'NLLS                      ';
optTypeDescrShort(LLS,:)     = 'LS                        ';
optTypeDescrShort(CLS,:)     = 'Constrained Nonlinear NLLS';
optTypeDescrShort(MIP,:)     = 'Mixed-Integer Programming ';
optTypeDescrShort(LP,:)      = 'Linear Programming        ';
optTypeDescrShort(GLB,:)     = 'GO, bounds                ';
optTypeDescrShort(GLC,:)     = 'GO, constraints           ';
optTypeDescrShort(MIQP,:)    = 'MIQP                      ';
optTypeDescrShort(MINLP,:)   = 'MINLP                     ';
optTypeDescrShort(SDP,:)     = 'SDP /LMI                  ';
optTypeDescrShort(BMI,:)     = 'Linear SDP /BMI           ';
optTypeDescrShort(CGO,:)     = 'CGO                       ';
optTypeDescrShort(ODE,:)     = 'ODEFit                    ';
optTypeDescrShort(EXP,:)     = 'ExpFit                    ';
optTypeDescrShort(MIQQ,:)    = 'MIQQ                      ';
optTypeDescrShort(GP,:)      = 'Geometric Programming     ';
optTypeDescrShort(MCO,:)     = 'Multi-Criteria Opt        ';
optTypeDescrShort(MCP,:)     = 'LCP, MCP                  ';
optTypeDescrShort(OC,:)      = 'Optimal Control           ';
optTypeDescrShort(end+1,:)   = '                          ';

% How to update with new solver
%
% Add new name to SolveList
% Add solvertype to SolvTypeList
% Add license number to whichTV
SolvList=str2mat(...
     'ucSolve','fmins','fminu' ...
    ,'qpe','qplm','qpSolve','qpBiggs','qpopt','qp','quadprog','qp-minos'...
    ,'qld' ...
    ,'nlpSolve','conSolve','sTrustr','constr','minos','npsol','PDCO'...
    ,'snopt','fmincon'...
    ,'leastsq','lsqnonlin' ...
    ,'lsei', 'lssol' ...
    ,'clsSolve','nlssol' ...
    ,'mipSolve','cutplane' ...
    ,'lpSimplex','akarmark','lp','linprog','lpsimp2','lpopt','lp-minos' ...
    ,'glbSolve','ego','glbFast','glcSolve','glcFast','glcCluster' ...
    ,'rbfSolve','bqpd','miqpBB','minlpBB','PENSDP','filterSQP','CPLEX' ...
    ,'PENBMI','PDSCO','nlpsolv' ...
    ,'knitro','conopt','nlpqlp','nlpjob','dfnlp','lgo' ...
    ,'misqp','miql','ql' ...
    ,'DualSolve','slsSolve','sqopt','modfit','coplgp','path' ...
    ,'multiMin','multiMinlp','ecpMINLP','stoaMINLP' ...
    ,'milpSolve','arbfmip');

% 70 are added to the solvers of less quality
SolvTypeList=[...
      UC UC UC ...
      QP+70 QP+70 QP QP+70 QP QP QP QP ...
      QP ...
      CON CON CON CON CON CON CON ...
      CON CON ...
      LS LS ...
      LLS LLS ...
      CLS CLS ...
      MIP MIP ...
      LP LP+70 LP LP LP+70 LP LP ...
      GLB CGO GLB GLC GLC GLC ...
      CGO QP MIQP MINLP SDP CON MIQP ...
      BMI CON CON+70 ...
      CON CON CON MCO CLS GLC ...
      MINLP MIQP QP ...
      LP CLS QP ODE GP MCP ...
      MINLP MINLP MINLP MINLP ...
      MIP CGO];

whichTV=[1 1 1 ...
         1 1 1 1 2 1 1 2 ...
         1 ...
         1 1 1 1 2 4 1 ...
         4 1 ...
         1 1 ...
         1 3 ...
         1 3 ...
         1 1 ...
         1 1 1 1 1 2 2 ...
         1 5 1 1 1 1 ...
         5 7 7 7 6 7 9 ...
         10 1 1 ...
         11 12 16 16 16 17 ...
         19 19 19 ...
         1 1 4 16 25 20 ...
         29 29 29 29 ...
         1 5];

if isempty(solvType) % No input argument was given
%   SolvList = SolvList(SolvTypeList < 70,:);
   if ~Silent

      SolvList = [];
      % Print information of possible input arguments
      fprintf('\nPlease give the TOMLAB optType name ');
      fprintf('or the solvType number, currently 1-%d, as input\n\n',nToUse);
      fprintf('optType  solvType  Description\n');
      for t = 1:nToUse
          fprintf(optType(t,:));
         for s = 1:(9-size(optType(t,:),2))
            fprintf(' ');
         end
         if t < 10
            fprintf('%d         ',t);
         else
            fprintf('%d        ',t);   
         end
         fprintf(optTypeDescr(t,:));
         fprintf('\n');
      end
      fprintf('\n');
   end
   return;
end

i = [];
if ischar(solvType)
   % Find recommended solver
   RecSolvLarge=GetSolver(solvType,1);
   RecSolvSmall=GetSolver(solvType,0);
   % Convert string to the solvType number
   solvType = checkType(solvType);
   if isempty(solvType)
      % Given solvType was given as a string but it did not match
      % any existing types. See checkType for the correct list
      SolvList = [];
      SolvTypeList = [];
      return;
   end
else
   solvType = min(nToUse,max(1,solvType));
   solvString = deblank(optType(solvType(1),:));
   %RecSolvLarge=[];
   %RecSolvSmall=[];
   %for i=1:size(solvstring,1)
   %    RecSolvLarge=[RecSolveLarge,GetSolver(solvString(i,:),1)];
   %    RecSolvSmall=[RecSolveSmall,GetSolver(solvString(i,:),0)];
   %end
   RecSolvLarge=GetSolver(solvString,1);
   RecSolvSmall=GetSolver(solvString,0);
end

% Read probTypeList that gives the list of what other types of problem
% a solver of type solvType will handle

ix = [];
for i=1:nToUse
    if i ~= solvType
       [DataFile,NameFile,DefFile,probTypeList]=nameprob(i,0);
       CanSolve = any(solvType==probTypeList);
       if CanSolve
          % Add the type i to the list of the types that can solve solvType
          ix = [ix i];
       end
   end
end

% Make a list of solvers of the types in ix,
% also able to solve solvType problem
iopt = [];
for i = 1:size(ix,2)
   iopt=cat(2,iopt,find(ix(i)==SolvTypeList));
end
OthSolvList=SolvList(iopt,:);
OthSolvTV=whichTV(iopt);

% Make a list of the solvers for the particular solvType
ipart=find(solvType==SolvTypeList);
PartSolvList=SolvList(ipart,:);
PartSolvTypeList=SolvTypeList(ipart);
PartSolvTV=whichTV(ipart);

% Check which ones there is a license for
% Check TV vector

SolvList = PartSolvList;
SolvTypeList = PartSolvTypeList;
if Silent, return; end

% Print recommended choices of solvers

if LargeScale
   fprintf('\nTomlab recommended choice for large scale ');
   fprintf(optTypeDescr(solvType,:));
   fprintf('\n\n');
   fprintf(RecSolvLarge);
   fprintf('\n');
end
fprintf('\nTomlab recommended choice for ');
if LargeScale
   fprintf('small scale ');
end
fprintf(optTypeDescr(solvType,:));
fprintf('\n\n');
fprintf(RecSolvSmall);
fprintf('\n\n');

fprintf('Other solvers for ');
fprintf(deblank(optTypeDescrShort(solvType,:)));
fprintf('\n\n');

fprintf('   Licensed:\n\n   ');
count = 0;
for i = 1:size(PartSolvList)   
   if TV(PartSolvTV(i)) & ...
   ~(strcmpi(deblank(PartSolvList(i,:)),RecSolvSmall) | ...
    (strcmpi(deblank(PartSolvList(i,:)),RecSolvLarge) &...
      LargeScale))
      fprintf(PartSolvList(i,:));
      fprintf('\n   ');
      count = 1;
   end
end
if ~count
   fprintf('NONE\n');
end
fprintf('\n');

fprintf('   Non-licensed:\n\n   ');
count = 0;
for i = 1:size(PartSolvList)   
   if ~TV(PartSolvTV(i))& ...
   ~(strcmpi(deblank(PartSolvList(i,:)),RecSolvSmall) | ...
    (strcmpi(deblank(PartSolvList(i,:)),RecSolvLarge) &...
      LargeScale))
      fprintf(PartSolvList(i,:));
      fprintf('\n   ');
      count = 1;
   end
end
if ~count
   fprintf('NONE\n');
end
fprintf('\n');

% Print other solvers also capable of solving problem type

fprintf('Solvers also handling ');
fprintf(optTypeDescrShort(solvType,:));
fprintf('\n\n');

fprintf('   Licensed:\n\n   ');
count = 0;
for i = 1:size(OthSolvList)   
   if TV(OthSolvTV(i))
      fprintf(OthSolvList(i,:));
      fprintf('\n   ');
      count = 1;
   end
end
if ~count
   fprintf('NONE\n');
end
fprintf('\n');

fprintf('   Non-licensed:\n\n   ');
count = 0;
for i = 1:size(OthSolvList)   
   if ~TV(OthSolvTV(i))
      fprintf(OthSolvList(i,:));
      fprintf('\n   ');
      count = 1;
   end
end
if ~count
   fprintf('NONE\n');
end
fprintf('\n');

% MODIFICATION LOG
%
% 981129  hkh  Written
% 981203  mbk  Check if isempty(solvType) then solvType was given as input
%              argument but it did not match any of 'uc','qp', ...
% 981209  hkh  Add lpsimp2 among LP solvers
% 990906  hkh  Add glcSolve and glcRun
% 990913  hkh  Add MIP as number 7, and use use exp, not ef for number 5.
% 001004  hkh  Adding lpopt
% 010714  hkh  Adding xpress-mp and qp-xpress-mp
% 010715  hkh  Adding glbFast
% 010815  hkh  Adding glcFast
% 011111  hkh  Adding glcCluster and rbfSolve
% 020701  hkh  Adding six new solvers
% 020702  hkh  Update optType with new problem types
% 020708  bjo  Change output to recommended solvers + other solvers.
%              Also print which solvers are licensed and which are not.
% 030117  hkh  Change type miqq to bmi
% 030123  hkh  Correct the list whichTV, add PDCO, PDSCO
% 030213  ango Correct names (Dundee) and filterSQP type
% 040126  hkh  Now only returns the solvers for the particular type
% 040126  hkh  Added flag Silent, if true totally silent
% 040126  hkh  Empty argument displays all available solvers
% 040413  med  Added more solvers
% 040414  hkh  Correct ego, now type glc (10)
% 050502  hkh  Add CGO,ODE,OC. Only display nToUse many
% 050602  hkh  Add GP. Make more general solution, using checkType
% 050801  med  isstr replaced by ischar
% 060212  hkh  Add arbfmip
% 090813  med  mlint check
% 141003  bjo  Remove xa, xpress, oqnlp, socs and dido.
% 141003  bjo  Add misqp, miql, ql, multimin, multiminlp, ecpminlp
%              and stoaminlp.
% 141006  bjo  Replace all numerals representing problem types with 
%              variables obtained from checkType(strType).
% 141007  bjo  Remove call to strmatch, just use checkType(solvType).
% 141008  bjo  Simplify help text.