% TOMLAB gateway routine for computation of the objective function 
% to the MINLP feasibility problem
%
%      m
% min sum w_j * g_j(x,y_k)+ 
%     j=1
%
%  where y_k are fixed integer variables
%  and g_j(x,y_k)+ = max(0,g_j(x,y_k))
%  measures the constraint violation.
%
% function f = feas_f(x,Prob)
function f = feas_f(x,Prob)

global n_c NARG
global mad_c
global NLP_xc NLP_c  % Communication nlp_c/dc

NARG(4) = 2;

x = x(:);
Ax = [];
if Prob.mLin > 0
   Ax = Prob.A*x;
end
ProbNLP = Prob.user.ProbNLP;
w = DefPar(Prob.user,'w',[]);
f = 0;
b_L = Prob.b_L;
b_U = Prob.b_U;
c_L = Prob.c_L;
c_U = Prob.c_U;
mLin = Prob.mLin;
mNonLin = Prob.mNonLin;
if isempty(w)
   w = ones(Prob.mLin+Prob.mNonLin,1);
end
c = nlp_c(x,ProbNLP);
% NARG = []; NLP_xc=[]; NLP_c=[];
for i = 1:mLin
   if Ax(i) < b_L(i)
      f = f + w(i)*(b_L(i)-Ax(i));
   elseif Ax(i) > b_U(i)
      f = f + w(i)*(Ax(i)-b_U(i));
   end
end
for i = 1:mNonLin
   if c(i) < c_L(i)
      f = f + w(mLin+i)*(c_L(i)-c(i));
   elseif c(i) > c_U(i)
      f = f + w(mLin+i)*(c(i)-c_U(i));
   end
end
% Must reset to avoid crash in endSolve!
n_c = 0;

% MODIFICATION LOG
%
% 121213   bjo   Written for usage by LP/NLP-BB solver (stoaMINLP).