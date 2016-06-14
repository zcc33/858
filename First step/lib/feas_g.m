% TOMLAB gateway routine for computation of the gradient of the 
% objective function to the MINLP feasibility problem
%
%      m
% min sum w_j * g_j(x,y_k)+ 
%     j=1
%
%  where y_k are fixed integer variables
%  and g_j(x,y_k)+ = max(0,g_j(x,y_k))
%  measures the constraint violation.
%
% function g = feas_g(x,Prob)
function g = feas_g(x,Prob)

global n_c n_dc NARG
global mad_c mad_dc
global NLP_xc NLP_xdc NLP_c NLP_dc  % Communication nlp_c/dc

NARG(5) = 2;

x = x(:);
Ax = [];
if Prob.mLin > 0
   Ax = Prob.A*x;
end
ProbNLP = Prob.user.ProbNLP;
w = DefPar(Prob.user,'w',[]);
f_k = 0;
b_L = Prob.b_L;
b_U = Prob.b_U;
c_L = Prob.c_L;
c_U = Prob.c_U;
mLin = Prob.mLin; 
mNonLin = Prob.mNonLin;
N = length(x); % or Prob.N
if isempty(w)
   w = ones(Prob.mLin+Prob.mNonLin,1);
end
g = zeros(N,1);
c = nlp_c(x,ProbNLP);
dc = nlp_dc(x,ProbNLP);
for i = 1:mLin
   if Ax(i) < b_L(i)
      g = g - w(i)*A(i,:)';
   elseif Ax(i) > b_U(i)
      g = g + w(i)*A(i,:)';
   end
end
for i = 1:mNonLin
   if c(i) < c_L(i)
      g = g - w(mLin+i)*dc(i,:)';
   elseif c(i) > c_U(i)
      g = g + w(mLin+i)*dc(i,:)';
   end
end
% Must reset to avoid crash in endSolve!
n_c = 0;
n_dc = 0;

% MODIFICATION LOG
%
% 121213   bjo   Written for usage by LP/NLP-BB solver (stoaMINLP).