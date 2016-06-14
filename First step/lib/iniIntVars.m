% Init routine. Called before solving an optimization problem
%
% The routine identifies the variables which are not continuously differentiable
% and therefore require a special treatment when estimating derivatives.
% The following fields are set in Prob.MIP:
%
% .CIV        Logical index vector of all integer variables which are
%             continuously differentiable.
%
% .SIV_f      Logical vector of size N x 1, where N is the number of variables.
%             SIV_f(i) == true indicates that the variable x(i) is not 
%             contiuouosly differentiable in the objective function.
%
% .SIV_c      Logical matrix of size m x N, where m is the number of nonlinear
%             constraints and N the number of variables. SIV_c(i,j) == true
%             indicates  that variables x(j) is not continuously differentiable
%             in the constraints function c(i).
%             
% .ConIntVars Logical index vector of all variables which are continuously
%             differentiable.
%
% function Prob = iniIntVars(Prob,DerLvlObj,DerLvlCon);
%
% DerLvlObj  Level of derivatives needed for the objective, 0,1,2
% DerLvlCon  Level of derivatives needed for the constraints, 0,1,2

% Bjorn Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2014 by Tomlab Optimization Inc., $Release: 8.1.0$
% Written November 17, 2014.  Last modified March 2, 2015.
function Prob = iniIntVars(Prob,DerLvlObj,DerLvlCon)

if nargin < 2
   DerLvlCon = 0;
   if nargin < 1
      DerLvlObj = 0;
   end
end

if ~any(checkType({'mip','miqp','miqq','minlp','glb','glc'},Prob.probType)) 
   return
end

% Dimensions
N = Prob.N;
mNonLin = Prob.mNonLin;

% Functions
f = DefPar(Prob.FUNCS,'f',[]);
g = DefPar(Prob.FUNCS,'g',[]);
H = DefPar(Prob.FUNCS,'H',[]);
c = DefPar(Prob.FUNCS,'c',[]);
dc = DefPar(Prob.FUNCS,'dc',[]);
d2c = DefPar(Prob.FUNCS,'H',[]);

% Small step for testing continuous differentiability
eps = 5.4774e-07;
optParam = DefPar(Prob,'optParam',[]);
if ~isempty(optParam)
   eps = DefPar(Prob.optParam,'DiffInt',eps);
end

% Error flags
f_ERROR = false; c_ERROR = false;

% Variable type indices
CIV = DefPar(Prob.MIP,'CIV',[]); % Continuously differentiable integer variables
SIV_f = DefPar(Prob.MIP,'SIV_f',[]); % Strict integer variables in objective function
SIV_c = DefPar(Prob.MIP,'SIV_c',[]); % Strict integer variables in constraint functions
IV = false(N,1); % All types of integer variables
if ~isempty(Prob.MIP.IntVars)
   IV(Prob.MIP.IntVars) = true;
end
CV = ~IV; % Continuous variables

% If both are empty, determine the classes
SIV_set = false(N,1);
if isempty(CIV) 
   if (isempty(SIV_f) || isempty(SIV_c))
      % find variables which give empty, NaN, or unchanged objective and 
      % constraint values, those are strict integers
      SIV_f = false(N,1);
      SIV_c = false(mNonLin,N);
      for xi = 1:N
          if IV(xi)
             % try calling objective function
             if ~isempty(f)
                [SIV_f(xi),ERROR] = isNotDifferentiable(Prob,xi,f,g,1,eps);
                f_ERROR = ERROR || f_ERROR;
             end
             % try calling constraints function
             if ~isempty(c) && mNonLin > 0
                [SIV_c(:,xi),ERROR] = isNotDifferentiable(Prob,xi,c,dc,mNonLin,eps);
                c_ERROR = ERROR || c_ERROR;
             end
          end
          if any([SIV_f(xi); SIV_c(:,xi)])
             SIV_set(xi) = true;
          end
      end
      CIV = ~SIV_set & IV;
   elseif isempty(SIV_f)
      CIV = IV & ~(any(SIV_c)');
    elseif isempty(SIV_c)
      CIV = IV & ~SIV_f;
   else
      CIV = IV & ~(SIV_f | any(SIV_c)');
   end
end
Prob.MIP.CIV = CIV;
Prob.MIP.SIV_f = SIV_f;
Prob.MIP.SIV_c = SIV_c;
Prob.MIP.ConIntVars = CV | CIV;

% Store original functions for use in callbacks
Prob.MIP.FUNCS = Prob.FUNCS;
Prob.MIP.FUNCS.narg = zeros(10,1);
% Prevent calls with noninteger values to functions returning errors
if f_ERROR
   Prob.FUNCS.f = 'nlp_f_iv';
   Prob.MIP.FUNCS.narg(1) = xnargin(Prob.MIP.FUNCS.f);
end
if f_ERROR
   Prob.FUNCS.c = 'nlp_c_iv';
   Prob.MIP.FUNCS.narg(4) = xnargin(Prob.MIP.FUNCS.c);
end

if ~isempty(SIV_f)
   % Change callback functions to support the strict variables
   if DerLvlObj > 0
      Prob.FUNCS.g = 'nlp_g_iv';
      if ~isempty(g)
         Prob.MIP.FUNCS.narg(2) = xnargin(g);
      end
      if DerLvlObj > 1
         Prob.FUNCS.H = 'nlp_H_iv';
         if ~isempty(H)
            Prob.MIP.FUNCS.narg(3) = xnargin(H);
         end
      end
   end
end
if ~isempty(SIV_c)
   if DerLvlCon > 0
      Prob.FUNCS.dc = 'nlp_dc_iv';
      if ~isempty(dc)
         Prob.MIP.FUNCS.narg(5) = xnargin(dc);
      end
      if DerLvlCon > 1
         Prob.FUNCS.d2c = 'nlp_d2c_iv';
         if ~isempty(d2c)
            Prob.MIP.FUNCS.narg(5) = xnargin(d2c);
         end
      end
   end
end

% Update patterns to avoid accidental noninteger evaluations in numerical
% differentiation routines
nsiv = sum(SIV_f);
if nsiv < N
   % update HessPattern
   Pattern = DefPar(Prob,'HessPattern',[]);
   if isempty(Pattern)
      % Estimate if a sparse pattern should be used
      nnz = N*N - (nsiv*N) - nsiv*(N-nsiv);
      if nnz < (N*(N - 1) - 1)/2
         ncd = N-nsiv;
         cdi = find(~SIV_f);
         i = repmat(cdi(:),ncd,1);
         j = repmat(cdi(:)',ncd,1);
         j = j(:);
         e = ones(length(i),1);
         Pattern = sparse(i,j,e,N,N);
      else
         Pattern = ones(N,N);
         Pattern(SIV_f,:) = 0;
         Pattern(:,SIV_f) = 0;
      end
   end
   Prob.HessPattern = Pattern;
end
anySIV_c = any(SIV_c);
nsiv = sum(anySIV_c);
if nsiv < N
   % update ConsPattern
   Pattern = DefPar(Prob,'ConsPattern',[]);
   if isempty(Pattern)
      % Estimate if a sparse pattern should be used
      Pattern = ~SIV_c;
      nnz = sum(sum(Pattern));
      if nnz < (mNonLin*(N - 1) - 1)/2
         Pattern = sparse(Pattern);
      end
   end
   Prob.ConsPattern = Pattern;
end

end

function [tf,ERROR] = isNotDifferentiable(Prob,xi,func,grad,func_size,eps)
   % Start with assuming that all functions are differentiable in the variable
   ERROR = false; % Pass back a flag if an error occured while evaluating the 
                  % function func using non-integer values.
   tf = false(func_size,1);
   x = get_x(Prob);
   x(Prob.MIP.IntVars) = round(x(Prob.MIP.IntVars));

   % A variable with equal lower and upper bound cannot be differentiated,
   % instead it should be removed from the variable set and replaced in all 
   % functions by a constant or fixed parameter.
   if Prob.x_L(xi) == Prob.x_U(xi) 
      tf = true(func_size,1); 
      return;
   end
   
   % First call to see if all values are returned when setting x to a noninteger
   x_org = x(xi);
   halfstep = 0.49999999; % If function rounds variables, 0.5 results in a whole 
                          % integer step. Avoid this.
   x(xi) = x(xi) + halfstep;
   if x(xi) >= Prob.x_U(xi)
      x(xi) = x(xi) - 2*halfstep;
   end
   try
      func_value = feval(func,x,Prob);
   catch
      ERROR = true;
      func_value = [];
   end
   if isempty(func_value)
       tf = true(func_size,1);
       return;
   end
   tf(isnan(func_value)) = true;

   % Second call to see if the function changes when taking a small step.
   x(xi) = x_org+eps;
   if x(xi) >= Prob.x_U(xi)
      x(xi) = x(xi) - 1;
      if(xi) <= Prob.x_L(xi)
         error(['iniIntVars: Variable x(' num2str(xi) ') should be a constant parameter ' ]);
      end
   end
   try
      func_value1 = feval(func,x,Prob);
   catch
      ERROR = true;
      func_value1 = [];
   end
   x(xi) = x_org;
   try
      func_value2 = feval(func,x,Prob);
   catch
      ERROR = true;
      func_value2 = [];
   end
   if isempty(func_value1) || isempty(func_value2)
      tf = true(func_size,1);
      return;
   end
   func_diff_smallstep = func_value1 - func_value2;
   func_diff_halfstep = func_value - func_value2;
   
   tf(isnan(func_value1) | isnan(func_value2)) = true;
   zerochange = (func_diff_smallstep == 0 & func_diff_halfstep == 0);
   
   % Check if a whole step changes the function values. If not, the variable
   % does not participate in the function
   if any(zerochange)
      if x(xi) - Prob.x_L(xi) < Prob.x_U(xi) - x(xi) 
         x(xi) = x(xi) + 1;
      else
         x(xi) = x(xi) - 1;
      end
      % Decision table:
      %
      %      step  half small integer  
      % func_diff  0    0     nz      => strict integer
      % ---||----  0    0     0       => not participating
      %
      func_value3 = feval(func,x,Prob);
      func_diff_integerstep = func_value - func_value3;
      isstrict = (func_diff_integerstep ~= 0) & zerochange;
      tf(isstrict) = true;
      % Note: The evaluation of 
      % >> (func_diff_integerstep == 0) && zerochange;
      % could be used to set zeros in HessPattern and ConsPattern
   end
   
   
   % Third call to see if all values are returned from gradient estimation
   if ~isempty(grad)
      try
         grad_value = feval(grad,x,Prob);
      catch
         grad_value = [];
      end
      if isempty(grad_value)
         tf = true(func_size,1); 
         return;
      end
      if isvector(grad_value)
         if isnan(grad_value(xi))
            tf = true;
            return
         end
      elseif length(size(grad_value))==2
         tf(isnan(grad_value(:,xi))) = true;
      else
         error(['Array of higher dimension than 2 returned from call to ' grad]);
      end
   end
end

function x = get_x(Prob)
   x = Prob.x_0;
   x_L = Prob.x_L;
   x_U = Prob.x_U;
   if isempty(x)
      infxL = isinf(x_L);
      infxU = isinf(x_U);
      x( infxL &&  infxU) = 0;
      x( infxL && ~infxU) = x_L + 1;
      x(~infxL &&  infxU) = x_U - 1;
      x( infxL &&  infxU) = (Prob.x_U - Prob.x_L);
      x = max(min(x_0,x_U),x_L);
   end
end

% MODIFICATION LOG
%
% 141117  bjo  Created.
% 150208  bjo  Correct initialization of CIV for different cases of SIV_f and
%              SIV_c.
% 150302  bjo  Check of problem type was disabled
% 150313  ango No ismatrix in older Matlab.