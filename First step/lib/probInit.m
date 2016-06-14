% probInit.m :
%
% General initialization routine for TOMLAB optimization problems,
% when a file in the TOMLAB Init File format is used
%
% function [Prob] = probInit (probFile, P, dims, Prob);
%
% INPUT:
% probFile Name of problem definition Init File, as a string, e.g. 'con_prob'.
% P        Problem number
%          If P=0, call Init File to get the predefined problems and let the user
%          choose a problem by displaying a menu
%          if isempty(P), just return the number of problems
%
% dims     Problems dimension, number of residuals and other variable parameters for some of
%          the predefined problems. 
%          if isempty(dims), default values are used.
%          Do help on the problem file, e.g. help mgh_prob, to see what dim values are possible
%
%          Example: Prob = probInit('glbv_prob',1,50);
%          Problem 1 in glbv_prob is defined with dimension 50
%
% Prob     Optional Tomlab problem structure with all parameters defining the problem
%
% OUTPUT:
% Prob     Problem structure or
%          Number of problems, if isempty(P)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1998-2011 by Tomlab Optimization Inc., Sweden. $Release: 7.9.0$
% Written June 16, 1999. Last modified Dec 15, 2011.

function Prob = probInit (probFile, P, dims, Prob)

if nargin < 4
   Prob=[];
   if nargin < 3
      dims=[];
      if nargin < 2
         P=[];
         if nargin < 1
            probFile='con_prob';
         end
      end
   end
end

if isempty(P)   % Only return the number of predefined problems
   probList=feval(probFile);
   nProbs=size(probList,1);
   Prob=nProbs;
   return
end

if P(1)<=0  % Give menu to choose problem from
   P=strmenu('Choice of test function ',feval(probFile,[]));
end


[probList, Prob] = feval(probFile, P, dims, Prob);

global probType

probType=Prob.probType;

Prob.probFile=probFile;

% MODIFICATION LOG:
%
% 980825  hkh  Changed error in comments. Changed call to usr_prob-routine.
% 980910  mbk  Do not set Prob to [] on line 72.
% 981026  hkh  Changed ProbFile to probFile. Set field Prob.probFile to probFile
% 990617  hkh  Change comments
% 990622  hkh  Add first argument (empty) to initial call to probFile,usr_prob
% 001012  hkh  Delete usr_prob possibility, not needed any longer
% 050228  hkh  Change prob.uP to Prob.uP
% 050302  hkh  Prob or uP as 4th input, to enable Prob.uP to be set
% 080624  hkh  Ask optionally the problem dimension, used e.g. by glbv_prob
% 111215  hkh  Removed parameter ask
% 111215  hkh  Revised. 3rd parameter dims setting dimension, # residuals etc
