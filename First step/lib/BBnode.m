% -------------------------------------------------------------------------
function [Node,L,Lix] = BBnode(NodeSelMethod,BacktrackMethod,BACKTRACK,L,...
                               pred,Depth,Icomp,nI,f_min,fIPMin,Lx,s0,E)
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% BBnode implements node selection strategies in Branch and Bound algorithms
%
% INPUT:
%   NodeSelMethod    Node selection method in branch and bound
%           = 1 Depth First. Priority on nodes with more integer components.
%           = 2 Breadth First. Priority on nodes with more integer components.
%           = 3 Pure LIFO (Last in, first out) Depth First
%           = 4 Pure FIFO (First in, first out) Breadth First
%           = 5 Maximize integer components.
%           = 6 Best bound
%           = 7 Best estimate using pseudo-costs
%           = 8 Best projection
%   BacktrackMethod  Node selection used when backtrack criterion fulfilled. 
%              Methods available are the same as listed for NodeSelMethod above.
%   BACKTRACK Flag, set to nonzero when the backtrack criterion is fulfilled.
%   L         Active node list
%   pred      Preceding node in the tree for each node
%   Depth     Depth in the tree for each node
%   Icomp     Number of integer components in each computed node
%   nI        Number of integer variables in problem
%   f_min     Lower bound for each node
%   fIPmin    Lowest obtained integer solution objective
%   Lx        Integer variables values from latest relaxed LP solution
%              of the parent nodes to the nodes in L.
%   s0        Sum of integer variable fractions in initial relaxed solution
%   E         Best estimate of nodes in L
%
% OUTPUT:
%   Node      Selected node in L
%   L         Updated node list (Node removed)
%   Lix       Index of selected node in L

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Oct 15, 2009.   Last modified Jun 25, 2013.

Lix = [];

if length(L) == 1
   Node = L(1);
   L = [];
   return
end

if BACKTRACK
   NodeSelMethod = BacktrackMethod;
end

% Problem selection and relaxation
if NodeSelMethod == 1
   % Find nodes with max depth, deepest down in the tree
   % Select the node with most integer components, first one if ties
   Lev    = Depth(L);
   LevMax = max(Lev);
   ix     = find(Lev==LevMax);
   if length(ix) > 1
      IC  = Icomp(pred(L(ix)));
      [ICmax, ICidx] = max(IC);
      i   = L(ix(ICidx));
      Lix = ix(ICidx);
   else
      i   = L(ix);
      Lix = ix;
   end

   L = setdiff(L,i);
elseif NodeSelMethod == 2
   % Find nodes with min depth, highest up in the tree
   % Select the node with most integer components, first one if ties
   Lev    = Depth(L);
   LevMin = min(Lev);
   ix     = find(Lev==LevMin);
   if length(ix) > 1
      IC = Icomp(pred(L(ix)));
      [ICmax, ICidx] = max(IC);
      i  = L(ix(ICidx));
   else
      i  = L(ix);
   end
   % fprintf('i %d OLD FIFO i %d\n',i,L(1));
   L = setdiff(L,i);
elseif NodeSelMethod==3
   % Old LIFO depth search, last in, first out
   i = L(length(L));
   L = L(1:length(L)-1);
elseif NodeSelMethod==4
   % Old FIFO breadth search, first in, first out
   i = L(1);
   L = L(2:length(L));
elseif NodeSelMethod==5
   % Maximize integer components
   IC                 = Icomp(pred(L));
   [ICmax, Lix]     = max(IC);
   if ICmax == nI-1
      i               = L(Lix);
   else
      ix = find(IC==ICmax);
      if length(ix) == 1
         i            = L(Lix);
      else
         [Fmin, Lix] = min(f_min(L));
         i            = L(Lix);
      end
   end
   L = setdiff(L,i);
elseif NodeSelMethod==6
   % Best bound
   [Fmin,Lix] = min(f_min(L));
   i = L(Lix);
   L = setdiff(L,i);
elseif NodeSelMethod==7
   % Best estimate using pseudo-costs
   [Hn,lix] = min(E);
   i = L(lix);
   L = setdiff(L,i);
elseif NodeSelMethod==8
   % Best projection
   s = sum(Lx,1)';
   E = f_min(L) + s*((fIPMin) - pred(1))/s0;
   [notused,Lix] = min(E);
   i = L(Lix);
   L = setdiff(L,i);
end

Node = i;

% MODIFICATION LOG
%
% 091015  hkh  Written
% 130626  bjo  Removed methods with strategy switching, handled by new input
%              arguments BacktrackMethod and BACKTRACK flag instead.
% 130704  bjo  Added best bound and best projection.
% 130717  bjo  Added best estimate using pseudo-costs
