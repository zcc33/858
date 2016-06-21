function [discrete_rxns_vector]=getDiscreteRxns(model, source, target)
%Input: the non-discretized source/target gene expression levels

%output will be -1 if decrease, 0 if same, 1 if increase
%change the prc (quantile) until get around 75-150 non-zero output entries
%uses a binary search for the prc, assuming that the higher the prc, the
%fewer non-zero entries we have, and the lower the prc, the more
 
%what i want to get is around 100 non-zero output entries

lower = 0;
upper = 1;
prc = 0.98;


discrete_rxns_vector = zeros(size(model.rxnNames));
non_zero = sum(~(discrete_rxns_vector ==0));

diff = abs(source - target);

while non_zero < 75 || non_zero > 150
    thr = quantile(diff , prc);
    
    geneState = zeros(size(source));
    geneState(source - target>thr) = -1;
    geneState(target - source>thr) = 1;
    
    discrete_rxns_vector = zeros(size(model.rxnNames));
    for i = 1:length(model.rules)
        discrete_rxns_vector(i) = evalExpRule(model.rules{i}, geneState);
    end
    
    non_zero = sum(~(discrete_rxns_vector ==0));
    
    if non_zero <75
        upper = prc;
        prc = lower + (upper-lower)/2;
    elseif non_zero > 150
        lower = prc;
        prc = lower + (upper-lower)/2;
    end
end

   
    
