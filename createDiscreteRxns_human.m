function [discrete_rxns_vector]=createDiscreteRxns_human(source, target, smean, tmean, model, stat)
%make stat 0
%output will be -1 if decrease, 0 if same, 1 if increase
%change the prc (quantile) until get around 75-150 non-zero output entries

 %%source = smean, target=tmean
 
pval = 0.00000000005
prc = 0.9805
geneState = zeros(size(model.genes_unique));

if stat ==1
    ignoreG = zeros(size(model.genes_unique));

%     for i = 1:length(model.genes_unique)
%         if length(find(model.rxnGeneMat(:,i)))>20
%             ignoreG(i)=1;
%         end
%     end

    for i = 1:length(model.genes_unique)
        [h1(i),p1(i)] = ttest(target(i,:),smean(i),0.05,'left');
    end

    for i = 1:length(model.genes_unique)
        [h2(i),p2(i)] = ttest(target(i,:),smean(i),0.05,'right');
    end
    
    discrete_rxns_vector = zeros(size(model.rxns));

    for i = 1:length(model.genes_unique)
        if p2(i)<=pval
            geneState(i) = 1;
        end
        if p1(i)<=pval
            geneState(i) = -1;
        end
    end

end

if stat ==0

    diff = abs(source - target);
    thr = quantile(diff , prc);
    
    geneState((source - target)>thr) = -1;
    geneState((target - source)>thr) = 1;
    
end
    




geneState2 = zeros(size(model.genes));

for i=1:length(model.genes_unique)
    ind = model.genes_unique_map==i;
    geneState2(ind) = geneState(i);
end
    
for i = 1:length(model.rules)
    discrete_rxns_vector(i) = evalExpRule(model.rules{i}, geneState2);
end

discrete_rxns_vector=discrete_rxns_vector.';

