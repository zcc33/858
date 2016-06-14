function [getReactions]
geneIndex=[];
expressions=[];
[n,m] = size(model.genes_unique);
for i = 1:n
    a = find(Data.gene_name == model.genes_unique(i));
    if length(a) > 0
        geneIndex(i) = find(Data.gene_name == model.genes_unique(i));
        
    end
end