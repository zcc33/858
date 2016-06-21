function [healthy, cancer, healthy_discrete, cancer_discrete] = getExpressions(model, Data)

%which are the relevant samples we want

%this empty matrix will come to store all the column indices of the healthy
%samples we want
[~, m] = size(Data.GE);
relevant_healthy=zeros(1, m);
for i = 1:m-1
    
    %this is to make sure that the healthy samples are followed by a cancer
    %sample, and moreover, that this cancer sample corresponds to the ID of
    %the healthy sample
    
    %b returns the string that tells us whether the index is for healthy or
    %cancer
    b=Data.sample{i};
    b2=Data.sample{i+1};
    
    if strcmp(b(end-1:end),'11') && strcmp(b2(end-1:end), '01') && strcmp(b(end-6:end-3), b2(end-6:end-3))
        relevant_healthy(i) = 1;
    else
        relevant_healthy(i) = 0;
    
    end
end

%get the indices of relevant_healthy in the Data.GE
indices = zeros(1, sum(relevant_healthy));
iterator = 1;
for j = 1:m-1
    if relevant_healthy(j) ==1
        indices(iterator) = j;
        iterator = iterator +1;
    end
end

%get their gene expression data
%supposing the cancer ones immediately follow the healthy samples
num_unique_genes = length(model.genes_unique);
num_samples = length(indices);

healthy = zeros(num_unique_genes, num_samples);
cancer = zeros(num_unique_genes, num_samples);

for i = 1:num_unique_genes
    a = find(Data.gene_name == model.genes_unique(i));
    if ~isempty(a)
        for j = 1:num_samples
            healthy(i, j) = Data.GE(a, indices(j));
            cancer(i, j) = Data.GE(a, indices(j)+1);
        end
    end
end




%map from unique genes to non-unique
num_nonunique_genes = length(model.genes_unique_map);
[n2, m2] = size(healthy);

healthy1 = zeros(num_nonunique_genes, m2);
cancer1 = zeros(num_nonunique_genes, m2);

it = 1;
for i = 1:n2
    a = find(model.genes_unique_map == i);
    for j = 1:length(a)
        healthy1(it, 1:end) = healthy(i, 1:end);
        cancer1(it, 1:end) = cancer(i, 1:end);
        it = it+1;
    end
end

%these are the non-discrete versions that will be returned
healthy = healthy1;
cancer = cancer1;


%time to discretize, which will also be returned
[n, m] = size(cancer);
cancer_discrete = zeros(n,m);
for i = 1:m
    twentyfifth = prctile(cancer(1:end, i), 25);
    seventyfifth = prctile(cancer(1:end, i), 75);
    for j = 1:n
        if cancer(j, i) < twentyfifth
            cancer_discrete(j, i) = -1;
        elseif cancer(j,i) < seventyfifth
            cancer_discrete(j,i) = 0;
        else
            cancer_discrete(j,i) = 1;
        end
    end
end


[n, m] = size(healthy);
healthy_discrete = zeros(n,m);
for i = 1:m
    twentyfifth = prctile(healthy(1:end, i), 25);
    seventyfifth = prctile(healthy(1:end, i), 75);
    for j = 1:n
        if healthy(j, i) < twentyfifth
            healthy_discrete(j, i) = -1;
        elseif healthy(j,i) < seventyfifth
            healthy_discrete(j,i) = 0;
        else
            healthy_discrete(j,i) = 1;
        end
    end
end