function [healthy1, cancer1] = getExpressions(model, Data)

%which are the relevant samples we want
relevant_healthy=[];
[~, m] = size(Data.GE);
for i = 1:m-1
    
    %this is to make sure that the healthy samples are followed by a cancer
    %sample, and moreover, that this cancer sample corresponds to the ID of
    %the healthy sample
    a=Data.sample(i);
    b=a{1};
    
    a2=Data.sample(i+1);
    b2=a2{1};
    
    if strcmp(b(end-1:end),'11') & strcmp(b2(end-1:end), '01') & strcmp(b(end-6:end-3), b2(end-6:end-3))
        relevant_healthy(i) = 1;
    else
        relevant_healthy(i) = 0;
    
    end
end

%get their gene expression data

healthy = zeros(length(model.genes_unique), sum(relevant_healthy));
cancer = zeros(length(model.genes_unique), sum(relevant_healthy));
[n, ~] = size(model.genes_unique);

it = 1;
for j = 1:m-1
    
    if relevant_healthy(j) == 1
        for i = 1:n
            a = find(Data.gene_name == model.genes_unique(i));
            if length(a) > 0
                healthy(i, it) = Data.GE(a, j);
            end
        end
        it = it+1;
    end
end

%supposing the cancer ones immediately follow the healthy samples
it = 1;
for j = 1:m-1
    
    if relevant_healthy(j) == 1
        for i = 1:n
            a = find(Data.gene_name == model.genes_unique(i));
            if length(a) > 0
                cancer(i, it) = Data.GE(a, j+1);
            end
        end
        it = it+1;
    end
end



%time to discretize
[n m] = size(cancer);
for i = 1:m
    twentyfifth = prctile(cancer(1:end, i), 25);
    seventyfifth = prctile(cancer(1:end, i), 75);
    for j = 1:n
        if cancer(j, i) < twentyfifth
            cancer(j, i) = -1;
        elseif cancer(j,i) < seventyfifth
            cancer(j,i) = 0;
        else
            cancer(j,i) = 1;
        end
    end
end


[n m] = size(healthy);
for i = 1:m
    twentyfifth = prctile(healthy(1:end, i), 25);
    seventyfifth = prctile(healthy(1:end, i), 75);
    for j = 1:n
        if healthy(j, i) < twentyfifth
            healthy(j, i) = -1;
        elseif healthy(j,i) < seventyfifth
            healthy(j,i) = 0;
        else
            healthy(j,i) = 1;
        end
    end
end


%map from unique genes to non-unique
[n m] = size(model.genes_unique_map);
[n2 m2] = size(healthy);

healthy1 = zeros(n, m2);
cancer1 = zeros(n, m2);

it = 1;
for i = 1:n2
    a = find(model.genes_unique_map == i);
    for j = 1:length(a)
        healthy1(it, 1:end) = healthy(i, 1:end);
        it = it+1;
    end
end

it = 1;
for i = 1:n2
    a = find(model.genes_unique_map == i);
    for j = 1:length(a)
        cancer1(it, 1:end) = cancer(i, 1:end);
        it = it+1;
    end
end