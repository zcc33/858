function [rxns] = getRxnsToDelete(model, source, target)
%for use in MTA. get the set of reactions that are relevant, that are
%candidates for deletion in the MTA run. 

%Returns: a list of indices of relevant reactions
%NEED TO INPUT THE NON-DISCRETIZED GENE EXPRESSIONS

rxns1 = zeros(size(model.rxns));
for i = 1:length(model.rxns)
    rxns1(i) = evalExpRule2(model.rules{i}, source) | evalExpRule2(model.rules{i}, target);
end

rxns = zeros(sum(rxns1), 1);
it = 1;
for i = 1:length(model.rxns)
    if rxns1(i) == 1
        rxns(it) = i;
        it = it+1;
    end
end