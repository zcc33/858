function[model] = defineHumanMediaNCI60(model,flag) %RPMI-1640

amino_idxes = [
    % Amino acids
    3457 % L-Alanine
    165 % L-Arginine
    1590 % L-Asparagine
    2387 % L-Aspartate
    1089 % L-Cystine
    1372 % L-Glutamine
    2567 % L-Glutamate
    2510 % Glycine
    3719 % L-Histidine
    2655 % L-Isoleucine
    889 % L-Leucine
    2721 % L-Lysine
    2978 % L-Methionine
    997 % L-Phenylalanine
    2262 % L-Proline
    1422 % L-Serine
    2363 % L-Threonine
    1993 % L-Tryptophan
    3206 % L-Tyrosine
    2763 % L-Valine
    ]';
vitamin_idxes = [
    % Vitamins
    1156 % Biotin
    662 % Choline
    3471 % Pantothenate
    2228 % Folate
    110 % myo-Inositol
    1457 % Nicotinamide
    1676 % Pyridoxal
    3385 % Riboflavin
    420 % Thiamine
    % B12 (not in model)
    ]';
mineral_idxes = [
    % Minerals
    3165 % Calcium
    2039 % Iron
    % Magnesium (not in model)
    1210 % Potassium
    3670 % Sodium
    163 % Chloride
    1943 % Phosphate
    % Sulfate
    %zinc
    ]';
glc_idx = 968; % Glucose
o2_idx = 3446; % Oxygen
gln_idx = 1372; % Glutamine
gtrd_idx = 463; % Gluthatione reduced

if flag
    ex_rxns = strmatch('EX_',model.rxns);
    model.lb(ex_rxns) = 0;
    model.lb(amino_idxes) = -0.05;
    model.lb(vitamin_idxes) = -0.005;
    model.lb(mineral_idxes) = -5;
    model.lb(glc_idx) = -5;
    model.lb(gtrd_idx) = -0.005;
    model.lb(gln_idx) = -0.5;
    model.lb(o2_idx) = -10;
end
model.rowlb = zeros(length(model.mets),1);
model.rowub = zeros(length(model.mets),1);
model.int_vars = zeros(length(model.rxns),1);
model.c(3745) = 1;
