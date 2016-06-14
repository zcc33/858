mfunction[model] = defineHumanMediaRPMI(model)

amino_idxes1 = 165; % L-Arginine

amino_idxes2 = [
    % Amino acids
    3457 % L-Alanine
    1590 % L-Asparagine
    2387 % L-Aspartate
    1089 % L-Cystine
    2567 % L-Glutamate
    2510 % Glycine
    3719 % L-Histidine
    2655 % L-Isoleucine
    889 % L-Leucine
    2721 % L-Lysine
    2978 % L-Methionine
    997 % L-Phenylalanine
    2269 % L-Proline
    1422 % L-Serine
    2363 % L-Threonine
    3206 % L-Tyrosine
    2763 % L-Valine
    ]';
amino_idxes3 = 1993; % L-Tryptophan

% Vitamins
vitamin_idxes1 = 110; % myo-Inositol

vitamin_idxes2 = [
    662 % Choline
    2228 % Folate
    1457 % Nicotinamide
    1676 % Pyridoxal
    420 % Thiamine
    ]';
vitamin_idxes3 = [
    1156 % Biotin
    3471 % Pantothenate
    3385 % Riboflavin
    ]';

% Minerals
mineral_idxes1 = [
    3670 % Sodium
    163 % Chloride
    ]';
mineral_idxes2 = [
    3165 % Calcium
    1210 % Potassium
    1943 % Phosphate
    3279 % H2O
    ]';

%Nucleotides
nucleotides_idx1 = [
    2913 % CMP
    2652 % UMP
    ]';
nucleotides_idx2 = [
    1753 % Adenine
    1715 % AMP
    91 % GMP
    3384 %IMP
    ]';

nucleotides_idx3 = [
    1680 %Cytidine
    956 % Deoxycytidine
    2026 % Deoxyurridine
    3075 % Thymidine
    3177 % Uracil
    877 % Uridine
    ]';

nucleotides_idx4 = 2241; % Urate

lipids = [
520 % linoleic acid
1639 % Hexadecanoate (n-C16:0) (Palmitic acid)
2805 % 'Hexadecenoate (n-C16:1) (palmitoleic acid)
]';

Others1 = 557; % Triiodothyonine

Others2 = [
    2253 % L-Thyroxine
    1553 % Thymine
    77 % Taurocholic acid
    3337 % Oxidized glutathione
    1768 % glycocholate
    1632 % 4-Pyridoxate
    3674 % 4-Aminobutanoate
    ]';

Others3 = [
    2627 % Acetoacetate
    2683 % 2-Oxoglutarate
    2598 % L-Carnitine
    1892 % Lactose
    675 % Serotonin
    1741 % Sucrose
    ]';

Others4 = [
    1803 % Bilirubin
    1496 % Citrate
    1514 % Glycerol
    1001 % Hypoxanthine
    3361 % Oxalate
    3174 % Taurine
    ]';

Others5 = [
    1841 % Creatine
    357 % L-Homoserine
    3153 % Ornithine
    2701 % Succinate
    ]';

Others6 = 941; % L-Lactate

glc_idx = 968; % Glucose
o2_idx = 3446; % Oxygen
gln_idx = 1372; % Glutamine
gtrd_idx = 463; % Reduced gluthatione


ex_rxns = strmatch('EX_',model.rxns);
model.lb(ex_rxns) = 0;
model.lb(amino_idxes1) = -1;
model.lb(amino_idxes2) = -0.1;
model.lb(amino_idxes3) = -0.01;

model.lb(vitamin_idxes1) = -0.1;
model.lb(vitamin_idxes2) = -0.01;
model.lb(vitamin_idxes3) = -0.001;

model.lb(mineral_idxes1) = -100;
model.lb(mineral_idxes2) = -1;

model.lb(nucleotides_idx1) = -1e-6;
model.lb(nucleotides_idx2) = -1e-5;
model.lb(nucleotides_idx3) = -1e-4;
model.lb(nucleotides_idx4) = -0.001;

model.lb(lipids) = -0.1;

model.lb(Others1) = -1e-6;
model.lb(Others2) = -1e-5;
model.lb(Others3) = -1e-4;
model.lb(Others4) = -0.001;
model.lb(Others5) = -0.01;
model.lb(Others6) = -0.1;

model.lb(glc_idx) = -10;
model.lb(gtrd_idx) = -0.001;
model.lb(gln_idx) = -2;
model.lb(o2_idx) = -100;
