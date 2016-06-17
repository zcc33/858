%function [] = main()

addpath 'First step'
addpath 'Second_step'
load LUAD_data.mat
load model.mat

[healthy, cancer]= getExpressions(model, Data);
[~, testing] = sampleiMAT(model, healthy(1:end, 1))