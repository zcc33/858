function [] = initial()
%run this once to do the initialization work

%add the relevant folder paths for all the sub-files we may need
addpath 'First step';
addpath 'Second_step';

%load the model, and then run it through the DMEM media
load('model.mat');
model = defineHumanMediaDMEM(model, 1);

%load the lung cancer data
load('LUAD_data.mat');

%get the healthy/cancer (and their discretized) samples we want
[master_healthy, master_cancer, master_healthy_d, master_cancer_d] = getExpressions(model, Data);

%how many relevant samples we got
[~, num_samples] = size(master_healthy);

%create the data directory, a folder for each sample, and a log file for
%all the samples telling us if it's not done (0) or done (1) processing
mkdir('data')
log_progress = zeros(num_samples,1);
save('data/log_progress.mat', 'log_progress')
for i = 1:num_samples
    mkdir('data',int2str(i));
end

%now save the sample data in the data file for the main function to use
save('data/master_healthy.mat', 'master_healthy')
save('data/master_cancer.mat', 'master_cancer')
save('data/master_healthy_d.mat', 'master_healthy_d')
save('data/master_cancer_d.mat', 'master_cancer_d')