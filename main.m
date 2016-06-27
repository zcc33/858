function [] = main()
%remember, we are going from cancer (source) to healthy (target)

%add the relevant folder paths for all the sub-files we may need
addpath 'First step';
addpath 'Second_step';

%load the model, and then run it through the DMEM media
load('model.mat');
model = defineHumanMediaDMEM(model, 1);

%load the lung cancer data
load('LUAD_data.mat');

%load the sample healthy/cancer data from the initialization
load('data/master_healthy.mat');
load('data/master_cancer.mat');
load('data/master_healthy_d.mat');
load('data/master_cancer_d.mat');

%load the log file to see which samples already completed
load('data/log_progress.mat');

%how many samples we need to go through
num_samples = length(log_progress);

%go through the log and process each sample that hasn't already been
%processed
i = 1;
num_errors = 0;

while i < (num_samples+1)
    
    %if the sample is completed, say so
    if log_progress(i,1) ==1
        msg = ['Sample ' int2str(i) ' (out of ' int2str(num_samples) ') completed'];
        disp(msg);
        
        %continue iteration
        i = i+1;
    
    %if the sample was too error-riddled to continue, say so
    elseif log_progress(i,1) ==2
        msg = ['Sample ' int2str(i) ' (out of ' int2str(num_samples) ') had too many ERRORS.'];
        disp(msg);
        
        %continue iteration
        i = i+1;
        break;
        
    
    %if the sample isn't completed, process it, save all data, then say so
    elseif log_progress(i,1) == 0
        
        %get the sample-specific data from the master data files
        healthy = master_healthy(:,i);
        cancer = master_cancer(:,i);
        healthy_d = master_healthy_d(:,i);
        cancer_d = master_cancer_d(:,i);
        
        %save those sample-specific data to the sample's folder
        save(['data/' int2str(i) '/healthy.mat'], 'healthy')
        save(['data/' int2str(i) '/cancer.mat'], 'cancer')
        save(['data/' int2str(i) '/healthy_d.mat'], 'healthy_d')
        save(['data/' int2str(i) '/cancer_d.mat'], 'cancer_d')
        
        %(1) get discrete reactions vector
        discrete_rxns_vector = getDiscreteRxns(model, cancer, healthy);
        save(['data/' int2str(i) '/discrete_rxns_vector.mat'], 'discrete_rxns_vector')
        
        %(2) get reactions to delete vector
        rxns_to_delete = getRxnsToDelete(model, cancer, healthy);
        save(['data/' int2str(i) '/rxns_to_delete.mat'], 'rxns_to_delete')
        
        %print message that we are now beginning iMAT
        msg = ['Sample ' int2str(i) ' has gotten (rxns to delete) and (discrete rxns). Now beginning iMAT...'];
        disp(msg);
        
        %(3) run iMAT to get v_ref
        %remember that the input is the discrete version of the SOURCE
        %(cancer) gene expressions
        [~, v_ref] = sampleiMAT(model, cancer_d);
        save(['data/' int2str(i) '/v_ref.mat'], 'v_ref');
        
        %print message that we are done iMAT, beginning MTA
        msg = ['Sample ' int2str(i) ' has gotten (v_ref). Now beginning MTA...'];
        disp(msg);
        
        %print how many iterations this sample has undergone due to errors
        msg = ['There have been ' int2str(num_errors) ' errors with this sample. Hoping this works...'];
        disp(msg);
        
        %THERE MAY BE ERRORS HERE
        %"try" to do this if no errors
        try
            %(4) run MTA to get the score and stat vectors
            [score, stat, score_placebo, stat_placebo] = MTA(model, v_ref, discrete_rxns_vector, rxns_to_delete);
        
            %transpose to get row vectors and then save data
            score = score';
            stat = stat';
            save(['data/' int2str(i) '/score.mat'], 'score');
            save(['data/' int2str(i) '/stat.mat'], 'stat');
            save(['data/' int2str(i) '/score_placebo.mat'], 'score_placebo');
            save(['data/' int2str(i) '/stat_placebo.mat'], 'stat_placebo');
        
            %now that the sample is done, log it and save log file
            log_progress(i,1) = 1;
            save('data/log_progress.mat', 'log_progress');
        
            %print that this sample is done
            msg = ['Sample ' int2str(i) ' (out of ' int2str(num_samples) ') completed'];
            disp(msg);
            
            %continue iteration
            i = i+1;
            
            %error tracker back to 0
            num_errors = 0;
            
        %if there are errors, then do this
        catch
            %no iterating, go back and do this sample AGAIN (hopefully no
            %errors next time, since iMAT is non-deterministic
            num_errors = num_errors +1;
            
            if num_errors < 4
                
                %print that this sample is not yet done
                msg = ['Sample ' int2str(i) ' encountered error #' int2str(num_errors) '. We will retry...'];
                disp(msg);
                
            else
                
                %BUT if there have been 3 or more errors, just continue and
                %label that sample as having errors
                log_progress(i,1) = 2;
                save('data/log_progress.mat', 'log_progress');
        
                %print that this sample is done
                msg = ['Sample ' int2str(i) ' has too many errors. Moving on...'];
                disp(msg);
                
                %continue iteration
                i = i+1;
                %reset error tracker
                num_errors = 0;
            end
        end
    end
end