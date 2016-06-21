function[sample_points,v_ref] = sampleiMAT(model,discrete_data)
Private_DefineParameters;

[milp,v_act_fwr,v_act_bck,v_inact] = formiMATMilp(model, discrete_data);

rxns = [v_act_fwr; v_act_bck; v_inact];
integers = milp.start_sol(length(model.rxns)+1:end);

int_act_fwr = integers(1:length(v_act_fwr));
int_act_bck = integers(length(v_act_fwr)+1:length(v_act_fwr)+length(v_act_bck));
int_inact = integers(length(v_act_fwr)+length(v_act_bck)+1:end);

active_fwr  = v_act_fwr(find(int_act_fwr==1));
active_bck  = v_act_bck(find(int_act_bck==1));
inactive  = v_inact(find(int_inact==1));
inactive_rev = find(model.lb(inactive)<0);
inactive_non_rev = find(model.lb(inactive)==0);

model.lb(active_fwr) = ACTIVE_FLUX;
model.ub(active_bck) = -ACTIVE_FLUX;
model.lb(inactive(inactive_rev)) = -INACTIVE_FLUX;
model.lb(inactive(inactive_non_rev)) = 0;
model.ub(inactive) = INACTIVE_FLUX;

[sample_points]= model_sample(model, 1000);
v_ref = mean(sample_points,2);