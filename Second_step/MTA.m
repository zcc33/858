%MTA:
%Input - 1. The Model iwht a given media
%                2. v_ref -  a reference vector describing the source state as derived from an iMAT-based analysis
%                3. discrete_rxns_vector - A discrete vector (1,-1,0), describing reactions that their
%                    flux should increase/decrease (1,-1, respectively) or
%                    remain constant (0) in order to transform to the
%                    target state
%                4. rxns_to_delete - A set of reactions to KO
%Output - 1. score - The score obtained by MTA following the KO of each of
%                    the reactions in rxns_to_delete
%         2. stat -  The solver's returned status for each KO

function[score,stat] = MTA(model,v_ref,discrete_rxns_vector,rxns_to_delete)

DEFINE_PARAM;
thr = 0.01; %threshold defining the the minimum required change for the integer constraints
e_rxns = find(discrete_rxns_vector==1);
r_rxns = find(discrete_rxns_vector==-1);
s_rxns = find(discrete_rxns_vector==0);
s_rxns(find(s_rxns==(find(model.c==1)))) = [];

fwd = [intersect(find(v_ref >= 0),e_rxns)  intersect(find(v_ref < 0),r_rxns)];
bck = [intersect(find(v_ref <= 0),e_rxns)  intersect(find(v_ref > 0),r_rxns)];

%Building the constraints matrix
cons_rxns = [];
[m_s,n_s] = size(model.S);
for i=1:length(fwd)
    model = addFWDCons(fwd(i),v_ref(fwd(i)),model,thr);
    cons_rxns = [cons_rxns fwd(i)];
end


for i=1:length(bck)
    if v_ref(bck(i))==0 & model.lb(bck(i))==0
        continue;
    end
    model = addBCKCons(bck(i),v_ref(bck(i)),model,thr);
    cons_rxns = [cons_rxns bck(i)];
end


%Constructing the optimization integer vector
[m,n] = size(model.S);
model.int_vars(n_s+1:n) = 1;
model.c = zeros(n,1);
model.c([n_s+2:2:n]) = 1;

%Constructing the optimization QP vector and matrix
vec = zeros(n,1);
vec(s_rxns) = 2;
model.F = diag(vec);
model.c(s_rxns) = -2*v_ref(s_rxns);

%%%
model.c([n_s+2:2:n]) = alpha/2;%discrete_rxns_p(cons_rxns);

vec = zeros(n,1);
vec(s_rxns) = 2*(1-alpha);
model.F = diag(vec);
model.c(s_rxns) = -2*v_ref(s_rxns).*(1-alpha);
%%%

LB = model.lb;
UB = model.ub;

for i=1:length(rxns_to_delete)
    i
    model.lb = LB;
    model.ub = UB;
    
    %Performing the KO
    model.lb(rxns_to_delete(i)) = 0;
    model.ub(rxns_to_delete(i)) = 0;
    
    %Running the optimization
    Res = RunTomlabMIQP(model,0);
    
    stat(i)= Res.result_status;
    
    %Calculating the score
    [diff_change] = calculateDiff(v_ref,fwd,bck,cons_rxns,Res.result_vector);
    diff_steady = sum(abs(Res.result_vector(s_rxns)-v_ref(s_rxns)));
    score(i) = diff_change/diff_steady;
end











