function [sample_points warmup_points state]= model_sample(model, samples_count, warmup_points, state)
    if nargin < 3
        warmup_points = model_sampling_warmup(model, 5000, 1);
    end;
    if nargin < 4
        [throw_away_points, state] = model_sampling_ACHR(model, warmup_points,  1000, 400);
    end
    [sample_points state] = model_sampling_ACHR(model, warmup_points,  samples_count, 400, state);
end
