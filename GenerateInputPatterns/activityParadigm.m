%% THREE STIMULUS PAIRS
%This function will create a 2 x (n_trials) size matrix. 
%Row 1: vector of 0's,1's, or 2's which correspond to whether the trial is
%for sequence AB, BC, or CA, respectively. 
%Row 2: vector of 0's and 1's for whether the trial is a positive or
%negative outcome (0 == -, 1 == +)

function [stabilize] = activityParadigm(data, n_trials)
%%
stimuli = [zeros(1, n_trials)];

%Approximately 85% positive outcomes in sequences
outcome = rand(1,n_trials)<0.85;

%Concatenation of two row vectors into 2xn_trials matrix describing staged
%paradigm
%length(stimuli)
%length(outcome)
stabilize = [stimuli;outcome];
end