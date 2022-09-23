%% THREE STIMULUS PAIRS
%This function will create a 2 x (n_trials) size matrix. 
%Row 1: vector of 0's,1's, or 2's which correspond to whether the trial is
%for sequence AB, BC, or CA, respectively. 
%Row 2: vector of 0's and 1's for whether the trial is a positive or
%negative outcome (0 == -, 1 == +)

function [data] = progressiveParadigm1(n_trials, pm, data)
%%
blocks = pm.block_reps;
padding = pm.padding;
stimuli = [zeros(1, padding)];
for blocki = 1:length(blocks)
    %a, b, c are used to control in which block the stimuli pair AB, BC, AC
    %appears. They are used on line
    
    %AB-pair
    a = 1; 
    
    %BC-pair
    b = 0;
    
    c = 0;
    
    % AB will always appear; 
    % BC will appear after the first third blocks;
    % AC will appear after the second third blocks;
    stimuli = [stimuli a*zeros(1,blocks(blocki)) b*ones(1,blocks(blocki)) c*2*ones(1,blocks(blocki))]; 
end
num_testing_trials = pm.testing_reps_per_stim;
testing = [a*zeros(1,num_testing_trials)];
testing = testing(randperm(length(testing)));
data.training_length = length(stimuli);
disp(data.training_length)
stimuli = [stimuli testing];

%save testing sequences for countDecision
data.testing_seq = testing;
data.stimuli = stimuli;
%%
%Approximately 85% positive outcomes in sequences
outcome = rand(1,n_trials)<0.85;

%Concatenation of two row vectors into 2xn_trials matrix describing staged
%paradigm
%stimuli = zeros(1,n_trials);
length(stimuli)
length(outcome)
training_sequence = [stimuli;outcome];
data.staged = training_sequence;
figure
hold on
for i = 1:n_trials
    rectangle('Position',[(i-1) stimuli(i) 1 1],'FaceColor','black');
end
yticks([0.5 1.5]);
yticklabels({'AB','BC', 'AC'})
xlabel('Trial Number')
ylabel('Stimulus Pair')

title("Progressive Paradigm")
hold off
end
