%% THREE STIMULUS PAIRS
%This function will create a 2 x (n_trials) size matrix. 
%Row 1: vector of 0's,1's, or 2's which correspond to whether the trial is
%for sequence AB, BC, or CA, respectively. 
%Row 2: vector of 0's and 1's for whether the trial is a positive or
%negative outcome (0 == -, 1 == +)

function [data,pm] = progressiveParadigm(pm, data)
%%
n_trials = pm.number_of_trials;
blocks = pm.block_reps;
padding = pm.padding;
stimuli = [zeros(1, padding)];
a = 1;
BC = 0;
AC = 0;
phase1 = 100;
phase2 = 100;
if pm.pairs_to_learn == 2
    BC = 1;
    phase1 = 2;
elseif pm.pairs_to_learn == 3
    BC = 1;
    AC = 1;
    phase2 = 4;
end

for blocki = 1:length(blocks)
    %a, b, c are used to control in which block the stimuli pair AB, BC, AC
    %appears. They are used on line
        
    %BC-pair
    if (blocki > phase1) && BC
        b = 1;
    else
        b = 0;
    end
    
    %AC-pair
    if (blocki > phase2) && AC
        c = 1;
    else
        c = 0;
    end
    disp("0.6*length(blocks)")
    disp(0.6*length(blocks))
    disp("blocki")
    disp(blocki)
    disp("B")
    disp(b)
    disp("C")
    disp(c)
    
    
    % AB will always appear; 
    % BC will appear after ~the first third blocks;
    % AC will appear after ~the second third blocks;
    stimuli = [stimuli a*zeros(1,blocks(blocki)) ones(1,b*blocks(blocki)) c*2*ones(1,c*blocks(blocki))]; 
end
num_testing_trials = pm.testing_reps_per_stim;
testing = [zeros(1,num_testing_trials) ones(1,b*num_testing_trials) 2*ones(1,c*num_testing_trials)];
testing = testing(randperm(length(testing)));
data.training_length = length(stimuli);
stimuli = [stimuli testing];

%save testing sequences for countDecision
data.testing_seq = testing;
data.stimuli = stimuli;
%%
%Approximately 85% positive outcomes in sequences
outcome = rand(1,n_trials)<0.85;
disp(length(testing));
disp("actual length of stimuli: ");
disp(length(stimuli));
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

title([num2str(round((blocks(1)/n_trials)*100)),'% Initial Block'])
hold off

end
