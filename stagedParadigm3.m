%% THREE STIMULUS PAIRS
%This function will create a 2 x (n_trials) size matrix. 
%Row 1: vector of 0's,1's, or 2's which correspond to whether the trial is
%for sequence AB, BC, or CA, respectively. 
%Row 2: vector of 0's and 1's for whether the trial is a positive or
%negative outcome (0 == -, 1 == +)

function training_sequence = stagedParadigm3(pm, n_trials)
stimuli = [];


%The type of sequence generated on each trial is broken into 5 stages:
%Stage 1: 
%Stage 2:
%Stage 3: 
%Stage 4: 
%Stage 5: 
%Calculating percentages of n_trials that are used multiple times.

blocks = pm.block_reps;
padding = pm.padding;
stimuli = [zeros(1, padding)];
for blocki = 1:length(blocks)
    stimuli = [stimuli zeros(1,blocks(blocki)) ones(1,blocks(blocki)) 2*ones(1,blocks(blocki))];
end
%disp(size(stimuli))
num_testing_trials = pm.testing_reps_per_stim;
testing = [zeros(1,num_testing_trials) ones(1,num_testing_trials) 2*ones(1,num_testing_trials)];

stimuli = [stimuli testing(randperm(length(testing)))];
%% First 30 percent = [15% AB; 15% BC];
% p1 = round(0.15*n_trials);
% %2 x 15
% stimuli = [stimuli zeros(1,p1) ones(1,p1)];
% %2 x 10
% stimuli = [stimuli zeros(1,p2) ones(1,p2)];
% %2 x 5
% stimuli = [stimuli zeros(1,p3) ones(1,p3)];
% %6 x 2.5
% stimuli = [stimuli zeros(1,p4) ones(1,p4)];
% stimuli = [stimuli zeros(1,p4) ones(1,p4)];
% stimuli = [stimuli zeros(1,p4) ones(1,p4)];
% %4 x 1.25 =
% stimuli = [stimuli zeros(1,p5) ones(1,p5)];
% stimuli = [stimuli zeros(1,p5) ones(1,p5)];
% %4 x 0.75 = 
% stimuli = [stimuli zeros(1,p6) ones(1,p6)];
% stimuli = [stimuli zeros(1,p6) ones(1,p6)];
% %remaining ~17 percent
% indexes_left = n_trials-(2*p1+2*p2+2*p3+6*p4+4*p5+4*p6);
% final_block = floor(2*rand(1,indexes_left));
% stimuli = [stimuli final_block(randperm(indexes_left))];

%% First 25 percent = [12.5% AB; 12.5% BC];
% p1 = round(0.1250*n_trials);
% %2 x 12.5
% stimuli = [stimuli zeros(1,p1) ones(1,p1)];
% %2 x 10
% stimuli = [stimuli zeros(1,p2) ones(1,p2)];
% %4 x 5
% stimuli = [stimuli zeros(1,p3) ones(1,p3)];
% stimuli = [stimuli zeros(1,p3) ones(1,p3)];
% %4 x 2.5
% stimuli = [stimuli zeros(1,p4) ones(1,p4)];
% stimuli = [stimuli zeros(1,p4) ones(1,p4)];
% %4 x 1.25
% stimuli = [stimuli zeros(1,p5) ones(1,p5)];
% stimuli = [stimuli zeros(1,p5) ones(1,p5)];
% %4 x 0.75
% stimuli = [stimuli zeros(1,p6) ones(1,p6)];
% stimuli = [stimuli zeros(1,p6) ones(1,p6)];
% %remaining ~17 percent
% indexes_left = n_trials-(2*p1+2*p2+4*p3+4*p4+4*p5+4*p6);
% final_block = floor(2*rand(1,indexes_left));
% stimuli = [stimuli final_block(randperm(indexes_left))];
%% First 20 percent = [10% AB; 10% BC];
% %p2 = 10%
% p1 = p2;
% %2 x 10
% stimuli = [stimuli zeros(1,p2) ones(1,p2)];
% %8 x 5
% stimuli = [stimuli zeros(1,p3) ones(1,p3)];
% stimuli = [stimuli zeros(1,p3) ones(1,p3)];
% stimuli = [stimuli zeros(1,p3) ones(1,p3)];
% stimuli = [stimuli zeros(1,p3) ones(1,p3)];
% %6 x 2.5
% stimuli = [stimuli zeros(1,p4) ones(1,p4)];
% stimuli = [stimuli zeros(1,p4) ones(1,p4)];
% stimuli = [stimuli zeros(1,p4) ones(1,p4)];
% %4 x 1.25
% stimuli = [stimuli zeros(1,p5) ones(1,p5)];
% stimuli = [stimuli zeros(1,p5) ones(1,p5)];
% %4 x 0.75
% stimuli = [stimuli zeros(1,p6) ones(1,p6)];
% stimuli = [stimuli zeros(1,p6) ones(1,p6)];
% %remaining ~17 percent
% indexes_left = n_trials-(0*p1+2*p2+8*p3+6*p4+4*p5+4*p6);
% final_block = floor(2*rand(1,indexes_left));
% stimuli = [stimuli final_block(randperm(indexes_left))];
%%
%Approximately 85% positive outcomes in sequences
outcome = rand(1,n_trials)<0.85;

%stacking of two row vectors into 2xn_trials matrix describing staged
%paradigm
%stimuli = zeros(1,n_trials);
length(stimuli)
length(outcome)
training_sequence = [stimuli;outcome];
%disp(training_sequence);
figure
hold on
for i = 1:n_trials
    rectangle('Position',[(i-1) stimuli(i) 1 1],'FaceColor','black');
end
yticks([0.5 1.5, 2.5]);
yticklabels({'AB','BC','AC'})
xlabel('Trial Number')
ylabel('Stimulus Pair')

title([num2str(round((blocks(1)/n_trials)*100)),'% Initial Block'])
if pm.toggle_post_then_pre
    rules = "TT";
else rules = "TF";
end
figurename = strcat("/Volumes/Phyl/MATLAB/Levy Lab/Data/",rules,"/StagedFigure_",datestr(now, 'mm-dd-yy_HH:MM'),"_Boxcar",num2str(pm.boxcar_window),".fig");
saveas(gcf, figurename)
hold off
end
