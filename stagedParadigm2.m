%% TWO STIMULUS PAIRS
%This function will create a 2 x (n_trials) size matrix. 
%Row 1: vector of 0's,1's, or 2's which correspond to whether the trial is
%for sequence AB, BC, or CA, respectively. 
%Row 2: vector of 0's and 1's for whether the trial is a positive or
%negative outcome (0 == -, 1 == +)

function training_sequence = stagedParadigm2(n_trials)
stimuli = [];


%The type of sequence generated on each trial is broken into 5 stages:
%Stage 1: 
%Stage 2:
%Stage 3: 
%Stage 4: 
%Stage 5: 
%Calculating percentages of n_trials that are used multiple times.
p1 = round(0.15*n_trials);
p2 = round(0.1*n_trials);
p3 = round(0.05*n_trials);
p4 = round(0.025*n_trials);
p5 = round(0.0125*n_trials);
p6 = round(0.0075*n_trials);

%% First 40 percent = [20% AB; 20% BC];
p1 = round(0.2*n_trials);
%2 x 20
stimuli = [stimuli zeros(1,p1) ones(1,p1)];
%2 x 10
stimuli = [stimuli zeros(1,p2) ones(1,p2)];
%2 x 5
stimuli = [stimuli zeros(1,p3) ones(1,p3)];
%2 x 2.5
stimuli = [stimuli zeros(1,p4) ones(1,p4)];
%4 x 1.25
stimuli = [stimuli zeros(1,p5) ones(1,p5)];
stimuli = [stimuli zeros(1,p5) ones(1,p5)];
%4 x 0.75
stimuli = [stimuli zeros(1,p6) ones(1,p6)];
stimuli = [stimuli zeros(1,p6) ones(1,p6)];
%remaining ~17 percent
indexes_left = n_trials-(2*p1+2*p2+2*p3+2*p4+4*p5+4*p6);
final_block = floor(2*rand(1,indexes_left));
stimuli = [stimuli final_block(randperm(indexes_left))];
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
% %4 x 1.25
% stimuli = [stimuli zeros(1,p5) ones(1,p5)];
% stimuli = [stimuli zeros(1,p5) ones(1,p5)];
% %4 x 0.75
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

%Concatenation of two row vectors into 2xn_trials matrix describing staged
%paradigm
%stimuli = zeros(1,n_trials);
length(stimuli);
length(outcome);
training_sequence = [stimuli;outcome];
% figure(1)
% hold on
% for i = 1:n_trials
%     rectangle('Position',[(i-1) stimuli(i) 1 1],'FaceColor','black');
% end
% yticks([0.5 1.5]);
% yticklabels({'AB','BC'})
% xlabel('Trial Number')
% ylabel('Stimulus Pair')
% title([num2str(round((p1/n_trials)*100)),'% Initial Block'])
% hold off



end
