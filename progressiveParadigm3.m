%% THREE STIMULUS PAIRS
%This function will create a 2 x (n_trials) size matrix. 
%Row 1: vector of 0's,1's, or 2's which correspond to whether the trial is
%for sequence AB, BC, or CA, respectively. 
%Row 2: vector of 0's and 1's for whether the trial is a positive or
%negative outcome (0 == -, 1 == +)

function [data,pm] = progressiveParadigm3(n_trials, pm)
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
    if blocki < 0.3*length(blocks)
        b = 0;
    else
        b = 1;
    end
    
    %AC-pair
    if blocki < 0.6*length(blocks)
        c =0;
    else
        c = 1;
    end
    
    % AB will always appear; 
    % BC will appear after the first third blocks;
    % AC will appear after the second third blocks;
    stimuli = [stimuli a*zeros(1,blocks(blocki)) b*ones(1,blocks(blocki)) c*2*ones(1,blocks(blocki))]; 
end
num_testing_trials = pm.testing_reps_per_stim;
testing = [zeros(1,num_testing_trials) ones(1,num_testing_trials) 2*ones(1,num_testing_trials)];
testing = testing(randperm(length(testing)));
data.training_length = length(stimuli);
stimuli = [stimuli testing];

%save testing sequences for countDecision
data.testing_seq = testing;
data.stimuli = stimuli;
%% 
% p1 = round(0.27*n_trials);
% p2 = round(0.07*n_trials);
% p3 = round(0.04*n_trials);
% p4 = round(0.02*n_trials);
% p5 = round(0.01*n_trials);
%  
% %1 x 30 = 27 AB
% stimuli = [stimuli zeros(1,p1)];
% 
% %2 x 11 = 22 BC & AB
% stimuli = [stimuli ones(1,p2) zeros(1,p3)];
% stimuli = [stimuli ones(1,p2) zeros(1,p3)];
% 
% %3 x 4 = 12 BC & AB
% stimuli = [stimuli ones(1,p4) zeros(1,p4)];
% stimuli = [stimuli ones(1,p4) zeros(1,p4)];
% stimuli = [stimuli ones(1,p4) zeros(1,p4)];
% 
% %1 x 7 = 7 CA
% stimuli = [stimuli 2*ones(1,p2)];
% 
% %7 x 2 = 14 BC AB BC CA AB BC CA
% stimuli = [stimuli ones(1,p4) zeros(1,p4) ones(1,p4) 2*ones(1,p4) zeros(1,p4) ones(1,p4) 2*ones(1,p4)];
% 
% %remaining ~82 percent
% indexes_left = n_trials-(1*p1+3*p2+2*p3+13*p4);
% disp("indexes_left")
% disp(indexes_left);
% final_block = [zeros(1,indexes_left/3) ones(1, indexes_left/3) 2*ones(1, indexes_left-(2*indexes_left/3))]; % this regerates a 1 x indexes_left matrix of random ints between 0 and 2 
% disp(final_block);
% stimuli = [stimuli final_block(randperm(indexes_left))];
% disp("REAL PROGRESSIVE PARADIGM")
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
