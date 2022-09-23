%% ONE STIMULUS PAIR
%This function will create a 2 x (n_trials) size matrix. 
%Row 1: vector of 0's correspond to AB stimulus pair. 
%Row 2: vector of 0's and 1's for whether the trial is a positive or
%negative outcome (0 == -, 1 == +)

function training_sequence = stagedParadigm1(n_trials)
stimuli = zeros(1,n_trials);


%%
%Approximately 85% positive outcomes in sequences
outcome = rand(1,n_trials)<0.85;

%Concatenation of two row vectors into 2xn_trials matrix describing staged
%paradigm
%stimuli = zeros(1,n_trials);
length(stimuli);
length(outcome);
training_sequence = [stimuli;outcome];
% figure
% hold on
% for i = 1:n_trials
%     if(outcome(i)==0)
%         rectangle('Position',[(i-1) 0 1 1],'FaceColor','black');
%     else
%         rectangle('Position',[(i-1) 0 1 1],'FaceColor','red');
%     end
% end
% yticks([0.5]);
% yticklabels({'AB'})
% xlabel('Trial Number')
% ylabel('Stimulus Pair')
% 
% title('Red = +, Black = -')
% hold off


end
