% Colton Bogucki
% Purpose: DetermineTrialSuccessTP is a success detector for transverse
% patterning sequences. A success is determined by comparing the number of
% repeated fires for neurons within the correct and incorrect decision
% block. If there are more repeated fires in the correct decison block, the
% trial was a success. Otherwise, the trial was a failure.
% Returns: result - 'success' or 'failure'
%          b1,b2 - correct and incorrect decision blocks for visualization
%          within ScrollThroughTP
%          right_decision,wrong_decision - number of repeated fires within
%          the correct and incorrect decision blocks

%[block1,block2,right_decision,wrong_decision,
function [block1,block2,right_decision,wrong_decision,result] = DetermineTrialSuccessTP(parameters,data,trial,stim)
%constants
stim_nrns           = parameters.stim_nrns; %number of stimulus neurons
dec_nrns            = parameters.dec_nrns; %number of decision neurons
stim_dwell          = parameters.stim_dwell; %number of timesteps stimulus neurons are on
dec_dwell           = parameters.dec_dwell; %number of timesteps decision neurons are on
outcome_dwell       = parameters.outcome_dwell; %number of timesteps correctness neurons are on
% current trial test matrix
test_mat = data.z_test(:,:,trial);

%AB
if stim == 0
    %'a' block
    start1 = 3*stim_nrns;
    finish1 = start1+dec_nrns;
    %'b' block
    start2 = finish1;
    finish2 = start2+dec_nrns;
%BC
elseif stim == 1
    %'b' block
    start1 =  3*stim_nrns+dec_nrns;
    finish1 = start1 + dec_nrns;
    %'c' block
    start2 = finish1;
    finish2 = start2+dec_nrns;
%CA
else
    %'c' block
    start1 = 3*stim_nrns + 2*dec_nrns;
    finish1 = start1 + dec_nrns;
    %'a' block
    start2 = 3*stim_nrns;
    finish2 = start2+dec_nrns;
end
pre_train = parameters.pre_train_length;
train_len = parameters.length_of_each_trial;
%matrix of the right decision
block1 = test_mat(start1:finish1,pre_train+1:train_len);
%matrix of the wrong decision
block2 = test_mat(start2:finish2,pre_train+1:train_len);

% %repeated fires for block 1 (the correct stimulus decision)
% [nrns,time] = size(block1);
% fires = zeros(nrns,1);
% for i = 1:time-1
%     %sum row of neurons at timestep i and i+1, if the sum is 2, then the
%     %neuron fired repeatedly
%     summed = (block1(:,i)+block1(:,i+1) == 2);
%     fires = fires+summed;
% end
% right_decision = sum(fires>0);
%  
% %repeated fires for block 2 (the incorrect stimulus decision)
% [nrns,time] = size(block2);
% fires = zeros(nrns,1);
% for i = 1:time-1
%     %sum row of neurons at timestep i and i+1, if the sum is 2, then the
%     %neuron fired repeatedly
%     summed = (block2(:,i)+block2(:,i+1) == 2);
%     fires = fires+summed;
% end
% wrong_decision = sum(fires>0);
right_decision = sum(sum(block1));
wrong_decision = sum(sum(block2));

if right_decision > wrong_decision
    result = 'success';
    %disp(['Trial ',num2str(trial),': ',result,', ',num2str(right_decision),' > ',num2str(wrong_decision)])
else
    result = 'failure';
    %disp(['Trial ',num2str(trial),': ',result,', ',num2str(right_decision),' =< ',num2str(wrong_decision)])
end

end