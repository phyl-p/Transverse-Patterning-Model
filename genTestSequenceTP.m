function test_sequence = genTestSequenceTP(pm,stimulus)
%This function creates the testing sequence for transverse patterning,
%specifically using the method of induced attractors.

    % PARAMETERS 
    %-------------------
    % Parameters from InitialiseParameters.m:  
    n_nrns          = pm.number_of_neurons;
    timescale       = pm.boxcar_window; %for boxcar adjustments
    %Must change test_length_of_each_trial in InitialiseParameters.m
    len_trial       = pm.test_length_of_each_trial;
    
    % Parameters to be set by user:
    stim_nrns           = pm.stim_nrns; %number of stimulus neurons
    dec_nrns            = pm.dec_nrns; %number of decision neurons
    outcome_nrns        = pm.outcome_nrns; %number of correctness neurons
    stim_dwell          = pm.stim_dwell; %number of timesteps stimulus neurons are on
    dec_dwell           = pm.dec_dwell; %number of timesteps decision neurons are on
    outcome_dwell       = pm.outcome_dwell; %number of timesteps correctness neurons are on
    extra_steps         = 0; %padding zeros at the end of training sequence
    frac_outcome_test = pm.frac_outcome_test;
    
    %     BOXCAR RESCALE
    stim_dwell = stim_dwell * timescale; 
    dec_dwell = dec_dwell * timescale;
    outcome_dwell = outcome_dwell * timescale;
    extra_steps = extra_steps * timescale;
    
test_sequence = genTest(pm.pre_train_length, frac_outcome_test,outcome_nrns,n_nrns,len_trial,stim_nrns,stim_dwell,dec_nrns,stimulus);

function prototype = genTest(pre_train_length, frac_outcome_test, outcome_nrns,n_nrns,len_trial,stim_nrns,stim_dwell,dec_nrns,stimulus)
    % Initialize an empty matrix with dimensions number of neurons by
    % length of trial
    prototype = zeros(n_nrns, len_trial);
    disp("length of trial")
    disp(len_trial)
    r_picks = randi([1,2000],2*dec_nrns,pre_train_length);
    for pre_train_timestep = 1: pre_train_length
        for i = 1:2*dec_nrns
            neuron = r_picks(i, pre_train_timestep);
            prototype(i, neuron) = 1;
        end
    end
    
    % Generating the stimulus matrices for each stimulus pair
    %(AB)
    if(stimulus == 0)
        stim_block = [ones(2*stim_nrns, stim_dwell);zeros(stim_nrns,stim_dwell)];
    %(BC)
    elseif(stimulus == 1)
        stim_block = [zeros(stim_nrns,stim_dwell);ones(2*stim_nrns, stim_dwell)];
    % (CA)
    elseif(stimulus == 2)
        stim_block = [ones(stim_nrns, stim_dwell);zeros(stim_nrns, stim_dwell);ones(stim_nrns,stim_dwell)];
    end  
    
    % these variables act as moving pointers to indicate
    % where the next pattern of firing neurons is located
    initial_stim_nrn = 1;
    final_stim_nrn = 3*stim_nrns;
    start = pre_train_length + 1;
    stop = pre_train_length + stim_dwell;
    
    % Fills in stimulus neurons
    % x-axis (start:stop) is the duration of stimulus neurons firing
    % y-axis (initial_stim_nrn:final_stim_nrn) are the neurons firing
    prototype(initial_stim_nrn:final_stim_nrn, start:stop) = stim_block;
    
    %Sets up outcome neuron firing position (y-axis)
    initial_outcome_nrn = (final_stim_nrn) + (3*dec_nrns) + 1;
    %disp(initial_outcome_nrn)
    final_outcome_nrn = initial_outcome_nrn + outcome_nrns-1;
    %disp(final_outcome_nrn)
    %Creates a block of ones for the portion of the positive outcome neurons that
    %get turned on
    induced_att_block = ones(outcome_nrns,len_trial-pre_train_length);
    off_prob = induced_att_block .* rand(size(induced_att_block)); 
   % disp(off_prob)
    noise_indices =  off_prob > frac_outcome_test & off_prob ~= 0;
   % disp(noise_indices)
    induced_att_block(noise_indices) = 0;
    
    % Fills in portion of positive outcome 
    prototype(initial_outcome_nrn:final_outcome_nrn, start:len_trial) = induced_att_block;
end 
end