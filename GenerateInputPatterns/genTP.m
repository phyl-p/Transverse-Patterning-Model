% Colton Bogucki (2/14/21)
% This file produces a training input sequence to simulate transverse
% patterning.
% Arguments:
% pm = parameters from InitialiseParameters
% stimulus = either 0,1, or 2. 
%   0 means the stimulus pair is AB
%   1 means the stimulus pair is BC
%   2 means the stimulus pair is CA
% outcome = either 0 or 1
%   0 means the outcome is negative (-)
%   1 means the outcome is positive (+)

function input_sequence = genTP(pm,stimulus,outcome)

%Phyl (09-27-21): decalred these Global variables in InitialiseParameters
%Phyl (11-11-21): changed these Global variables to params in InitialiseParameters

%Global variables
%stim_nrns = number of neurons in stimuli (A,B,C)
%dec_nrns = number of neurons in decision (a,b,c)
%outcome_nrns = number of neurons in yes/no (+,-)

% stim_nrns = pm.stim_nrns;
% dec_nrns = pm.dec_nrns
% outcome_nrns = pm.outcome_nrns;

%stim_duration = dwell time per stimulus
%dec_duration = dwell time per decision
%outcome_duration = dwell time per outcome

% stim_dwell = pm.stim_dwell;
% dec_dwell = pm.dec_dwell;
% outcome_dwell = pm.outcome_dwell;

%Logicals for whichever subtask is being generated
%AB = (AB)a or (AB)b
%BC = (BC)b or (BC)c
%CA = (CA)c or (CA)a

%global AB
%global BC
%global CA

%Logicals for whichever decision is made
%global plus
%global minus

%The combination of the two sets of logicals above create a matrix of
%possibilities for the externally activated sequence:
%     |--------|--------|
%     |   +    |   -    | 
%|----|--------|--------|
%| AB | (AB)a+ | (AB)b- |
%|----|--------|--------|
%| BC | (BC)b+ | (BC)c- |
%|----|--------|--------|
%| CA | (CA)c+ | (CA)a- |
%-----------------------|


%-------------------
    % PARAMETERS 
    %-------------------
    % Parameters from InitialiseParameters.m:  
    n_nrns          = pm.number_of_neurons;
    timescale       = pm.boxcar_window; %for boxcar adjustments
    len_trial       = pm.length_of_each_trial;
    
    % Parameters to be set by user:
    %IF THESE ARE CHANGED, SIMILAR CHANGES MUST BE MADE IN:
    %DetermineTrialSuccessTP - beginning
    %RunModel - RecordTrial and TestNetwork functions - lines on spy diagram
    %ScrollThroughData - StandardDisplay function - lines on spy diagram
    %MainSingleTask - nrn_viewing_rng = [1 stim_nrns+dec_nrns+outcome_nrns]
    stim_nrns           = pm.stim_nrns; %number of stimulus neurons
    dec_nrns            = pm.dec_nrns; %number of decision neurons
    outcome_nrns        = pm.outcome_nrns; %number of correctness neurons
    stim_dwell          = pm.stim_dwell; %number of timesteps stimulus neurons are on
    dec_dwell           = pm.dec_dwell; %number of timesteps decision neurons are on
    outcome_dwell       = pm.outcome_dwell; %number of timesteps correctness neurons are on
    extra_steps         = pm.extra_steps; %padding zeros at the end of training sequence
    pre_train_length    = pm.pre_train_length; 
%disp("pm.outcom_nrns in genTP")
%disp(outcome_nrns)
    % Noise options 
    % ***Right now, noise is the same as the simple sequence*** 
       
%     % training noises
%     off_noise_during_stim       = false; 
%     off_noise_during_dec        = false; 
%     off_noise_during_correct    = false; 
%     on_noise_during_stim        = false; 
%     on_noise_during_dec         = false; 
%     on_noise_during_correct     = false; 
%     
%     % testing noises
%     off_noise_during_stim_test      = false; 
%     off_noise_during_dec_test       = false; 
%     off_noise_during_correct_test   = false; 
%     on_noise_during_stim_test       = false;
%     on_noise_during_dec_test        = false; 
%     on_noise_during_correct_test    = false; 
    
    toggle_produce_visualisation_of_input = false;
    toggle_save_input_sequence_to_file    = false;
    
%     BOXCAR RESCALE
    stim_dwell = stim_dwell * timescale; 
    dec_dwell = dec_dwell * timescale;
    outcome_dwell = outcome_dwell * timescale;
    extra_steps = extra_steps * timescale;
%     len_trial = len_trial * timescale;
    

    %-------------------------------
    % Function Preconditions
    %-------------------------------
    if n_nrns < (3*stim_nrns + 3*dec_nrns + 2*outcome_nrns)
        disp('n_nrns must be greater than or equal to 3*stim_nrns + 3*dec_nrns + 2*outcome_nrns. Try increasing n_nrns, decreasing stim_nrns, decreasing dec_nrns, or decreasing correct_nrns')
    end 
%     if len_trial ~= pm.length_of_each_trial
%         error(['Check InitialiseParameters.m to ensure parameter.length_of_each_trial is correct. Trial length should equal stim_dwell + dec_dwell + correct_dwell. The correct length is ', num2str(len_trial), '. Find or Ctrl+F: parameter.length_of_each_trial']);
%     end
%     if len_trial ~= pm.test_length_of_each_trial
%         error(['Check InitialiseParameters.m to ensure test_length_of_each_trial is correct. The test trial length should equal stim_dwell + dec_dwell + correct_dwell. The correct length is ', num2str(len_trial), '. Find or Ctrl+F: test_length_of_each_trial']);
%     end
%     if stim_dwell ~= pm.first_stimulus_length
%         error('Check InitialiseParameters.m under testing to ensure the first stimulus length is correct. first_stimulus_length should equal stim_dwell. Find or Ctrl+F: first_stimulus_length');
%     end
    
    %-----------------------------------------------
    % Implementation (calls function genPrototypes) 
    %-----------------------------------------------
    input_sequence = genPrototypes(pre_train_length, n_nrns,len_trial,stim_nrns,stim_dwell,dec_nrns,dec_dwell,outcome_nrns,outcome_dwell,stimulus,outcome);
    if extra_steps > 0
        % padding of zeros at the end
        input_sequence = [input_sequence, zeros(n_nrns,extra_steps)];
    end
    
    %-------------------------------
    % Visualize Sequence
    %-------------------------------
     if toggle_produce_visualisation_of_input == true
        figure
        spy(input_sequence)
     end
     
    %-------------------------------
    % Save Sequence
    %-------------------------------
    if toggle_save_input_sequence_to_file == true
        save('Input.mat','input_sequence')
    end

end
function prototypes = genPrototypes(pre_train_length, n_nrns,len_trial,stim_nrns,stim_dwell,dec_nrns,dec_dwell,outcome_nrns,outcome_dwell,stimulus,outcome)
%GENPROTOTYPES: This function does the main implementation of the
%transverse patterning sequence. It is called by genTP.  
    
    % Initialize an empty matrix with dimensions number of neurons by
    % length of trial
    prototypes = zeros(n_nrns, len_trial);
    
    %Here are the concatenations of matrices which produce the correct
    %firing blocks of neurons:
    %AB stimulus = [ones(2*stim_nrns, stim_dwell);zeros(stim_nrns,stim_dwell)];
    %BC stimulus = [zeros(stim_nrns,stim_dwell);ones(2*stim_nrns, stim_dwell)];
    %CA stimulus = [ones(stim_nrns, stim_dwell);zeros(stim_nrns, stim_dwell);ones(stim_nrns,stim_dwell)];
    %Decision a = [ones(dec_nrns, dec_dwell);zeros(2*dec_nrns,dec_dwell)];
    %Decision b = [zeros(dec_nrns, dec_dwell);ones(dec_nrns,dec_dwell);zeros(dec_nrns, dec_dwell)];
    %Decision c = [zeros(2*dec_nrns, dec_dwell);ones(dec_nrns,dec_dwell)];
    %+ outcome = [ones(outcome_nrns, outcome_dwell);zeros(outcome_nrns,outcome_dwell)];
    %- outcome = [zeros(outcome_nrns, outcome_dwell);ones(outcome_nrns, outcome_dwell)];
    
    % Generating the matrices for each input sequence
    %(AB) a +
    r_picks = randi([1,2000],2*dec_nrns,pre_train_length);
    for pre_train_timestep = 1: pre_train_length
        for i = 1:2*dec_nrns
            neuron = r_picks(i, pre_train_timestep);
            prototypes(i, neuron) = 1;
        end
    end
    
    if(stimulus == 0 && outcome == 1)
        stim_block = [ones(2*stim_nrns, stim_dwell);zeros(stim_nrns,stim_dwell)];
        dec_block = [ones(dec_nrns, dec_dwell);zeros(2*dec_nrns,dec_dwell)];
        outcome_block = [ones(outcome_nrns, outcome_dwell);zeros(outcome_nrns,outcome_dwell)];
        %disp("AB")
        %disp(outcome_block)
    %(AB) b -
    elseif(stimulus == 0 && outcome == 0)
        stim_block = [ones(2*stim_nrns, stim_dwell);zeros(stim_nrns,stim_dwell)];
        dec_block = [zeros(dec_nrns, dec_dwell);ones(dec_nrns,dec_dwell);zeros(dec_nrns, dec_dwell)];
        outcome_block = [zeros(outcome_nrns, outcome_dwell);ones(outcome_nrns, outcome_dwell)];
    %(BC) b +
    elseif(stimulus == 1 && outcome == 1)
        stim_block = [zeros(stim_nrns,stim_dwell);ones(2*stim_nrns, stim_dwell)];
        dec_block = [zeros(dec_nrns, dec_dwell);ones(dec_nrns,dec_dwell);zeros(dec_nrns, dec_dwell)];
        outcome_block = [ones(outcome_nrns, outcome_dwell);zeros(outcome_nrns,outcome_dwell)];
    % (BC) c -
    elseif(stimulus == 1 && outcome == 0)
        stim_block = [zeros(stim_nrns,stim_dwell);ones(2*stim_nrns, stim_dwell)];
        dec_block = [zeros(2*dec_nrns, dec_dwell);ones(dec_nrns,dec_dwell)];
        outcome_block = [zeros(outcome_nrns, outcome_dwell);ones(outcome_nrns, outcome_dwell)];
    % (CA) c +
    elseif(stimulus == 2 && outcome == 1)
        stim_block = [ones(stim_nrns, stim_dwell);zeros(stim_nrns, stim_dwell);ones(stim_nrns,stim_dwell)];
        dec_block = [zeros(2*dec_nrns, dec_dwell);ones(dec_nrns,dec_dwell)];
        outcome_block = [ones(outcome_nrns, outcome_dwell);zeros(outcome_nrns,outcome_dwell)];
    % (CA) a -
    elseif(stimulus == 2 && outcome == 0)
        stim_block = [ones(stim_nrns, stim_dwell);zeros(stim_nrns, stim_dwell);ones(stim_nrns,stim_dwell)];
        dec_block = [ones(dec_nrns, dec_dwell);zeros(2*dec_nrns,dec_dwell)];
        outcome_block = [zeros(outcome_nrns, outcome_dwell);ones(outcome_nrns, outcome_dwell)];
    end
    
    % these variables act as moving pointers to indicate
    % where the next pattern of firing neurons is located
    initial_stim_nrn = 1;
    final_stim_nrn = 3*stim_nrns;
    start = pre_train_length + 1;
    stop = pre_train_length + stim_dwell;
    
    
    % --- Main Training Sequence Generator ---

    % Fills in stimulus neurons
    % x-axis (start:stop) is the duration of stimulus neurons firing
    % y-axis (initial_cs_nrn:final_cs_nrn) are the neurons firing
    prototypes(initial_stim_nrn:final_stim_nrn, start:stop) = stim_block;

    % Stimulus -> Decision neurons
    % Sets up start and stop for decision neurons (along x-axis)
    start = start + stim_dwell;
    stop = stop + dec_dwell;
    
    % Sets up decision neurons firing location (along y-axis)
    initial_dec_nrn = final_stim_nrn + 1;
    final_dec_nrn = initial_dec_nrn + 3*dec_nrns - 1;
    
    % Fills in decision neurons
    % x-axis (start:stop) is the duration of decision neurons firing
    % y-axis (initial_dec_nrn:final_dec_nrn) are the neurons firing
    prototypes(initial_dec_nrn:final_dec_nrn, start:stop) = dec_block;
    
    % Decision -> outcome neurons
    % Sets up start and stop for outcome neurons (along x-axis)
    start = start + dec_dwell;
    stop = stop + outcome_dwell;

    % Sets up outcome neurons firing location (along y-axis)
    initial_outcome_nrn = final_dec_nrn + 1;
    final_outcome_nrn = initial_outcome_nrn + 2*outcome_nrns - 1;

    % Fill in outcome neurons  
    % x-axis (start:stop) is the duration of outcome neurons firing
    % y-axis (initial_outcome_nrn:final_outcome_nrn) are the neurons firing
    prototypes(initial_outcome_nrn:final_outcome_nrn, start:stop) = outcome_block;
%     figure
%     spy(prototypes)
    end
    

