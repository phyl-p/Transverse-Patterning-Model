clear
% Parameters to be set by user:
    stim_nrns           = 16; %number of stiumulus neurons
    dec_nrns            = 32; %number of decision neurons
    outcome_nrns        = 32; %number of correctness neurons
    stim_dwell          = 5; %number of timesteps stimulus neurons are on
    dec_dwell           = 3; %number of timesteps decision neurons are on
    outcome_dwell       = 10; %number of timesteps correctness neurons are on
    extra_steps         = 0; %padding zeros at the end of training sequence
    len_trial           = stim_dwell + dec_dwell + outcome_dwell + extra_steps;
    n_nrns              = 3*stim_nrns + 3*dec_nrns + 2*outcome_nrns;
    AB                  = false; %subtask AB
    BC                  = false; %subtask BC
    CA                  = true; %subtask CA
    a                   = false; %decision a
    b                   = false; %decision b
    c                   = true; %decision c
    

    

input_sequence = genPrototypes(n_nrns,len_trial,stim_nrns,stim_dwell,dec_nrns,dec_dwell,outcome_nrns,outcome_dwell,AB,BC,CA,a,b,c);
figure
subplot(1,3,1)
spy(input_sequence)

AB=false;
BC=true;
input_sequence = genPrototypes(n_nrns,len_trial,stim_nrns,stim_dwell,dec_nrns,dec_dwell,outcome_nrns,outcome_dwell,AB,BC,CA,a,b,c);
subplot(1,3,2)
spy(input_sequence)

BC=false;
CA=true;
input_sequence = genPrototypes(n_nrns,len_trial,stim_nrns,stim_dwell,dec_nrns,dec_dwell,outcome_nrns,outcome_dwell,AB,BC,CA,a,b,c);
subplot(1,3,3)
spy(input_sequence)

     
function prototypes = genPrototypes(n_nrns,len_trial,stim_nrns,stim_dwell,dec_nrns,dec_dwell,outcome_nrns,outcome_dwell,AB,BC,CA,a,b,c)
%GENPROTOTYPES: This function does the main implementation of the trace
%conditioning sequence. It is called by genTrace.  
    
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
    if(AB*a)
        stim_block = [ones(2*stim_nrns, stim_dwell);zeros(stim_nrns,stim_dwell)];
        dec_block = [ones(dec_nrns, dec_dwell);zeros(2*dec_nrns,dec_dwell)];
        outcome_block = [ones(outcome_nrns, outcome_dwell);zeros(outcome_nrns,outcome_dwell)];
    %(AB) b -
    elseif(AB*b)
        stim_block = [ones(2*stim_nrns, stim_dwell);zeros(stim_nrns,stim_dwell)];
        dec_block = [zeros(dec_nrns, dec_dwell);ones(dec_nrns,dec_dwell);zeros(dec_nrns, dec_dwell)];
        outcome_block = [zeros(outcome_nrns, outcome_dwell);ones(outcome_nrns, outcome_dwell)];
    %(AB) c -
    elseif(AB*b)
        stim_block = [ones(2*stim_nrns, stim_dwell);zeros(stim_nrns,stim_dwell)];
        dec_block = [zeros(2*dec_nrns, dec_dwell);ones(dec_nrns,dec_dwell)];
        outcome_block = [zeros(outcome_nrns, outcome_dwell);ones(outcome_nrns, outcome_dwell)];
    %(BC) b +
    elseif(BC*b)
        stim_block = [zeros(stim_nrns,stim_dwell);ones(2*stim_nrns, stim_dwell)];
        dec_block = [zeros(dec_nrns, dec_dwell);ones(dec_nrns,dec_dwell);zeros(dec_nrns, dec_dwell)];
        outcome_block = [ones(outcome_nrns, outcome_dwell);zeros(outcome_nrns,outcome_dwell)];
    % (BC) a -
    elseif(BC*a)
        stim_block = [zeros(stim_nrns,stim_dwell);ones(2*stim_nrns, stim_dwell)];
        dec_block = [ones(dec_nrns, dec_dwell);zeros(2*dec_nrns,dec_dwell)];
        outcome_block = [zeros(outcome_nrns, outcome_dwell);ones(outcome_nrns, outcome_dwell)];
    % (BC) c -
    elseif(BC*c)
        stim_block = [zeros(stim_nrns,stim_dwell);ones(2*stim_nrns, stim_dwell)];
        dec_block = [zeros(2*dec_nrns, dec_dwell);ones(dec_nrns,dec_dwell)];
        outcome_block = [zeros(outcome_nrns, outcome_dwell);ones(outcome_nrns, outcome_dwell)];
    % (CA) c +
    elseif(CA*c)
        stim_block = [ones(stim_nrns, stim_dwell);zeros(stim_nrns, stim_dwell);ones(stim_nrns,stim_dwell)];
        dec_block = [zeros(2*dec_nrns, dec_dwell);ones(dec_nrns,dec_dwell)];
        outcome_block = [ones(outcome_nrns, outcome_dwell);zeros(outcome_nrns,outcome_dwell)];
    % (CA) a -
    elseif(CA*a)
        stim_block = [ones(stim_nrns, stim_dwell);zeros(stim_nrns, stim_dwell);ones(stim_nrns,stim_dwell)];
        dec_block = [ones(dec_nrns, dec_dwell);zeros(2*dec_nrns,dec_dwell)];
        outcome_block = [zeros(outcome_nrns, outcome_dwell);ones(outcome_nrns, outcome_dwell)];
    % (CA) b -
    elseif(CA*b)
        stim_block = [ones(stim_nrns, stim_dwell);zeros(stim_nrns, stim_dwell);ones(stim_nrns,stim_dwell)];
        dec_block = [zeros(dec_nrns, dec_dwell);ones(dec_nrns,dec_dwell);zeros(dec_nrns, dec_dwell)];
        outcome_block = [zeros(outcome_nrns, outcome_dwell);ones(outcome_nrns, outcome_dwell)];
    end

    % these variables act as moving pointers to indicate
    % where the next pattern of firing neurons is located
    initial_stim_nrn = 1;
    final_stim_nrn = 3*stim_nrns;
    start = 1;
    stop = stim_dwell;
    
    
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

    % Fill in correctnes neurons  
    % x-axis (start:stop) is the duration of outcome neurons firing
    % y-axis (initial_outcome_nrn:final_outcome_nrn) are the neurons firing
    prototypes(initial_outcome_nrn:final_outcome_nrn, start:stop) = outcome_block;
   
end
