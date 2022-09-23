%Dana Schultz, Jul 21 2020 old initialise

%IMPORTANT
%before running, check noise, activity control
%make sure constant activity control with externals

function parameters = InitialiseParameters(learningRule, boxcar, length)%network_construction_seed,z_0_seed,varargin
%constant parameters
%for TP
if learningRule == 1
    post_pre_rule_status = 0;
else
    post_pre_rule_status = 1;
end 

if boxcar == 1
    boxcar_scales = [1]; % should sum to 1
    boxcar_scales_inh = [1];
elseif boxcar == 2
    boxcar_scales = [0.5 0.5];%[0.5 0.5],... % should sum to 1
    boxcar_scales_inh = [0.5 0.5];
    
elseif boxcar == 4
    boxcar_scales = [0.25 0.25 0.25 0.25];%[0.5 0.5],... % should sum to 1
    boxcar_scales_inh = [0.25 0.25 0.25 0.25];
end    
pairs_to_learn = 1; % later incorporated into the params. struct, TP section
padding = 0; % # of trials for stabilization

testing_reps_per_stim = 5;

%for progressive, feed in the 
switch pairs_to_learn
    case 1
        training_block = [25]; %NEED TO CHANGE block length for progressive paradigm
    case 2
        training_block = [15 8 5 5 3 1];%[15 8 5 3 3 1]; %NEED TO CHANGE block length for progressive paradigm
    case 3
        training_block = [20 12 9 9 4 4 2 1 1];
end
total_trials = padding + pairs_to_learn*sum(training_block) + pairs_to_learn*testing_reps_per_stim;
%disp(total_trials)

pre_train_length = 15;

parameters = struct(...
    %{
    'cycle_limit',... [10; ... %length of each trial 100; ... %number of
    trials 1; ... 1],... 'lpmi',                         1,... loop
    position of main iteration; i.e. if external input changes to next time
    step when inc_vect is [0,1,0] then lpmi is 2.
    %}
    ...%************************
    ...% transverse patterning *
    ...%************************
    'pre_train_length',             pre_train_length,...
    'n_training_trials',            total_trials-pairs_to_learn*testing_reps_per_stim,...
    'toggle_tp',                    true,...%if this is true, TP will be implemented instead of the simple sequence/trace conditioning
    'staged',                       false,...%this will add a field to the data structure which specifies the training/testing sequences based on the staged paradigm approach
    'toggle_progressive',           true,...%set training paradigm to Progressive
    'testing_reps_per_stim',        testing_reps_per_stim,...
    'pairs_to_learn',               pairs_to_learn,...%if 1, only AB; if 2, AB & BC; if 3, AB, BC & AC
    'stim_nrns',                    length,...
    'dec_nrns',                     2*length,...
    'outcome_nrns',                 2*length,...
    'stim_dwell',                   5,...
    'dec_dwell',                    5,...
    'outcome_dwell',                5,...
    'frac_outcome_test',            0.2,... %need to be implemented, possibly apply noise
    'extra_steps',                  0,...
    'length_of_each_trial',         15 + pre_train_length,...
    ...
    ... % activity stablizer toggle
    'toggle_stablize_activity',     0,...   %pre run network without synaptic modification
    'stabilization_period',         30,... %number of trials for network activity stabilization with no learning
    ... % staged paradigm block specification
    'padding',                      padding,... %padding trials for activity to stablize
    'block_reps',                   training_block,...
    ... % visualization toggles
    'toggle_figures',               true,...
    'toggle_display_failure_mode',  false,...
    'toggle_disp_param',            false,...
    ...
    ... % data recording settings
    'toggle_plot_activity',         true,...          
    'toggle_save_all_data',         true,...
    'toggle_record_success',        true,...
    'toggle_record_training',       true,...
    'toggle_record_y',              true,... % y and z dimensions: # neurons x timesteps x trials
    'toggle_record_z_train',        true,... 
    'toggle_record_z_test',         true,...
    'toggle_record_weights_exc',    true,...
    'toggle_record_weights_inh',    true,...
    'toggle_record_k0',             true,...
    'toggle_record_kff',            true,...
    ...
    ... % stopwatch toggle
    'toggle_fxn_stopwatch',         true,... % toggle true to record elapsed time
    ...
    ... % general settings
    'toggle_training',              true,...
    'toggle_testing' ,              true,...
    'number_of_trials',             total_trials,... 
    'consecutive_successes_to_halt',5,... % after this many successes-in-a-row, halt the program (this saves time); will only run if the previous toggle is on.
    ...
    ... % boxcar settings
    'toggle_divisive_refractory',   true,...
    'boxcar_refractory_period',     1,... %when set to 1 there is no refractory period
    'boxcar_window',                boxcar,... %excitation boxcar
    'boxcar_scales',                boxcar_scales,...%[0.5 0.5],... % should sum to 1
    'boxcar_window_inh',            boxcar,...
    'boxcar_scales_inh',            boxcar_scales_inh,...
    ...
    ... % random number seeds
    'network_construction_seed',    randi(100,1),... %network_construction_seed
    'z_0_seed',                     randi(100,1),... %z_0_seed won't be used
    ...
    ... % network topology settings
    'number_of_neurons',            uint16(2000),...
    'connectivity',                 0.1,... %fan in is exact i.e. each postsyn. nrn is enervated by same number of nrns
    'desired_mean_z',               0.1,...
    ...
    ... % straight sequence input settings
    'toggle_external_input',        true,... 
    'ext_activation',               uint16(30),... %default 30
    'stutter',                      uint16(1),...
    'shift',                        uint16(15),... %default 15
    'n_patterns',                   uint16(4),... %!!!!!!!
    ...
    ... % externals and noise
    'percent_ext_of_act',           0.1, ... %indicates the percentages of external activation among the percentage of total activations
    'extra_timesteps_train',        0,... % add this many timesteps with no external activation to the end of each trial; use test_length_of_each_trial to add timesteps during testing
    'on_noise',                     0.0,... %probability between 0-1
    'off_noise',                    0.22,... %0.15 %probability between 0-1
    'test_off_noise',               0.1,... %added by phyl, just for the attractor,"fractional_outcome_on"
    'test_on_noise',                0.0,... %new on noise has not been implemented, 
    ...
    ... % testing parameters
    'test_length_of_each_trial',    15 + pre_train_length,... %set this equal to stutter * number of patterns (for simple sequence)
    'first_stimulus_length',        uint16(10),... %in general, set this equal to stutter * boxcar window
    ...
    ... % trace conditioning
    'trace_interval',               uint16(0),... % (this can be deleted because user sets trace interval length in genTrace.m) 
    'toggle_trace',                 false,... %true trace, false normal; 
    ...                                         %If trace is toggled, there are special options, including noise, success 
    ...                                         %See genTrace.m, genTraceNoise.m, DetermineTrialSuccessTrace.m 
    ... % modification rates
    'toggle_k0_preliminary_mod',    false,...%adjust k0 before training begins to attune to desired activity 
    'toggle_k0_training_mod',       true,...%modify k0 at the end of each trial to maintain desired activity
    'epsilon_pre_then_post',        0.04,....005,...0.0015
    'epsilon_post_then_pre',        0.04,...
    'epsilon_feedback',             0.15,...%changes per timestep %11-11-21, set to 0.04, nov-22, 0.048 for TF
    'epsilon_k_0',                  0.2,... %modification rate used to adapt k0 before and/or during training; make sure appropriate toggles above are turned on
    'epsilon_k_ff',                 0.2,...%kff will not be updated unless pre then post is off
    ... 
    ... % synaptic modification settings
    'toggle_pre_then_post',         true,...
    'toggle_post_then_pre',         post_pre_rule_status,... %!!!!!!!
    'toggle_stutter_e_fold_decay',  false,... %set to false when changing decay rates
    'fractional_mem_pre',           0.78,... %1101 0.78 %nmda; used with pre_then_post, controls Zbar pre decay, -1/log(fractional) = exponential time constant
    'fractional_mem_post',          0.78,... %maybe: make third quad less effective %.78 1101 0.78 spiketiming; used with post_then_pre, controls Zbar post decay
    'offset_pre_then_post',         0,...
    'offset_post_then_pre',         0,...
    ...
    ... % starting values for inhibitory variables
    'k_0_start',                    0.52,... %0.47,...  %future: might program k_0/k_ff = constant
    'k_fb_start',                   0.0475,...
    'k_ff_start',                   0.0009,...%0.0083,... 0.0066 ...
    ...
    ... % weight settings
    'toggle_rand_weights',          false,... %toggle between random and uniform excitatory weights
    'weight_start',                 0.4,... %starting value if using uniform weights
    'weight_high',                  0.2,... %upper limit for weight values if using random distribution
    'weight_low',                   0.8,... %lower limit for weight values if using random distribution
    'weight_inhib_start',           1 ... %choose starting weight value for inhibitory synapse; default 1
    ...
);
if nargin>0
    parameters = WriteToStruct(parameters);%parameters,varargin{:}
end
parameters = AdjustParameters(parameters);
end

function parameters = AdjustParameters(parameters)
    
    if parameters.toggle_stutter_e_fold_decay
        stutter = double(parameters.stutter);
        parameters.fractional_mem_pre = exp(-1/(stutter-2));
        parameters.fractional_mem_post = exp(-1/(stutter-2));
    end
    
    if parameters.toggle_trace == true
        % User must manually set trial length here (see genTrace.m)
        parameters.length_of_each_trial = 43; 
    elseif parameters.toggle_tp == true
        %User will manually specify trial length here
        parameters.length_of_each_trial = 15 + parameters.pre_train_length;
    else
        % Default (simple sequence): sets trial length = patterns * stutter
        parameters.length_of_each_trial = parameters.n_patterns * ...
            (parameters.stutter + parameters.trace_interval); %ignore trace_interval here (not used)
    end

    %network and topology settings
    n_nrns = double(parameters.number_of_neurons);
    parameters.n_fanin = round(n_nrns*parameters.connectivity);
    %set range of neurons that are to be viewed during the simulation
            
    if (parameters.toggle_tp == true)
        stim = parameters.stim_nrns;
        dec = parameters.dec_nrns;
        out = parameters.outcome_nrns;
        parameters.nrn_viewing_range = (3*stim+3*dec+2*out);
    else
        parameters.nrn_viewing_range = parameters.number_of_neurons;
    end
       
end






