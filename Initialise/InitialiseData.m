function [data,pm] = InitialiseData(pm)
    if pm.toggle_tp == 1 && pm.toggle_progressive == 1
        block = pm.block_reps;
        phases = length(block);
        switch pm.pairs_to_learn                           
                case 1     % Up
                    training = sum(block);
                    testing = 1*pm.testing_reps_per_stim;
                case 2     % Down
                    training = sum(block) + sum(block(0.3*phases+1:end));
                    testing = 2*pm.testing_reps_per_stim;
                case 3     % Left
                    training = sum(block) + sum(block(0.6*phases+1:end)) + ...
                            sum(block(0.6*phases+1:end));
                    testing = 3*pm.testing_reps_per_stim;           
        end
        n_trials = training + testing;
        disp(training)
        disp(testing)
        disp("n_trials calculated")
        disp(n_trials)
    end
    
    n_nrns = pm.number_of_neurons;
    len_trial = pm.length_of_each_trial;
    n_trials = pm.number_of_trials;
    test_len_trial = pm.test_length_of_each_trial;

    data = struct(...
        'z_train_last_trial',           zeros(n_nrns,len_trial,'logical')   , ...
        'z_test_last_trial',            zeros(n_nrns,len_trial,'logical')   , ...
        'mean_z_train_last_trial' ,     zeros(len_trial,1               )   , ...
        'mean_z_train_all_trials' ,     zeros(len_trial,n_trials        )   , ...
        'mean_z_test_all_trials' ,      zeros(len_trial,n_trials         )   , ...
        'k0_history'     ,              zeros(n_trials,1                )   , ...
        'kff_history'    ,              zeros(n_trials,1                )   , ...
        'success_vect'   ,              zeros(n_trials,1                )   , ...
        'max_consecutive_successes' ,   0 ...
    );

    if pm.toggle_record_z_train == true
        data.z_train    =               zeros(n_nrns,len_trial,n_trials,'logical');
    end
    if pm.toggle_record_z_test == true
        data.z_test     =               zeros(n_nrns,test_len_trial,n_trials,'logical');
    end
    if pm.toggle_record_y == true
        data.y_train    =               zeros(n_nrns,len_trial,n_trials);
        data.y_train_last_trial =       zeros(n_nrns,len_trial,1);
    end
    if pm.toggle_record_weights_exc == true
        data.w_excite_trial     =       zeros(n_nrns,pm.n_fanin,n_trials);
    end
    if pm.toggle_record_weights_inh == true
        data.w_inhib_trial      =       zeros(n_nrns,n_trials);
    end
    %adds a matrix to the data structure which controls the staged paradigm
    %implementation
    if pm.toggle_stablize_activity == true && pm.toggle_tp == true
        %disp("gen stabilize")
        data.stabilize       = activityParadigm(data, pm.stabilization_period);
    end

     %initialize stimuli setup for the training and testing sequences
    if pm.toggle_tp == true && pm.staged == true
        data.staged             =       stagedParadigm3(pm, n_trials); 
        %disp('training with Staged Paradigm')
    elseif pm.toggle_tp == true && pm.toggle_progressive == true && pm.staged == false
        data            =       progressiveParadigm(pm, data); 
        %disp('training with Progressive Paradigm')
    end
    %disp(data)
end