function [parameters,variables,data] = Initialise_PresetBoxcar(parameters,retrieve_k)
    timescale=parameters.boxcar_window;
    parameters.boxcar = ones(1,timescale);
%     parameters.pre_then_post_offset = timescale;  %only affects the number of columns of z_prev which in turn changes the z_nmda term to z_nmda(t-nmdatimeshift)
    %%%%%%%%%%%%%%%%%%%%%%

    variables = InitialiseCoreVars(parameters);
    %NOTE: kff,k0 and kfb are not rescaled;
    %{
    InitialiseCoreVars already transfers kff,k0 and kfb from params to vars
    therefore rescaling has no effect if done on params.k_*_start
    %}
    [parameters,variables] = RescaleParameters(parameters,variables,timescale);

    %at each timestep, we want n_active_nrns/timescale inputs to be active,
    %therefore we must turn off n_active_nrns-n_active_nrns/timescale inputs
    
%     parameters.off_noise = floor(parameters.n_active_nrns*(1-1/timescale));
%     if parameters.off_noise >= parameters.n_active_nrns
%         error('offnoise too high relative to number of active neurons; reduce boxcar window')
%     end


    %initialise data struct used for data collection
    [data,parameters] = InitialiseData(parameters);

    %generate new K constants or retrieve old ones from file
    param_backup = parameters;
    
    if ~exist('retrieve_k','var')
        retrieve_k = false;
    end
    if class(retrieve_k) == "logical" && retrieve_k == false
        if parameters.epsilon_k_0 ~= 0 && parameters.toggle_k0_preliminary_mod
            [variables.k_0,variables.k_ff,~] = PresetK(parameters,variables,data);
            %[variables.k_0,variables.k_ff,variables.weights_feedback_inhib] = PresetK(parameters,variables,data);
        end
    else
        try
            if class(retrieve_k) == "logical" && retrieve_k == true
                [variables.k_0,variables.k_ff] = RetrievePresetK(parameters);
            else %otherwise assume retrieve_k is some string that specifies a filename
                [variables.k_0,variables.k_ff] = RetrievePresetK(parameters,retrieve_k);
            end
        catch ME
            rethrow(ME)
            %error('need to create appropriate preset K table (that matches both timescale and number of patterns) and save to PresetKs/ktable.mat')
        end
    end
    
    %restore original parameter settings
    parameters = param_backup;
    
    if ~parameters.toggle_k0_training_mod
        parameters.epsilon_k_0 = 0; %turn off k0 modification after preliminary adaptation is complete
    end

    %disp(variables.k_0)
    %disp(variables.k_ff)
    %turn off modification for k_0 and k_ff


    %variables = variables_backup

    %parameters.toggle_spiketiming = false;
end