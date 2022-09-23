
function vars = InitialiseCoreVars(params)
    
    %generate input prototypes (trace sequence, tp sequence, or simple sequence);
    %all inputs in the simulation will be noisy versions of the following prototype:
    if params.toggle_trace == true 
        vars.input_prototypes = genTrace(params); 
    elseif params.toggle_tp == true
        vars.input_prototypes = 0;
        %genTP(params);
    else
        vars.input_prototypes = genSequence(params); %default is simple sequence
    end
        
    %pre check if modifications to k_0 will be well behaved, ie no negative
    %values or divergence to k_0
    n_active_nrns = double(params.ext_activation); 
    n_nrns=params.number_of_neurons; 
    timescale = numel(params.boxcar);
%     if params.desired_mean_z<(double(n_active_nrns)/(timescale*double(n_nrns)))
%         error('desired_mean_z is too low, will result in infinite expansion of k_0')
%     end
    
    %mutable numeric values, such as inhibition and excitation    
    vars.k_0 = params.k_0_start;
    vars.k_fb = params.k_fb_start;
    vars.k_ff = params.k_ff_start;
    vars=initVarArrays(vars,params);
    vars.success = false;
        
    %weights for excitation and feedback inhibition
    vars.weights_feedback_inhib = ones(params.number_of_neurons,1) * params.weight_inhib_start;
    vars = getWeights(vars,params);
    
end    
    
    
function vars = initVarArrays(vars,params)
    timescale = numel(params.boxcar);
    number_of_nrns = params.number_of_neurons;
    len_trial = params.length_of_each_trial;
    vars.y=zeros(number_of_nrns,1);
    vars.z=zeros(number_of_nrns,1,'logical');
    if params.offset_pre_then_post < 0 || params.offset_post_then_pre < 0
        error('offset values must be positive integers')
    else
        z_memory = max(params.offset_pre_then_post+1,params.offset_post_then_pre+1); %changed the second argument in max() to params.offset_post_then_pre+1 (4/11/21)
        vars.z_prev = zeros(number_of_nrns,z_memory); %params.offset_pre_then_post);
    end
    vars.excite=zeros(number_of_nrns,1);
    vars.inhib=zeros(number_of_nrns,1);
    %BOXCAR ADJUSTMENTS
    bcvect = params.boxcar;
    vars.boxcar_exc = zeros(number_of_nrns,length(bcvect));
    vars.boxcar_inh = zeros(1,params.boxcar_window_inh);
    if (numel(params.boxcar_scales) == params.boxcar_window) && (sum(params.boxcar_scales) == 1)...
            && (numel(params.boxcar_scales_inh) == params.boxcar_window_inh) ...
            && (sum(params.boxcar_scales) == 1)
        vars.boxcar_scales = params.boxcar_scales;
        vars.boxcar_scales_inh = params.boxcar_scales_inh;
    else
        error('boxcar_scales must sum to 1 and must contain the same number of elements as boxcar_window')
    end
    vars.boxcar_refractory_period = uint16(floor(timescale * params.boxcar_refractory_period));
    %activity
    vars.mean_z=0;
    vars.mean_z_history_during_trial=zeros(1,len_trial);
    %time scale constants
    vars.z_pre_then_post_decay=zeros(number_of_nrns,1);
    vars.z_post_then_pre_decay=zeros(number_of_nrns,1);
    vars.spike_threshold = 0.5; %params.spike_threshold;
    vars.dw_fbexcite=zeros(number_of_nrns,params.n_fanin);
    vars.mean_dw=0;
    if params.toggle_testing == 1
        vars.z_test = zeros(number_of_nrns,1);
    end
end

function vars = getWeights(vars,params)
    number_of_neurons = params.number_of_neurons;
    if params.toggle_rand_weights == 0
        vars.weights_excite=ones(number_of_neurons,params.n_fanin)*params.weight_start;
    else
        %a matrix with n_nrn rows, n_fanin columns with entries that
        %are drawn from a uniform random distribution on the interval
        %between params.weight_low and params.weight_high
        vars.weights_excite=rand(number_of_neurons,params.n_fanin)*(params.weight_high-params.weight_low)+params.weight_low;
    end
    vars.connections_fanin=zeros(number_of_neurons,params.n_fanin);
    for row = 1:number_of_neurons
        c_row=randperm(number_of_neurons,params.n_fanin);
        vars.connections_fanin(row,:) = c_row;
    end
end

