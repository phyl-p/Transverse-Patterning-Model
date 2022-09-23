%Dana Schultz, 13 Jul 2020 
%Jin Lee, 15 Jan 2019, "Main.m"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [resultAnalysis, data,variables,parameters] = MainSingleTask(learningRule, boxcar, External_input_length, computer) %seed_pair for sweeps
% Initialise and organise variables that will be used by 'RunModel.m'

%important functions are stored in subfolders 
%for the sake of organisation. addpath() adds these folders to the list of
%directories accessed by Matlab.
addpath('./Initialise')
addpath('./GenerateInputPatterns')
addpath('./Helpers')
addpath('./Boxcar')
addpath('./PresetKs')
addpath('./DataAnalysis')


%{
'parameters' is a struct containing values that do not change;
'variables' is a struct containing values that might change.
  this distinction is arbitrary, and some variables such as theta were
placed in the 'variables' struct just in case they might need to be 
manipulated in the future.

'data' is a struct used to store the history of changes in values
  e.g.
data.z_history is a 3 dimensional array that contains the ith neuron's 
output in the jth timestep during the kth trial 
in the ith row, jth column, kth page.
%}
%{
A word on structs:
As far as the code in this demo is concerned, structs are data structures
that allow variables to be stored hierarchically. To illustrate, they can
be thought of as containers; all values that are contained in these
containers are referred to as fields.
one can store arrays and scalar values in
these structs, but just as a container can contain other subcontainers, it
is possible to use a struct to store another struct. The following code
illustrates.

>> a_struct = struct;
>> a_struct.a_scalar = 3.141582;
>> a_struct.a_matrix = [1,0;0,1];
>> a_struct
a_struct = 

  struct with fields:

    a_scalar: 3.1416
    a_matrix: [2�2 double]

>> b_struct = struct;
>> b_struct.b_subfield1 = 2.7183;
>> b_struct.b_subfield2 = [0,0];
>> a_struct.a_field = b_struct;
>> a_struct
a_struct = 

  struct with fields:

    a_scalar: 3.1416
    a_matrix: [2�2 double]
     a_field: [1�1 struct]
>> a_struct.a_field
ans = 

  struct with fields:

    b_subfield1: 2.7183
    b_subfield2: [0 0]


%}
%{
all of the values used in the simulation are stored in these structs to
make the code cleaner, and are acceessed as fields. For example, if I
wanted to access 'weights_excite' which is a matrix containing the
excitatory weights (for a network using a fanin connectivity scheme in the
case of this demo), I would first note that the excitatory weights are
changed by the network after each timestep and must therefore be in the
'variables' struct. The matrix can be accessed as 
'variables.weights_excite'
%}
% tic
disp(datestr(now, 'dd/mm/yy-HH:MM'))
% changes to boxcar_window will result in rescaling of the following: 
% desired_mean_z, stutter, epsilon_weight_nmda/spiketiming/feedback, 
% decay_nmda/spiketiming, boxcar_refractory_period

%Sets initial values for all parameters. Edit the file for this function to
%make desired changes.
parameters = InitialiseParameters(learningRule, boxcar, External_input_length); %seed_pair(1),seed_pair(2) for sweeps

%**************************
%parameters.off_noise = noise;
%parameters.test_off_noise = noise;
%**************************

%{
%     'boxcar_window'                 , 4,... %default 1
%     'number_of_neurons'             , 401,... %default 1000
%     'consecutive_successes_to_halt' , 0,...
%     'toggle_figures'                , true,...
%     'toggle_spiketiming'            , true,... %!!!!!!!
%     'k_ff_start'                    , 0.0066,... %for bx4, use 0.04
%     'desired_mean_z'                , 0.1,...
%     'epsilon_k_0'                   , 0,... %default 0.5
%     'n_patterns'                    , 4,...
%     'toggle_record_z_train'         , true,...
%     'toggle_record_z_test'          , true,...
%     'toggle_record_k0'              , true,...
%     'toggle_record_weights'         , false... %record excitatory weights after each update
%     );
%Initalise Core Vars and Data are run in Initialise_PresetBoxcar
%variables = InitialiseCoreVars(parasumeters);
%data = InitialiseData(parameters);

% npat = 28; k0 = 0.5; kff = 0.0085; ktable = table(npat,k0,kff); save PresetKs/ktable ktable timescale
%}
if parameters.toggle_fxn_stopwatch
    tic
end
seed = parameters.network_construction_seed;
disp(['seed 1 = ',num2str(seed)])
rng(seed)

retrieve_k = false;
%viewing range set to fit TP model.

%%
[parameters,variables,data] = Initialise_PresetBoxcar(parameters,retrieve_k);
%%
    if parameters.toggle_progressive
        paradigm = "progressive";
    else
        paradigm = "staged";
    end

    if parameters.toggle_post_then_pre
        rules = "TT";
    else rules = "TF";
    end
    
    if computer == "PC"
        path = "E:";
    elseif computer == "mac"
        path = "/Volumes/Phyl";
    end

    folder = strcat(path,"/MATLAB/Levy Lab/v9boxcar_trace_tp (2)/v9boxcar_trace_tp/Data/","boxcar", num2str(boxcar), ...
        "/externalInputLength",num2str(External_input_length),"/",paradigm,'/',rules,'/',num2str(sum(parameters.n_training_trials)), ...
        "/desiredActivity", num2str(parameters.desired_mean_z));
   disp(folder)
    if ~exist(folder, 'dir')
       mkdir(folder)
    end
    stab = "";
if parameters.toggle_stablize_activity == true
    [variables,data] = ActivityStabilizer(parameters,variables,data);
 
    if parameters.toggle_plot_activity
        figure;
        time = 1:parameters.stabilization_period;
        xlabel('Trial Number');
        ylabel('Activity Level (%)');
        plot(time, mean(data.mean_z_train_all_trials(:,1:parameters.stabilization_period),1));
        yline(parameters.desired_mean_z); %plot desired mean z as a horizontal line
        title(strcat("external input per timestep: ", num2str(2*External_input_length), " ext. fraction: ", sprintf('%.4f',2*External_input_length/parameters.number_of_neurons),'/',num2str(sum(parameters.n_training_trials))));
        filename = strcat(folder,"/",datestr(now, 'mm-dd-yy_HH:MM'),"_Boxcar",num2str(parameters.boxcar_window),"stabilzer",num2str(parameters.stabilization_period),".fig");
        
        dim = [0.5 0.5 0.3 0.3];
        rule = {rules,'desired activity level:', sprintf('%.4f',parameters.desired_mean_z) };
        annotation('textbox',dim,'String',rule,'FitBoxToText','on');

        saveas(gcf, filename)
    end
else 
    stab = "ON";
end
%%

%parameters.epsilon_k_0 = 0;
%variables.k_ff = 0.02;
%having prepared all the variables in the 3 structs, the information needed
%for running the simulation is now passed onto the function 'RunModel'
%[variables,data] = RunModel(parameters,variables,data);
if parameters.toggle_fxn_stopwatch
    toc
end

[variables,data] = RunModel(parameters,variables,data);

if parameters.toggle_plot_activity
    figure;
    time = 1:data.training_length;
    plot(time, mean(data.mean_z_train_all_trials(:,1:data.training_length),1));
    xlabel('Trial Number');
    ylabel('Activity Level (%)');
    yline(parameters.desired_mean_z); %plot desired mean z as a horizontal line
    title(strcat("external input per timestep: ", num2str(2*External_input_length), " ext. fraction: ", sprintf('%.4f',2*External_input_length/parameters.number_of_neurons)));
    filename = strcat(folder,"/",datestr(now, 'mm-dd-yy_HH:MM'),"_Boxcar",num2str(parameters.boxcar_window),".fig");
    saveas(gcf, filename)
end

resultAnalysis = countDecision(parameters,data,variables);

if parameters.toggle_save_all_data
   
    pm_filename = strcat(folder, "/",datestr(now, 'mm-dd-yy_HH:MM'),"_parameters_Boxcar",num2str(parameters.boxcar_window),".mat");
    save(pm_filename,'parameters');

    data_filename = strcat(folder,"/",datestr(now, 'mm-dd-yy_HH:MM'),"_data_Boxcar",num2str(parameters.boxcar_window),".mat");
    save(data_filename,'data');
    
    resultAnalysis_filename = strcat(folder,"/",datestr(now, 'mm-dd-yy_HH:MM'),"_results_Boxcar",num2str(parameters.boxcar_window),".mat");
    save(resultAnalysis_filename,'resultAnalysis')

    %save recent run to workspace
    save('recentrun')
try
    evalin('base','load(''recentrun.mat'')')
catch
    warning('could not load variables to workspace')
end
end


