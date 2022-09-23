f = figure;

trial = 1;


format=3;

%format options:
%1 --> StandardDisplay()
    %displays training and testing trials as is, along with activity levels
%2 --> ReorderDisplay()
    %displays training and testing trials with rows sorted by order of
    %first timestep of firing
%3 --> show both standard and reordered displays
if format == 3
    set(gcf,'position',[250,250,900,500])
end

while true %ie. loop forever, but can terminate by closing the figure window and inducing an error
    clf
    
    GeneratePlots(parameters,data,format,trial)
    
    waitforbuttonpress
    if f.CurrentCharacter == '' % is the symbol for left arrow
        trial = trial-1;
    else
        trial = trial + 1;
    end
    trial = mod(trial-1,size(data.mean_z_train_all_trials,2))+1;
end



function GeneratePlots(parameters,data,format,trial)
%{
    %{
    n_trials = parameters.number_of_trials;
    len_trial = parameters.length_of_each_trial;
    mean_z_test = zeros(len_trial,n_trials);
    mean_z_test(:,:) = mean(data.z_test,1);
    %}
    mean_z_train = data.mean_z_train;
    
    
    
    %subplot(1,2,1)
    %spy(data.z_test(:,:,trial))
    
    %subplot(1,2,2)
    plot(mean_z_train(:,trial),'b')
    %{
    hold on
    plot(mean_z_test (:,trial),'r')
    %}
%    ylim([0,ceil(max(max([mean_z_test,mean_z_train]))*100)/100])
    ylim([0,ceil(max(max(mean_z_train))*100)/100])
    title(['trial ',num2str(trial),' (test in red) success = ',num2str(data.success_vect(trial))])
%}
    %spy(data.z_test(:,:,trial))
    %subplot(1,2,1)
    
    
    switch format
        case 1
            StandardDisplay(parameters,data,trial)
        case 2
            ReorderDisplay(data,trial)
        case 3
            CombinedDisplay(parameters,data,trial)
            
    end
    %{
    
    %}
    
    
    %{
    [reordered_all,i_reorder] = reorder(data.z_train(:,:,trial));
    reordered_ext = zeros(size(reordered_all));
    input=variables.input_prototypes(1:size(reordered_all,1),:);
    reordered_ext(i_reorder,:) = input;
    
    display = zeros([size(reordered_all),3]);
    display(:,:,1) = reordered_ext;
    display(:,:,2) = reordered_ext+reordered_all;
    display(:,:,3) = reordered_all;
    imshow(1-display)
    %}
    
    if parameters.toggle_post_then_pre * parameters.epsilon_post_then_pre > 0
        rulename = 'Bid';
    else
        rulename = 'Uni';
    end
    
    result = data.success_vect(trial);
    
    n_pat = parameters.n_patterns;
    %len_trial = parameters.length_of_each_trial;
    
    result_name = char(FailureModeDictionary(result));
    %result_id = FailureModeDictionary(result_name);
    switch result_name
        case {'success','unknown failure mode', 'noise', 'overcompression'}
            message = result_name;
        case {'activity died out','attractor','context attractor'}
            last_pat = ResultCodeToPatternsLearned(result,n_pat);
            %attractor_pat = dual_decimal_code - last_pat;
            message = [result_name, ', learned to pattern',...
                                num2str(last_pat)];
        otherwise
            error(['exception occured, error name unfamiliar: ', result_name])
    end
    
    title({[rulename,'trial ',num2str(trial)], message})
       
end

function StandardDisplay(parameters,data,trial)
    subplot(2,3,[1,4])
    spy(data.z_train(:,:,trial))
    if(parameters.toggle_tp == true)
        hold on
        stim = 25;
        dec = 50;
        out = 50;
        yticks([1 stim 2*stim 3*stim 3*stim+dec 3*stim+2*dec 3*stim+3*dec 3*stim+3*dec+out 3*stim+3*dec+2*out])
        len = params.length_of_each_trial + 1;
        line([0 len],[stim+0.5 stim+0.5],'LineStyle','-','Color','k');
        line([0 len],[2*stim+0.5 2*stim+0.5],'LineStyle','-','Color','k');
        line([0 len],[3*stim+0.5 3*stim+0.5],'LineStyle','-','Color','k');
        line([0 len],[3*stim+dec+0.5 3*stim+dec+0.5],'LineStyle','-','Color','k');
        line([0 len],[3*stim+2*dec+0.5 3*stim+2*dec+0.5],'LineStyle','-','Color','k');
        line([0 len],[3*stim+3*dec+0.5 3*stim+3*dec+0.5],'LineStyle','-','Color','k');
        line([0 len],[3*stim+3*dec+out+0.5 3*stim+3*dec+out+0.5],'LineStyle','-','Color','k');
        line([0 len],[3*stim+3*dec+2*out+0.5 3*stim+3*dec+2*out+0.5],'LineStyle','-','Color','k');
    end
    title('training')
    hold off
    
    
    subplot(2,3,[2,5])
    spy(data.z_test(:,:,trial))
    if(parameters.toggle_tp == true)
        hold on
        stim = 32;
        dec = 32;
        out = 32;
        len = params.length_of_each_trial + 1;
        yticks([1 stim 2*stim 3*stim 3*stim+dec 3*stim+2*dec 3*stim+3*dec 3*stim+3*dec+out 3*stim+3*dec+2*out])
        line([0 len],[stim+0.5 stim+0.5],'LineStyle','-','Color','k');
        line([0 len],[2*stim+0.5 2*stim+0.5],'LineStyle','-','Color','k');
        line([0 len],[3*stim+0.5 3*stim+0.5],'LineStyle','-','Color','k');
        line([0 len],[3*stim+dec+0.5 3*stim+dec+0.5],'LineStyle','-','Color','k');
        line([0 len],[3*stim+2*dec+0.5 3*stim+2*dec+0.5],'LineStyle','-','Color','k');
        line([0 len],[3*stim+3*dec+0.5 3*stim+3*dec+0.5],'LineStyle','-','Color','k');
        line([0 len],[3*stim+3*dec+out+0.5 3*stim+3*dec+out+0.5],'LineStyle','-','Color','k');
        line([0 len],[3*stim+3*dec+2*out+0.5 3*stim+3*dec+2*out+0.5],'LineStyle','-','Color','k');
    end
    title('testing')
    
    subplot(2,3,3)
    plot(mean(data.z_train(:,:,trial),1))
    ylim([0,parameters.desired_mean_z*2])
    line([0,parameters.length_of_each_trial],...
        repmat(parameters.desired_mean_z,[1,2]),'Color','red')
    xlim([0,parameters.length_of_each_trial])
    title({'training trial mean activity,', 'with desired activity in red'})
    
    subplot(2,3,6)
    plot(mean(data.z_test(:,:,trial),1))
    ylim([0,parameters.desired_mean_z*2])
    line([0,parameters.length_of_each_trial],...
        repmat(parameters.desired_mean_z,[1,2]),'Color','red')
    xlim([0,parameters.length_of_each_trial])
    title({'testing trial mean activity,', 'with desired activity in red'})
end

function ReorderDisplay(data,trial)
    z_train = data.z_train(:,:,trial);
    z_test  = data.z_test(:,:,trial);
    [z_train_reordered,i_reorder] = reorder(z_train);
    subplot(1,2,1)
    spy(z_train_reordered)
    subplot(1,2,2)
    spy(z_test(i_reorder,:))
end

function CombinedDisplay(parameters,data,trial)
    subplot(2,5,[1,6])
    spy(data.z_train(1:256,:,trial)) %spy(data.z_train(:,:,trial))
    title('training')
    if(parameters.toggle_tp == true)
        stim = 32;
        dec = 32;
        out = 32;
        len = parameters.length_of_each_trial +1;
        yticks([1 stim 2*stim 3*stim 3*stim+dec 3*stim+2*dec 3*stim+3*dec 3*stim+3*dec+out 3*stim+3*dec+2*out])
        line([0 len],[stim stim],'LineStyle','-','Color','k');
        line([0 len],[2*stim 2*stim],'LineStyle','-','Color','k');
        line([0 len],[3*stim 3*stim],'LineStyle','-','Color','k');
        line([0 len],[3*stim+dec 3*stim+dec],'LineStyle','-','Color','k');
        line([0 len],[3*stim+2*dec 3*stim+2*dec],'LineStyle','-','Color','k');
        line([0 len],[3*stim+3*dec 3*stim+3*dec],'LineStyle','-','Color','k');
        line([0 len],[3*stim+3*dec+out 3*stim+3*dec+out],'LineStyle','-','Color','k');
        line([0 len],[3*stim+3*dec+2*out 3*stim+3*dec+2*out],'LineStyle','-','Color','k');
    end
    
    subplot(2,5,[2,7])
    spy(data.z_test(1:256,:,trial)) %spy(data.z_test(:,:,trial))
    title('testing')
    if(parameters.toggle_tp == true)
        stim = 32;
        dec = 32;
        out = 32;
        yticks([1 stim 2*stim 3*stim 3*stim+dec 3*stim+2*dec 3*stim+3*dec 3*stim+3*dec+out 3*stim+3*dec+2*out])
        len = parameters.length_of_each_trial +1;
        line([0 len],[stim stim],'LineStyle','-','Color','k');
        line([0 len],[2*stim 2*stim],'LineStyle','-','Color','k');
        line([0 len],[3*stim 3*stim],'LineStyle','-','Color','k');
        line([0 len],[3*stim+dec 3*stim+dec],'LineStyle','-','Color','k');
        line([0 len],[3*stim+2*dec 3*stim+2*dec],'LineStyle','-','Color','k');
        line([0 len],[3*stim+3*dec 3*stim+3*dec],'LineStyle','-','Color','k');
        line([0 len],[3*stim+3*dec+out 3*stim+3*dec+out],'LineStyle','-','Color','k');
        line([0 len],[3*stim+3*dec+2*out 3*stim+3*dec+2*out],'LineStyle','-','Color','k');
    end
    
    %%% trace viewing helpers for testing plot %%%
    global success_timestep_1;
    global success_timestep_2;
    global cs_duration;
    global trace_duration;
    if parameters.toggle_trace    %helpers for viewing trace conditioning    
        start_of_us_nrns = cs_duration + trace_duration + 1; %timestep of US onset
        xline(success_timestep_1);
        xline(success_timestep_2);
        xline(start_of_us_nrns, '-','US onset');
        trace_result = DetermineTrialSuccessTrace(parameters,data,trial);
        xlabel(trace_result);
    end   
    
    subplot(2,5,3)
    plot(mean(data.z_train(:,:,trial),1))
    ylim([0,parameters.desired_mean_z*2])
    line([0,parameters.length_of_each_trial],...
        repmat(parameters.desired_mean_z,[1,2]),'Color','red')
    xlim([0,parameters.length_of_each_trial])
    title({'training trial mean activity,', 'with desired activity in red'})
    
    subplot(2,5,8)
    plot(mean(data.z_test(:,:,trial),1))
    ylim([0,parameters.desired_mean_z*2])
    line([0,parameters.length_of_each_trial],...
        repmat(parameters.desired_mean_z,[1,2]),'Color','red')
    xlim([0,parameters.length_of_each_trial])
    title({'testing trial mean activity,', 'with desired activity in red'})
    
    %reordered display
    z_train = data.z_train(:,:,trial);
    z_test  = data.z_test(:,:,trial);
    [z_train_reordered,i_reorder] = reorder(z_train);
    subplot(2,5,[4,9])
    spy(z_train_reordered)
    subplot(2,5,[5,10])
    spy(z_test(i_reorder,:))
    
%     set(gcf,'position',[237,216,897,501])
end

