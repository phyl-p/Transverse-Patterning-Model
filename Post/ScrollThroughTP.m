%ScrollThroughTP - Colton Bogucki 4/3/21
%This file allows for the success detection process for transverse
%patterning to be better visualized by displaying the decision blocks for
%each trial.
f = figure;
trial = 1;
while true %ie. loop forever, but can terminate by closing the figure window and inducing an error
    clf
    
    StandardDisplay(parameters,data,trial)
    
    waitforbuttonpress
    if f.CurrentCharacter == '' % is the symbol for left arrow
        trial = trial-1;
    else
        trial = trial + 1;
    end
    trial = mod(trial-1,size(data.mean_z_train_all_trials,2))+1;
end

%Creates a 4x1 subplot 
%|----------------|---------------|------------------|----------------|
%| training trial | testing trial | correct decision | wrong decision |
%|----------------|---------------|------------------|----------------|
function StandardDisplay(parameters,data,trial)
    %constants for drawing lines
    stim = parameters.stim_nrns;
    dec = parameters.dec_nrns;
    out = parameters.outcome_nrns;
    len = parameters.length_of_each_trial +1;
    total_num_nrns = 3*stim+3*dec+2*out+100;
    
    sub_frame_num = 3;
    sub_frame = 0;
    %training trial
    title(trial)
    if trial <= sum(parameters.block_reps)
        sub_frame = sub_frame + 1;
        sub_frame_num = sub_frame_num + 1;  
        subplot(1,sub_frame_num, sub_frame)
        spy(data.z_train(1:total_num_nrns,:,trial)) %spy(data.z_train(:,:,trial))
        %ticks at the separations between stimuli, decisions, and outcome nrns
        yticks([1 stim 2*stim 3*stim 3*stim+dec 3*stim+2*dec 3*stim+3*dec 3*stim+3*dec+out 3*stim+3*dec+2*out])
        title('training')
        line([0 len],[stim+0.5 stim+0.5],'LineStyle','-','Color','k');
        line([0 len],[2*stim+0.5 2*stim+0.5],'LineStyle','-','Color','k');
        line([0 len],[3*stim+0.5 3*stim+0.5],'LineStyle','-','Color','k');
        line([0 len],[3*stim+dec+0.5 3*stim+dec+0.5],'LineStyle','-','Color','k');
        line([0 len],[3*stim+2*dec+0.5 3*stim+2*dec+0.5],'LineStyle','-','Color','k');
        line([0 len],[3*stim+3*dec+0.5 3*stim+3*dec+0.5],'LineStyle','-','Color','k');
        line([0 len],[3*stim+3*dec+out+0.5 3*stim+3*dec+out+0.5],'LineStyle','-','Color','k');
        line([0 len],[3*stim+3*dec+2*out+0.5 3*stim+3*dec+2*out+0.5],'LineStyle','-','Color','k');
    end 
    %testing trial
    sub_frame = sub_frame + 1;
    subplot(1,sub_frame_num,sub_frame)
    spy(data.z_test(1:total_num_nrns,:,trial)) %spy(data.z_test(:,:,trial))
    yticks([1 stim 2*stim 3*stim 3*stim+dec 3*stim+2*dec 3*stim+3*dec 3*stim+3*dec+out 3*stim+3*dec+2*out])
    title('testing')
    line([0 len],[stim+0.5 stim+0.5],'LineStyle','-','Color','k');
    line([0 len],[2*stim+0.5 2*stim+0.5],'LineStyle','-','Color','k');
    line([0 len],[3*stim+0.5 3*stim+0.5],'LineStyle','-','Color','k');
    line([0 len],[3*stim+dec+0.5 3*stim+dec+0.5],'LineStyle','-','Color','k');
    line([0 len],[3*stim+2*dec+0.5 3*stim+2*dec+0.5],'LineStyle','-','Color','k');
    line([0 len],[3*stim+3*dec+0.5 3*stim+3*dec+0.5],'LineStyle','-','Color','k');
    line([0 len],[3*stim+3*dec+out+0.5 3*stim+3*dec+out+0.5],'LineStyle','-','Color','k');
    line([0 len],[3*stim+3*dec+2*out+0.5 3*stim+3*dec+2*out+0.5],'LineStyle','-','Color','k');
    
    stimulus = data.staged(1,trial);
    [b1,b2,r,w] = DetermineTrialSuccessTP(parameters,data,trial,stimulus);
    r = num2str(r);
    w = num2str(w);
    
    %title for graphs
    switch stimulus
        %AB
        case 0
            t1 = {'a';['Neuron Firings: ',r]};
            t2 = {'b';['Neuron Firings: ',w]};
            n1 = ['Neuron ',num2str(3*stim+1)];
            n2 = ['Neuron ',num2str(3*stim+dec)];
            n3 = ['Neuron ',num2str(3*stim+dec+1)];
            n4 = ['Neuron ',num2str(3*stim+2*dec)];
        %BC
        case 1
            t1 = {'b';['Neuron Firings: ',r]};
            t2 = {'c';['Neuron Firings: ',w]};
            n1 = ['Neuron ',num2str(3*stim+dec+1)];
            n2 = ['Neuron ',num2str(3*stim+2*dec)];
            n3 = ['Neuron ',num2str(3*stim+2*dec+1)];
            n4 = ['Neuron ',num2str(3*stim+3*dec)];
        %CA
        case 2
            t1 = {'c';['Neuron Firings: ',r]};
            t2 = {'a';['Neuron Firings: ',w]};
            n1 = ['Neuron ',num2str(3*stim+2*dec+1)];
            n2 = ['Neuron ',num2str(3*stim+3*dec)];
            n3 = ['Neuron ',num2str(3*stim+1)];
            n4 = ['Neuron ',num2str(3*stim+dec)];
    end
    
    total_length = parameters.length_of_each_trial*parameters.boxcar_window; 
    sub_frame = sub_frame +1;
    %correct decsision block
    subplot(1,sub_frame_num,sub_frame)
    %size(b1)
    spy(b1)
    yticks([1 dec+1])
    yticklabels({n1,n2})
    xticks([1 total_length])
    xticklabels({'1',num2str(total_length)})
    title(t1)
    
    %incorrect decision block
    sub_frame = sub_frame + 1;
    subplot(1,sub_frame_num,sub_frame)
    %size(b2)
    spy(b2)
    yticks([1 dec+1])
    yticklabels({n3,n4})
    xticks([1 total_length])
    xticklabels({'1',num2str(total_length)})
    title(t2)
    
    success = data.success_vect(trial);
    if trial>0.84*parameters.number_of_trials
        if success == 0
            str = ['Trial ',num2str(trial),': Failure'];
        else
            str = ['Trial ',num2str(trial),': Success'];
        end
    else
        str = ['Trial ',num2str(trial)];
    end
    sgtitle(str);
    
    
end