function sequence = genNoise(pm, prototypes,on_noise,off_noise)
    
    %first_stimulus is passed in 
    %Dana Schultz, 26 May 2020
    %takes a noiseless prototype sequence and introduces noise
    %-------------------------------
    % Preconditions
    %-------------------------------
    if prod(prod((prototypes==0)+(prototypes==1)))==0
        error('error in myfunc genNoise: all elements of prototype must be either 0 or 1')
    elseif on_noise > 1 ||on_noise < 0 || off_noise > 1 || off_noise < 0
        error('error in myfunc genNoise: on_noise and off_noise must be numbers between 0 and 1')
%     elseif floor(on_noise)~=on_noise || floor(off_noise)~=off_noise
%         error('error in myfunc genNoise: on_noise or off_noise must be integers')
    end
    %-------------------------------
    % Implementation
    %-------------------------------

    len_train=size(prototypes,2);   %length of training
    sequence=prototypes;
    
  
    for t = 1:len_train
            %select firing pattern at time t
            current_prototype=prototypes(:,t);

            %% new implementation (probability values instead of ints)
            %if proportion of 0's is greater than the desired on_noise, use
            %probability distribution to turn on neurons at desired noise level
            
            %on noise for off indices (?)
            off_indices = current_prototype == 0; %complements (turns every 1 in current prototype  to 0, all 0s turn to 1s) What's the point?
            on_prob = off_indices .* rand(size(current_prototype));
            noise_indices =  on_prob < on_noise & on_prob ~= 0;
            sequence(noise_indices,t) = 1;


            %off noise. repeat same process as above, turning off neurons at desired level
            off_prob = current_prototype .* rand(size(current_prototype)); 
            noise_indices =  off_prob < off_noise & off_prob ~= 0;
            sequence(noise_indices,t) = 0;
            
            
%             %% Oct. 28 21 version of manual inplementation
%             num_on_neurons = pm.number_of_neurons*pm.desired_mean_z*pm.percent_ext_of_act*(1-on_noise);
%             num_off_neurons = 2*pm.stim_nrns - num_on_neurons;
%             %implement noise
%             vec = [ones(1,num_on_neurons), zeros(1,pm.stim_nrns-num_on_neurons)];
%             vec = vec(randperm(25));
%             
%             off_indices = current_prototype == 0; %complements (turns every 1 in current prototype  to 0, all 0s turn to 1s) What's the point?
%             on_prob = off_indices .* rand(size(current_prototype));
%             noise_indices =  on_prob < on_noise & on_prob ~= 0;
%             sequence(noise_indices,t) = 1;
% 
% 
%             %repeat same process as above, turning off neurons at desired level
%             off_prob = current_prototype .* rand(size(current_prototype)); 
%             noise_indices =  off_prob < off_noise & off_prob ~= 0;
%             sequence(noise_indices,t) = 0;
    end

    
end