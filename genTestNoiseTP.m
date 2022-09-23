%Colton Bogucki
%Similar to genNoiseTrace, this function applies the on/off noises during
%testing for TP to specific time durations to ensure stable activity. 
function sequence = genTestNoiseTP(prototypes,on_noise,off_noise)
    
    global stimulus_duration
    
    %---------------------
    % Preconditions
    %---------------------
    if prod(prod((prototypes==0)+(prototypes==1)))==0
        error('error in myfunc genNoise: all elements of prototype must be either 0 or 1')
    elseif on_noise > 1 ||on_noise < 0 || off_noise > 1 || off_noise < 0
        error('error in myfunc genNoise: on_noise and off_noise must be numbers between 0 and 1')
    end

    %-------------------------------
    % Implementation
    %-------------------------------
    
    len_train=size(prototypes,2); %length of training
    sequence=prototypes; 
    
    for t=1:len_train        
        current_prototype=prototypes(:,t); 
        
        %ON NOISE%
        if t<stimulus_duration
            
            %if proportion of 0's is greater than the desired on_noise, use
            %probability distribution to turn on neurons at desired noise level
            off_indices = current_prototype == 0; %complements (every 1 turns to 0, all 0s turn to 1s)
            on_prob = off_indices .* rand(size(current_prototype));
            noise_indices =  on_prob < on_noise & on_prob ~= 0;
            sequence(noise_indices,t) = 1;


end