function getWC()
% =====================================================================
%
% Running this takes the structural connectome (streamline density) data
%  and simulates brain activity via the Wilson-Cowan model and returns
%  the pair-wise cognitive system synchrony matrix (S).
%
% ======================================================================

% Cognitive Systems Assignment.
CS = containers.Map; % Matlab equivalent of Python's dictionary.
CS('Att') = [29,30,60,61,62,88,89,117,118,119]; % Attention Cognitive System
CS('Aud') = [63,64,69,120,121,126]; % Auditory
CS('CO')  = [16,38,39,52,53,54,55,65,66,73,98,99,110,111,112,122,123]; % Cingulo-Opercular
CS('FP')  = [17,22,23,33,40,67,70,71,74,79,80,81,92,93,100,124,127,128]; % FrontoParietal
CS('mDM') = [26,45,50,51,56,57,58,59,84,104,108,109,113,114,115,116]; % Medial Default Mode
CS('MS')  = [37,42,43,44,46,47,48,49,97,102,103,105,106,107]; % Motor and Somatosensory
CS('SC')  = [1,2,3,4,5,6,7,8,9,10,11,12,13,14]; % Subcortical
CS('VT')  = [15,19,20,21,24,25,34,35,36,68,72,76,77,78,82,83,94,95,96,125,126]; % Ventral Temporal Association
CS('V')   = [18,27,28,31,32,41,75,85,86,87,90,91,101]; % Visual

% Get research group codes.
setCodes = ["355581","487062","602004","634143","750313","807122"];

main_path = cd;
main_path = strsplit(main_path,'\\Code');
main_path = string(main_path{1});
for k = 1:6
    
    % Import data
    path = main_path+"/Data/"+setCodes(k)+"_Structural_Connectome.mat";
    Data = importdata(path); % All the data for the respective dataset.
    
    % Sizes
    N = size(Data.Streamline_Density,3); % Number of subjects for that research group.
    M = size(Data.Streamline_Density,2); % Number of brain regions (n=128).
            
    % Scaling
    scale = 110; % Fine-tune such that we "just before" the active state.
            
    % Store
    S = zeros(9,9,M,N);
            
    % Run loop over all subjects in the research group.
    for subject = 1:N % Subject Number
        
        % Retrieve structural connectome data.
        G = Data.Streamline_Density(:,:,subject); % Structural connectome.
        D = Data.Time_Delays(:,:,subject); % Time-delay matrix.
        c5 = scale*Data.c5_crit(subject) - 100; % coupling strength parameter (c_E) --> tuned to "just before" transition point.
        
        % Run loop over each stimulated brain region.
        for stimulated_oscillator = 1:M % Perturb each brain region
            % Run the Wilson-Cowan dynamics and store.
            Dynamics = get_WC_Dynamics(G/scale, D, stimulated_oscillator, c5, CS); % Function to initialize Wilson-Cowan Dynamics & get synchony matrix (S).
            S(:,:,stimulated_oscillator,subject) = Dynamics.S; % Pair-wise cognitive system synchrony matrix.
        end
        % Save
        Data.S = S; % Store to Data.
        save(path,"Data"); % Save newly updated Data.
    end
end

quit()
end
