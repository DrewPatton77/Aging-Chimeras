
% ======================================================================
% 
% Running this evaluates the classification type 
%  (synchronous, chimera, asynchronous) and the synchronous pattern by
%  first binarizing the pair-wise cognitive system synchrony matrices 
%  by applying a threshold synchrony value (S >= threshold), then running
%  the Louvain algorithm on these binarized matrices.
%
% ======================================================================
codeNums = [355581,487062,602004,634143,750313,807122];
main_path = cd;
main_path = strsplit(main_path,'\\Code');
main_path = string(main_path{1});

for k = 1:6
    setCode = codeNums(k)
    path = main_path + "/Data/"+setCodes(k)+"_Structural_Connectome.mat";
    Data = importdata(path);

    %% Gather all Synchronisation Matrices
    num_subjects = size(Data.S,4);
    num_regions = size(Data.S,3);
    num_Parcellations = size(Data.S,1);
    thresh = 0.65;
    A = {}; % Empty Cell Array
    count = 1;
    for subject = 1:num_subjects
        for k = 1:num_regions
            C = Data1.S(:,:,k,subject);
            A{count}=(C >= thresh); % Binarize the synchronisation matrices for some threshold, \rho_thres = 0.8
            count = count + 1;
        end
    end
    %%
    gamma = 0.7;
    N=length(A{1});
    T=length(A);
    B=spalloc(N*T,N*T,N*N*T+2*N*T);
    twomu=0;
    for s = 1:T
        k = sum(sum(A{s})) / size(A{s},1);
        twom = sum(sum(A{s}));
        indx=[1:N]+(s-1)*N;
        if twom == 0
            B(indx,indx) = 0;
        else
            B(indx,indx) = A{s}-gamma*(k^2/twom)*ones(size(A{s}));
        end
    end
    %cd(main_path+"/Code/Louvain-Algorithm");
    [S,~] = genlouvain(B);
    S = reshape(S,N,T);
    
    %%
    Patterns = zeros(num_Parcellations,num_regions,num_subjects);
    count = 1;
    for subject = 1:num_subjects
        for k = 1:num_regions
            Patterns(:,k,subject) = S(:,count);
            count = count + 1;
        end
    end
    
    %% Classify the results and save the classification
    Classification = zeros(num_regions,num_subjects);
    count = 1;
    for subject = 1:num_subjects
        for k = 1:num_regions
            if length(unique(S(:,count))) == 1 % If only one unique value, then set the classification as synchronous.
                Classification(k,subject) = 3;
            elseif length(unique(S(:,count))) == num_Parcellations % If all are unique values, then set the classification as asynchronous.
                Classification(k,subject) = 1;
            else
                Classification(k,subject) = 2; % Otherwise, set it as a chimera state.
            end
            count = count + 1;
        end
    end
    
    %% Save the data.
    Data.Patterns = Patterns;
    Data.Classification = Classification;
    
    save(path,'Data')
    
    %% Get the synchronous patterns in terms of 1's (belongs as a synchronous cognitive system) and 0's (belongs as an asynchronous cognitive system)
    Numbered_Patterns = Patterns;
    Patterns = zeros(size(Numbered_Patterns));
    num_Parcellations = size(Data1.S,1);
    for i = 1:size(Numbered_Patterns,2)
        for j = 1:size(Numbered_Patterns,3)
            P = Numbered_Patterns(:,i,j);
            [P_unique,~,cols] = unique(P);
            
            P_count = accumarray(cols,1); % Count the number occurrences.
            Pattern = zeros(num_Parcellations,1);
            for l = 1:length(P_unique)
                if P_count(l) > 1 % If occurrence is greater than 1, there's a synchronized grouping.
                    Pattern = Pattern + ismember(P,P_unique(l));
                end
            end
            Patterns(:,i,j) = Pattern;
        end
    
    end
    
    %% Save the data.
    Data.Binarized_Patterns = Patterns;
    save(path,'Data')
end
