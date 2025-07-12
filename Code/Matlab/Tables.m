
% =============================================================
%
% Takes each research group's data, evaluates the various network
%  measures, and organizes the results into a table for further
%  analysis.
%
% =============================================================


% Cognitive Systems Assignment
CS_names = ["Att","Aud","CO","FP","DM","MS","SC","V","VT"];
Att = [29,30,60,61,62,88,89,117,118,119];
Aud = [63,64,69,120,121,126];
CO = [16,38,39,52,53,54,55,65,66,73,98,99,110,111,112,122,123];
FP = [17,22,23,33,40,67,70,71,74,79,80,81,92,93,100,124,127,128];
mDM = [26,45,50,51,56,57,58,59,84,104,108,109,113,114,115,116];
MS = [37,42,43,44,46,47,48,49,97,102,103,105,106,107];
SC = [1,2,3,4,5,6,7,8,9,10,11,12,13,14];
V = [18,27,28,31,32,41,75,85,86,87,90,91,101];
VT = [15,19,20,21,24,25,34,35,36,68,72,76,77,78,82,83,94,95,96,125];

Cognitive_System = repmat("",128,1);
CS = [{Att} {Aud} {CO} {FP} {mDM} {MS} {SC} {V} {VT}];
for i = 1:length(CS)
    Cognitive_System(CS{i}) = CS_names(i);
end

% Lobal region assignment.
BL_names = ["Frontal","Parietal","Temporal","Occipital","Subcortical"];
Frontal = [37,97,110,52,16,73,56,57,58,59,113,114,115,116,111,112,53,54,55,74,17,38,98,100,40,99,39,88,89,29,30,33,92,93,105,106,107,46,47,48,49,67,124];
Parietal = [45,104,26,84,117,118,119,60,61,62,22,23,79,80,81,65,66,122,123,102,103,42,43,44,108,109,50,51];
Temporal = [128,127,70,71,120,121,63,64,34,35,94,95,82,83,24,25,15,72,20,21,77,78,69,126,76,19,125,68,36,96];
Occipital = [27,28,85,86,87,31,32,90,91,75,18,101,41];
%Cingulate = [110,52,16,73,45,104,26,84];
Subcortical = [1,2,3,4,5,6,7,8,9,10,11,12,13,14];

BL = [{Frontal},{Parietal},{Temporal},{Occipital},{Subcortical}];
Lobe = repmat("",128,1);
for i = 1:length(BL)
    Lobe(BL{i}) = BL_names(i);
end

% Load data.
main_path = cd;
main_path = strsplit(main_path,'\\Code');
main_path = string(main_path{1});

codeNums = [355581,487062,602004,634143,750313,807122]; % The research group number of each dataset

for k = 1:6
    path = main_path + "/Data/"+setCodes(k)+"_Structural_Connectome.mat";
    Data = importdata(data_path); % Load the dataset.
    W0 = Data.Streamline_Density; % The streamline density forming the structural connectome.
    W0(isnan(W0)) = 0; % For any NaN values -- there's a single NaN element at k = 5, ind = 355, i = j =76.
    c5 = Data.c5_crit; % Excitatory-to-excitatory coupling parameter value putting dynamics at critical regime.
    W0 = sum(W0.*reshape(c5,1,1,numel(c5)),4); % Weighted structural connectome by the critical c5.
    Age0 = Data.Age; % The mean age of the age-bracket each individual belongs to.
    Class0 = Data.Classification; % The synchrony class. 1>Async,2>Chimera,3>Sync.
    P = Data.Binarized_Patterns;
    for n = 1:size(P,2)
        for m = 1:size(P,3)
            if sum(P(:,n,m)) == 9
                if Class0(n,m) == 2
                    Class0(n,m) = 3;
                end
            end
        end
    end
    sub = size(W0,3); % The number of individuals.
    reg = size(W0,1); % The number of brain regions.
    rho_G0 = Data.rho_G;
    Chimera_Index0 = Data.CI;
    Metastability_Index0 = Data.MI;
    X0 = Data.X;
    Y0 = Data.Y;
    Z0 = Data.Z;

    % Find each dataset's age partition
    if k == 1
        Partition = [15,25,35,50,65];
    elseif k == 2
        Partition = [20,50,65,90];
    elseif k == 3
        Partition = [15,30,45,55,65,80];
    elseif k == 4
        Partition = [5,15,20,25,40,50,60,70,90];
    elseif k == 5
        Partition = [15,25,30,35,45,50,65];
    elseif k == 6
        Partition = [15,20,25,30,35,50,65];
    end
 
    Age0_Partitioned = zeros(length(Age0),1);
    for i = 1:length(Partition) 
        if i == 1
            continue
        end
        indx = (Age0 < Partition(i)) & (Age0 >= Partition(i-1));
        Age0_Partitioned(indx) = mean(Age0(indx));
    end
    
    Between = zeros(reg,sub); Within = zeros(reg,sub);
    Participate = zeros(reg,sub); Segregate = zeros(reg,sub);
    Betweenness = zeros(reg,sub);
    Local_Efficiency = zeros(reg,sub);
    Nodal_Efficiency = zeros(reg,sub);
    
    for i = 1:sub
        [Between(:,i),Within(:,i),Participate(:,i),Segregate(:,i)] = get_Connectivity(W0(:,:,i));
        Nodal_Efficiency(:,i) = efficiency_nodal(W0(:,:,i));
        Local_Efficiency(:,i) = efficiency_wei(W0(:,:,i),2);
        Betweenness(:,i) = betweenness_wei(W0(:,:,i));
    end
    Degree = reshape(degrees_und(W0),reg,sub);
    Strength = reshape(strengths_und(W0),reg,sub);
    
    % Create each column for the data table.
    if k == 1
        W = W0;
        Age = Age0;
        Age_Partitioned = Age0_Partitioned;
        Research_Group = ones(size(W0,3),1).*codeNums(k);
        Class = Class0;
        rho_G = rho_G0;
        Chimera_Index = Chimera_Index0;
        Metastability_Index = Metastability_Index0;
        Degrees = Degree;
        Strengths = Strength;
        Betweens = Between; Withins = Within;
        Participates = Participate; Segregates = Segregate;
        Nodal_Efficiencies = Nodal_Efficiency; 
        Local_Efficiencies = Local_Efficiency;
        Betweennesses = Betweenness;
        X=X0;
        Y=Y0;
        Z=Z0;
        
    else
        W = cat(3,W,W0);
        Age = cat(1,Age,Age0);
        Age_Partitioned = cat(1,Age_Partitioned,Age0_Partitioned);
        Research_Group = cat(1,Research_Group,ones(size(W0,3),1).*codeNums(k));
        X = cat(2,X,X0);
        Y = cat(2,Y,Y0);
        Z = cat(2,Z,Z0);
        Class = cat(2,Class,Class0);
        Degrees = cat(2,Degrees,Degree);
        Strengths = cat(2,Strengths,Strength);
        Betweens = cat(2,Betweens,Between); Withins = cat(2,Withins,Within);
        Participates = cat(2,Participates,Participate); Segregates = cat(2,Segregates,Segregate);
        Nodal_Efficiencies = cat(2,Nodal_Efficiencies,Nodal_Efficiency);
        Local_Efficiencies = cat(2,Local_Efficiencies,Local_Efficiency);
        Betweennesses = cat(2,Betweennesses,Betweenness);
      
    end
end

Subject = linspace(1,2018,2018)';
Measures = [Degrees,Strengths,Betweens,Withins,Participates,Segregates,Nodal_Efficiencies,Local_Efficiencies,Betweennesses];

% Create the array for the data table.
Subject = linspace(1,2018,2018)';
T_DEGREE_ = [Subject,Research_Group,Age,Age_Partitioned,Class(1,:)',X(1,:)',Y(1,:)',Z(1,:)',Degrees(1,:)',ones(size(Degrees(1,:)')).*1];
T_STRENGTH_ = [Subject,Research_Group,Age,Age_Partitioned,Class(1,:)',X(1,:)',Y(1,:)',Z(1,:)',Strengths(1,:)',ones(size(Degrees(1,:)')).*1];
T_PARTICIPATION_COEFFICIENT_ = [Subject,Research_Group,Age,Age_Partitioned,Class(1,:)',X(1,:)',Y(1,:)',Z(1,:)',Participates(1,:)',ones(size(Degrees(1,:)')).*1];
T_SYSTEM_SEGREGATION_ = [Subject,Research_Group,Age,Age_Partitioned,Class(1,:)',X(1,:)',Y(1,:)',Z(1,:)',Segregates(1,:)',ones(size(Degrees(1,:)')).*1];
T_NODAL_EFFICIENCY_ = [Subject,Research_Group,Age,Age_Partitioned,Class(1,:)',X(1,:)',Y(1,:)',Z(1,:)',Nodal_Efficiencies(1,:)',ones(size(Degrees(1,:)')).*1];
T_LOCAL_EFFICIENCY_ = [Subject,Research_Group,Age,Age_Partitioned,Class(1,:)',X(1,:)',Y(1,:)',Z(1,:)',Local_Efficiencies(1,:)',ones(size(Degrees(1,:)')).*1];
T_BETWEENNESS_ = [Subject, Research_Group,Age,Age_Partitioned,Class(1,:)',X(1,:)',Y(1,:)',Z(1,:)',Betweennesses(1,:)',ones(size(Degrees(1,:)')).*1];
T_BETWEEN_ = [Subject, Research_Group,Age,Age_Partitioned,Class(1,:)',X(1,:)',Y(1,:)',Z(1,:)',Betweens(1,:)',ones(size(Degrees(1,:)')).*1];
T_WITHIN_ = [Subject, Research_Group,Age,Age_Partitioned,Class(1,:)',X(1,:)',Y(1,:)',Z(1,:)',Withins(1,:)',ones(size(Degrees(1,:)')).*1];
Text_ = [repmat(Cognitive_System(1),size(Degrees(1,:)')),repmat(Lobe(1),size(Degrees(1,:)'))];
for i = 1:size(Degrees,1)
    if i == 1
        T_DEGREE = T_DEGREE_;
        T_STRENGTH = T_STRENGTH_;
        T_PARTICIPATION_COEFFICIENT = T_PARTICIPATION_COEFFICIENT_;
        T_SYSTEM_SEGREGATION = T_SYSTEM_SEGREGATION_;
        T_NODAL_EFFICIENCY = T_NODAL_EFFICIENCY_; 
        T_LOCAL_EFFICIENCY = T_LOCAL_EFFICIENCY_;
        T_BETWEENNESS = T_BETWEENNESS_;
        T_BETWEEN = T_BETWEEN_;
        T_WITHIN = T_WITHIN_;
        Text = Text_;
        continue
    end
    T_DEGREE_ = [Subject,Research_Group,Age,Age_Partitioned,Class(i,:)',X(i,:)',Y(i,:)',Z(i,:)',Degrees(i,:)',ones(size(Degrees(1,:)')).*i];
    T_STRENGTH_ = [Subject,Research_Group,Age,Age_Partitioned,Class(i,:)',X(i,:)',Y(i,:)',Z(i,:)',Strengths(i,:)',ones(size(Degrees(1,:)')).*i];
    T_PARTICIPATION_COEFFICIENT_ = [Subject,Research_Group,Age,Age_Partitioned,Class(i,:)',X(i,:)',Y(i,:)',Z(i,:)',Participates(i,:)',ones(size(Degrees(1,:)')).*i];
    T_SYSTEM_SEGREGATION_ = [Subject,Research_Group,Age,Age_Partitioned,Class(i,:)',X(i,:)',Y(i,:)',Z(i,:)',Segregates(i,:)',ones(size(Degrees(1,:)')).*i];
    T_NODAL_EFFICIENCY_ = [Subject,Research_Group,Age,Age_Partitioned,Class(i,:)',X(i,:)',Y(i,:)',Z(i,:)',Nodal_Efficiencies(i,:)',ones(size(Degrees(1,:)')).*i];
    T_LOCAL_EFFICIENCY_ = [Subject,Research_Group,Age,Age_Partitioned,Class(i,:)',X(i,:)',Y(i,:)',Z(i,:)',Local_Efficiencies(i,:)',ones(size(Degrees(1,:)')).*i];
    T_BETWEENNESS_ = [Subject, Research_Group,Age,Age_Partitioned,Class(i,:)',X(i,:)',Y(i,:)',Z(i,:)',Betweennesses(i,:)',ones(size(Degrees(1,:)')).*i];
    T_BETWEEN_ = [Subject, Research_Group,Age,Age_Partitioned,Class(i,:)',X(i,:)',Y(i,:)',Z(i,:)',Betweens(i,:)',ones(size(Degrees(1,:)')).*i];
    T_WITHIN_ = [Subject, Research_Group,Age,Age_Partitioned,Class(i,:)',X(i,:)',Y(i,:)',Z(i,:)',Withins(i,:)',ones(size(Degrees(1,:)')).*i];
    Text_ = [repmat(Cognitive_System(i),size(Degrees(1,:)')),repmat(Lobe(i),size(Degrees(1,:)'))];

    T_DEGREE = [T_DEGREE;T_DEGREE_];
    T_STRENGTH = [T_STRENGTH; T_STRENGTH_];
    T_PARTICIPATION_COEFFICIENT = [T_PARTICIPATION_COEFFICIENT; T_PARTICIPATION_COEFFICIENT_];
    T_SYSTEM_SEGREGATION = [T_SYSTEM_SEGREGATION; T_SYSTEM_SEGREGATION_];
    T_NODAL_EFFICIENCY = [T_NODAL_EFFICIENCY; T_NODAL_EFFICIENCY_];
    T_LOCAL_EFFICIENCY = [T_LOCAL_EFFICIENCY;T_LOCAL_EFFICIENCY_];
    T_BETWEENNESS = [T_BETWEENNESS;T_BETWEENNESS_];
    T_BETWEEN = [T_BETWEEN;T_BETWEEN_];
    T_WITHIN = [T_WITHIN;T_WITHIN_];
    Text = [Text;Text_];
end

% Create the data table.
Table_Str = array2table(Text,'VariableNames',{'Cognitive_System','Lobe'});
Measures = ["Degree","Strength","Participation_Coefficient","System_Segregation","Nodal_Efficiency","Local_Efficiency","Betweenness","Within","Between"];
TABLES = [T_DEGREE, T_STRENGTH, T_PARTICIPATION_COEFFICIENT, T_SYSTEM_SEGREGATION,T_NODAL_EFFICIENCY,T_LOCAL_EFFICIENCY,T_BETWEENNESS,T_WITHIN,T_BETWEEN];
TABLES = reshape(TABLES,size(TABLES,1),size(T_DEGREE,2),12);
for i = 1:length(Measures)
    measure = Measures(i);
    Table_Int = array2table(TABLES(:,:,i),'VariableNames',{'Subject','Research_Group','Age','Age_Partitioned','Classification','X','Y','Z','Result','Brain_Region'});
    Table_ = [Table_Int,Table_Str];
    writetable(Table_,main_path+"/Data/"+measure+".txt",'Delimiter',' ');
end

%%
function Community_assignment = get_Network_Communities(G,Runs)
    for i = 1:Runs
        [Se,~] = community_louvain(G,2);
        Community_Numbers(:,i) = Se;
    end
    N = size(G,1);
    Agreement_Matrix = zeros(N,N);
    for i = 1:N
        for j = 1:N
            Agreement_Matrix(i,j) = mean(Community_Numbers(i,:) == Community_Numbers(j,:));
        end
    end

    [Community_assignment,~] = community_louvain(Agreement_Matrix);
end

function [Between,Within,PC,Segregation] = get_Connectivity(G)
    gamma = 1; % Resolution parameter
    Se = get_Network_Communities(G,100);

    Within = zeros(size(G,1),1);
    Total = sum(G)';
    Community_number = unique(Se);
    for i = 1:length(Community_number)
        idx = (Se == Community_number(i));
        Within(idx) = sum(G(idx,idx))';
    end
    Between = Total - Within;

    Segregation = (Within - Between)./Within;
    PC = participation_coef(G,Se);
end

function E_nodal = efficiency_nodal(W)
    n = length(W); % Number of nodes.
    dij = distance_wei(1./W); % shortest distances.
    D = 1./dij; % inverse distance.
    D(1:n+1:end) = 0; % skip all diagonals.
    E_nodal = sum(D)./(n-1); % Get nodal efficiency.
    E_nodal = E_nodal'; % Inverse for set-up.
end
