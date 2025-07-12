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
Cingulate = [110,52,16,73,45,104,26,84];
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
    W0(isnan(W0)) = 0;
    sub = size(W0,3); % The number of individuals.
    reg = size(W0,1); % The number of brain regions.
    Class0 = Data.Classification; % The synchrony class. 1>Async,2>Chimera,3>Sync.
    P = Data.Binarized_Patterns;
    
    Age0 = Data.Age;    
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

    Res_Att0 = reshape(P(1,:,:),reg,sub);
    Res_Aud0 = reshape(P(2,:,:),reg,sub);
    Res_CO0 = reshape(P(3,:,:),reg,sub);
    Res_FP0 = reshape(P(4,:,:),reg,sub);
    Res_DM0 = reshape(P(5,:,:),reg,sub);
    Res_MS0 = reshape(P(6,:,:),reg,sub);
    Res_SC0 = reshape(P(7,:,:),reg,sub);
    Res_VT0 = reshape(P(8,:,:),reg,sub);
    Res_V0 = reshape(P(9,:,:),reg,sub);

    Degree = reshape(degrees_und(W0),reg,sub);
    if k == 1
        Age = Age0;
        Age_Partitioned = Age0_Partitioned;
        X = X0;
        Y = Y0;
        Z = Z0;
        Research_Group = ones(size(W0,3),1).*codeNums(k);
        Class = Class0;
        Degrees = Degree;

        Res_Att = Res_Att0;
        Res_Aud = Res_Aud0;
        Res_CO = Res_CO0;
        Res_FP = Res_FP0;
        Res_DM = Res_DM0;
        Res_MS = Res_MS0;
        Res_SC = Res_SC0;
        Res_VT = Res_VT0;
        Res_V = Res_V0;

    else
        Age = cat(1,Age,Age0);
        Age_Partitioned = cat(1,Age_Partitioned,Age0_Partitioned);
        Research_Group = cat(1,Research_Group,ones(size(W0,3),1).*codeNums(k));
        X = cat(2,X,X0);
        Y = cat(2,Y,Y0);
        Z = cat(2,Z,Z0);
        Class = cat(2,Class,Class0);
        Degrees = cat(2,Degrees,Degree);
        
        Res_Att = cat(2,Res_Att,Res_Att0);
        Res_Aud = cat(2,Res_Aud,Res_Aud0);
        Res_CO = cat(2,Res_CO,Res_CO0);
        Res_FP = cat(2,Res_FP,Res_FP0);
        Res_DM = cat(2,Res_DM,Res_DM0);
        Res_MS = cat(2,Res_MS,Res_MS0);
        Res_SC = cat(2,Res_SC,Res_SC0);
        Res_VT = cat(2,Res_VT,Res_VT0);
        Res_V = cat(2,Res_V,Res_V0);
    end

end

% Create the array for the data table.
Subject = linspace(1,2018,2018)';

T_SYNCHRONY_CS_PATTERNS_ = [Subject,Research_Group,Age,Age_Partitioned,Class(1,:)',X(1,:)',Y(1,:)',Z(1,:)',Res_Att(1,:)',Res_Aud(1,:)',Res_CO(1,:)',Res_FP(1,:)',Res_DM(1,:)',Res_MS(1,:)',Res_SC(1,:)',Res_VT(1,:)',Res_V(1,:)',ones(size(Degrees(1,:)')).*1];
Text_ = [repmat(Cognitive_System(1),size(Degrees(1,:)')),repmat(Lobe(1),size(Degrees(1,:)'))];
for i = 1:size(Degrees,1)
    if i == 1
        
        T_SYNCHRONY_CS_PATTERNS = T_SYNCHRONY_CS_PATTERNS_; 
        Text = Text_;
        continue
    end
    
    T_SYNCHRONY_CS_PATTERNS_ = [Subject,Research_Group,Age,Age_Partitioned,Class(i,:)',X(i,:)',Y(i,:)',Z(i,:)',Res_Att(i,:)',Res_Aud(i,:)',Res_CO(i,:)',Res_FP(i,:)',Res_DM(i,:)',Res_MS(i,:)',Res_SC(i,:)',Res_VT(i,:)',Res_V(i,:)',ones(size(Degrees(1,:)')).*i];
    Text_ = [repmat(Cognitive_System(i),size(Degrees(1,:)')),repmat(Lobe(i),size(Degrees(1,:)'))];

    T_SYNCHRONY_CS_PATTERNS = [T_SYNCHRONY_CS_PATTERNS; T_SYNCHRONY_CS_PATTERNS_];
    Text = [Text;Text_];
end

% Create the data table.
Table_Str = array2table(Text,'VariableNames',{'Cognitive_System','Lobe'});

Pattern_Table_Int = array2table(T_SYNCHRONY_CS_PATTERNS,'VariableNames',{'Subject','Research_Group','Age','Age_Partitioned','Classification','X','Y','Z','Att','Aud','CO','FP','DM','MS','SC','VT','V','Brain_Region'});
Pattern_Table = [Pattern_Table_Int,Table_Str];
writetable(Pattern_Table,main_path+"/Data/"+"Patterns.txt",'Delimiter',' ');


