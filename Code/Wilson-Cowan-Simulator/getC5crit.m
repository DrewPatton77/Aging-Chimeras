function getC5crit()
% =======================================================================
% Running this evaluates where the critical point is for the parameter 
%  c5 (or c_E) via a midpoint method.
%
% Each research group's individual's structural connectome has a c5 that is
%  evaluated for here.
% =======================================================================


setCodes = ["355581","487062","602004","634143","750313","807122"];

main_path = cd;
main_path = strsplit(main_path,'\\Code');
main_path = string(main_path{1});

v = 10; % Conductance velocity along an mylenated axon in m/s.
for k = 1:6
    Data = importdata(pwd+"/Data/"+setCodes(k)+"_rawData_healthy"+".mat"); % All the data for the respective dataset.
    path = pwd+"/Data/"+setCodes(k)+"_Structural_Connectome.mat";
    Data.Time_Delays = Data.Streamline_Length/v; % Retrieve the time delays between each brain region.
    get_c5crit(Data,path)
end

quit();

% Find critical c5 values for each subject in the dataset.
function get_c5crit(Data,path)
    
    % Constants.
    time = 300; % Where
    dt = 1e-3;
    v = 10;
    
    % Initial Guess.
    high_est = 80;
    low_est = 0;
    tol = 0.1;
    maxNum = 50;
    
    c5_crit = zeros(size(Data.Streamline_Density,3),1);
    for subject = 1:size(Data.Streamline_Density,3)
        G = Data.Streamline_Density(:,:,subject); % Structural Connectome.
        D = Data.Streamline_Length(:,:,subject)/v; % Time Delay.
        
        c5_crit(subject) = mid_point(low_est,high_est,tol,maxNum,G,D,time,dt);
    end

    % Save Data.
    Data.c5_crit = c5_crit;
    save(path,Data);


function guess = mid_point(low_est,high_est,tol,maxNum,G,D,time,dt)
    % A mid-point method in finding the critical c5 value for the
    % Wilson-Cowan Model.
    %
    % Inputs:
    % -------
    %         low_est: The lower estimate.
    %        high_est: The higher estimate.
    %             tol: The tolerance.
    %          maxNum: Maximum number of bisections.
    %               G: Connectivity matrix.
    %               D: Time Delay matrix.
    %            time: Simulation duration time.
    %              dt: Time step.
    % Outputs:
    % --------
    %           guess: The c5 critical value
    %                 (if converged in n < maxNum steps).

    n = 1; % Starting value.
    while n < maxNum % Start while loop.

        guess = (low_est+high_est)/2; % The mid-section.
        if high_est - low_est < tol % Tolerance condition.
            guess = high_est; % Opted for the activated state side (one could use the quiet state side).
            return
        end
       
        Dynamics = wc_coupled_stochastic(G,D,time,dt,guess,guess/4,[],[]);
        cond = ( mean(log10(mean(Dynamics.e(200/dt:end,:)))) > -5); % Checking activation.
        
        if cond == 1
            high_est = guess;
        else
            low_est = guess;
        end
    end
    disp('Did not Converge within tolerance.')
    guess = high_est;
