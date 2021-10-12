function resultsEx = performExperiments(requirements, runMonteCarlo)
% performExperiments solves rigid, linear, nonlinear, and robust SEAs.
%   C = performExperiments(req, runMC) evaluates the optimal solution over
%   the possible values of uncertainty if runMC = true.

if (isfield(requirements, 'unc'))
    tauUnc          = requirements.unc.tauUnc;
    tauSUnc_a       = requirements.unc.tauSUnc_a;
    tauSUnc_m       = requirements.unc.tauSUnc_m;
    qldUnc_a        = requirements.unc.qldUnc_a;
    qldUnc_m        = requirements.unc.qldUnc_m;
    qlddUnc_a       = requirements.unc.qlddUnc_a;
    qlddUnc_m       = requirements.unc.qlddUnc_m;
    comUncNonLinear = requirements.unc.comUncNonLinear;
end

%-Define the load trajectory and SEA motor and spring parameters
[PEAParam, trajectory] = load_PEA_Parameters(requirements);

%-SOLVE FOR THE NOMINAL CASE (NO UNCERTAINTY)
%-Define the motor requirements with no spring and no uncertainty
robotRigid    = fCostRigidCase(PEAParam, trajectory);
%-Define the motor requirements with linear spring and no uncertainty
linearPEA     = optimal_PEA_Spring(PEAParam, trajectory, 'linear');
%-Define motor requirements with nonlinear spring and no uncertainty
nonlinearPEA  = optimal_PEA_Spring(PEAParam, trajectory, 'nonlinear');

%-SOLVE FOR THE UNCERTAIN CASE.
if (isfield(requirements, 'unc'))
    % Uncertainty in torque (unmodeled dynam.)
    PEAParam.tauUnc = PEAParam.RatedTorque*tauUnc;
    % Additive and multiplicative uncertainty in torque (unmodeled dynam.)
    trajectory.tauSUnc_a = rms(trajectory.torque)*tauSUnc_a;
    trajectory.tauSUnc_m = tauSUnc_m;
    % Additive and multiplicative uncertainty in velocity
    trajectory.qldUnc_a = rms(trajectory.qld)*qldUnc_a;
    trajectory.qldUnc_m = qldUnc_m;
    % Additive and multiplicative uncertainty in acceleration
    trajectory.qlddUnc_a = rms(trajectory.qldd)*qlddUnc_a;
    trajectory.qlddUnc_m = qlddUnc_m;
    % Percentage of uncertainty stiffness (1+-comp)
    PEAParam.comUncNonLinear = linearPEA.stiffness(1)*comUncNonLinear;
    robustLinearPEA = optimal_PEA_Spring(PEAParam, trajectory, 'linear');
    robustNonlinearPEA  = optimal_PEA_Spring(PEAParam, trajectory, 'nonlinear');
    resultsEx.robustlinearPEA       = robustLinearPEA;
    resultsEx.robustNonlinearPEA    = robustNonlinearPEA;
end



resultsEx.PEAParam              = PEAParam;
resultsEx.trajectory            = trajectory;
resultsEx.robotRigid            = robotRigid;
resultsEx.linearPEA             = linearPEA;
resultsEx.nonlinearPEA          = nonlinearPEA;
resultsEx.requirements          = requirements;

if (runMonteCarlo)
    %------------------ MONTE CARLO SIMULATION ----------------------------
    N_trials  = requirements.MCTrials;%Number of trials for the Monte Carlo
    for i = 1:N_trials
        %- Variables to be modified during random experiment
        MCNonlinearRobust   = robustNonlinearPEA;
        MCLinearRobust      = robustLinearPEA;
        MCNonlinearNomina   = nonlinearPEA;
        MCLinearNomina      = linearPEA;
        MCRobot             = PEAParam;
        MCtrajectory        = trajectory;
        n                   = length(trajectory.torque);
        
        %-Sample random variables to implement uncertainty
        %-Uncertainty in the manufacturing of the spring
        MCNonlinearRobust.stiffness = MCNonlinearRobust.stiffness + ...
            (rand(1)*2-1)*PEAParam.comUncNonLinear;
        MCLinearRobust.stiffness    = MCLinearRobust.stiffness + ...
            (rand(1)*2-1)*PEAParam.comUncNonLinear;
        MCNonlinearNomina.stiffness   = MCNonlinearNomina.stiffness + ...
            (rand(1)*2-1)*PEAParam.comUncNonLinear;
        MCLinearNomina.stiffness   = MCLinearNomina.stiffness + ...
            (rand(1)*2-1)*PEAParam.comUncNonLinear;
        %-Uncertainty in the load velocity
        MCtrajectory.qld    = trajectory.qld.*(1+(rand(n,1)*2-1)*...
            trajectory.qldUnc_m) + (rand(n,1)*2-1)*trajectory.qldUnc_a;
        %-Uncertainty in the load acceleration
        MCtrajectory.qldd   = trajectory.qldd*(1+(rand(1)*2-1)*...
            trajectory.qlddUnc_m)...
            + (rand(n,1)*2-1)*trajectory.qlddUnc_a;
        %-Uncertainty in the load torque
        MCtrajectory.torque =trajectory.torque.*(1+(rand(n,1)*2-1)*...
            trajectory.tauSUnc_m) + (rand(n,1)*2-1)*trajectory.tauSUnc_a;
        %-Uncertainty in the motor dynamics
        torqueUncertain     = (rand(n,1)*2-1)*PEAParam.tauUnc;
        
        % Perform stochastic experiment
        rigidDel.stiffness = zeros(n, 1);
        trialRigid              = energy_constraints_nonlinear_spring(...
            MCRobot, MCtrajectory, rigidDel, torqueUncertain);
        
        trialLinearNominal      = energy_constraints_nonlinear_spring(...
            MCRobot, MCtrajectory, MCLinearNomina, torqueUncertain);
        
        trialNonlinearNominal    = energy_constraints_nonlinear_spring(...
            MCRobot, MCtrajectory, MCNonlinearNomina, torqueUncertain);
        
        trialLinearRobust        = energy_constraints_nonlinear_spring(...
            MCRobot, MCtrajectory, MCLinearRobust, torqueUncertain);
        
        trialNonlinearRobust     = energy_constraints_nonlinear_spring(...
            MCRobot, MCtrajectory, MCNonlinearRobust, torqueUncertain);
        
        %-Storing the experiments
        RigidExperiment(i)              = trialRigid;
        LinearNominalExperiment(i)      = trialLinearNominal;
        NonlinearNominalExperiment(i)   = trialNonlinearNominal;
        LinearRobustExperiment(i)       = trialLinearRobust;
        NonlinearRobustExperiment(i)    = trialNonlinearRobust;
        RobotExperiment(i)              = MCRobot;
        TrajectoryExperiment(i)         = MCtrajectory;
    end
    % Checking statistics
    %     nRigidFeas           = zeros(N_trials,1);
    %     nLinearNominalFeas   = zeros(N_trials,1);
    %     nNonlinearNominalFeas= zeros(N_trials,1);
    %     nLinearRobustFeas    = zeros(N_trials,1);
    %     nNonlinearRobustFeas = zeros(N_trials,1);
    for i = 1:N_trials
        nRigidFeas(i)           = RigidExperiment(i).feasible;
        nLinearNominalFeas(i)   = LinearNominalExperiment(i).feasible;
        nNonlinearNominalFeas(i)= NonlinearNominalExperiment(i).feasible;
        nLinearRobustFeas(i)    = LinearRobustExperiment(i).feasible;
        nNonlinearRobustFeas(i) = NonlinearRobustExperiment(i).feasible;
        
        costRigidFeas(i)           = RigidExperiment(i).objFun;
        costLinearNominalFeas(i)   = LinearNominalExperiment(i).objFun;
        costNonlinearNominalFeas(i)= NonlinearNominalExperiment(i).objFun;
        costLinearRobustFeas(i)    = LinearRobustExperiment(i).objFun;
        costNonlinearRobustFeas(i) = NonlinearRobustExperiment(i).objFun;
    end
    %Taking advantage that the entries to the function are logical values.
    percentFeas = @(x) numel(x(x))/N_trials;
    resultsEx.rigidMC.percFeas          = percentFeas(nRigidFeas);
    resultsEx.linearNomMC.percFeas      = percentFeas(nLinearNominalFeas);
    resultsEx.nonlinearNomMC.percFeas = percentFeas(nNonlinearNominalFeas);
    resultsEx.linearRobustMC.percFeas   = percentFeas(nLinearRobustFeas);
    resultsEx.nonlinearRobustMC.percFeas=percentFeas(nNonlinearRobustFeas);
    
    resultsEx.rigidMC.cost          = costRigidFeas;
    resultsEx.linearNomMC.cost      = costLinearNominalFeas;
    resultsEx.nonlinearNomMC.cost   = costNonlinearNominalFeas;
    resultsEx.linearRobustMC.cost   = costLinearRobustFeas;
    resultsEx.nonlinearRobustMC.cost= costNonlinearRobustFeas;
    
    resultsEx.rigidMC.costRela      = (costRigidFeas-costRigidFeas)./...
        costRigidFeas*100;
    resultsEx.linearNomMC.costRela  = (costLinearNominalFeas - ...
        costRigidFeas)./costRigidFeas*100;
    resultsEx.nonlinearNomMC.costRela = (costNonlinearNominalFeas -...
        costRigidFeas)./costRigidFeas*100;
    resultsEx.linearRobustMC.costRela = (costLinearRobustFeas - ...
        costRigidFeas)./costRigidFeas*100;
    resultsEx.nonlinearRobustMC.costRela = (costNonlinearRobustFeas - ...
        costRigidFeas)./costRigidFeas*100;
    
    resultsEx.rigidMC.experiment            = RigidExperiment;
    resultsEx.linearNomMC.experiment        = LinearNominalExperiment;
    resultsEx.nonlinearNomMC.experiment     = NonlinearNominalExperiment;
    resultsEx.linearRobustMC.experiment     = LinearRobustExperiment;
    resultsEx.nonlinearRobustMC.experiment  = NonlinearRobustExperiment;
    resultsEx.robotMC.experiment            = RobotExperiment;
    resultsEx.trajExMC                      = TrajectoryExperiment;
end
end