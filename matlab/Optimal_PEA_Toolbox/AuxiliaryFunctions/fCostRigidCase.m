function [robotRigid] = fCostRigidCase(robot, trajectory)
%-- FUNCTION TO CALCULATE MOTOR REQUIREMENTS FOR THE GIVEN TRAJECTORY(ES)
%-- If varargin{1} is not empty consider the case of multiple trajectories

%-- Analysis in only one trajectory
resultsRigid = locFunPreprocessing(robot, trajectory);
%-- Post processing results
robotRigid = fOptimizationPostprocessing(robot, trajectory, resultsRigid);
%-- Computing optimization objective
if(strcmp(trajectory.objFun, 'RMS Velocity'))
    costCVX = rms(robotRigid.qmd);
elseif(strcmp(trajectory.objFun, 'RMS Torque'))
    costCVX = rms(robotRigid.tauM);
elseif(strcmp(trajectory.objFun, 'energy'))
    costCVX = robotRigid.energy_total;
elseif  (strcmp(trajectory.objFun, 'RMS Elongation'))
    costCVX = rms(robotRigid.elong);
end
robotRigid.cost = costCVX;
%-- Checking constraints
[feasible, violatedConstraint, SpeedTorquePlot] = fIsMotorFeasible(...
    robot, resultsRigid, 'Rigid Actuator');
robotRigid.feasible = feasible;
robotRigid.violatedConstraint = violatedConstraint;
robotRigid.plot = SpeedTorquePlot;
end

function [resultsRigid] = locFunPreprocessing(robot, trajectory)
%--------------Definitions
iM = robot.iM;
r = robot.r;
bm = robot.bm;
% km = robot.km;
eta = robot.eta;
%-------------Loading the reference trajectory
% time = trajectory.time;
ql = trajectory.ql;
qld = trajectory.qld;
qldd = trajectory.qldd;
torque = trajectory.torque;
%------------Calculating the motion of the motor
qm = ql.*r;
qmd = qld.*r;
qmdd = qldd.*r;
tauM = iM.*qmdd+ bm.*qmd+ torque./(eta*r);
%--CREATE POST PROCESSING FUNCTION
resultsRigid.qm = qm;
resultsRigid.qmd = qmd;
resultsRigid.qmdd = qmdd;
resultsRigid.tauM = tauM;
resultsRigid.torqueElastic = 0;
resultsRigid.elong = 0;
resultsRigid.solverStatus = 'none';

end