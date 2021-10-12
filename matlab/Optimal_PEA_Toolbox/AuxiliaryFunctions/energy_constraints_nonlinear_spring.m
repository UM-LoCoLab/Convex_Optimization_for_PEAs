function result = energy_constraints_nonlinear_spring(robot, ...
    trajectory, inSolution, varargin)
%ENERGY_CONSUMPTION_LINEAR_SPRING Summary of this function goes here
%   Varargin for the case that uncertain motor torque acts on the motor.

comp    = inSolution.compliance;

%-- Loading trajectory
ql      = trajectory.ql(1:end-1);
qld     = trajectory.qld(1:end-1);
qldd    = trajectory.qldd(1:end-1);
torS    = -trajectory.torque(1:end-1);
torSd   = -trajectory.torqued(1:end-1);
% torSdd  = -trajectory.torquedd(1:end-1);
comp    = comp(1:end-1);
time    = trajectory.time(1:end-1);

%-Index when torque is the closest to zero
minTorqueIndx = find(abs(torS) == min(abs(torS)));


%-Numerical difference
D1 = fnumerical_derivative_matrix_FixedSampleRate(time(2)-time(1), ...
    length(comp), 1);

%-- Derivative of elongation
elon_d  = comp.*torSd;
elon_dd = D1*elon_d;
elon    = cumtrapz(time, elon_d);
elon    = elon - elon(minTorqueIndx);

%--Loading robot parameters
r = robot.r;

qm = (ql - elon)*r;
qmd = (qld - elon_d)*r;
qmdd = (qldd - elon_dd)*r;

qm(end+1)   = qm(1);
qmd(end+1)  = qmd(1);
qmdd(end+1) = qmdd(1);

result = ...
fPostProcessing(qm, qmd, qmdd, 'none', robot, trajectory, varargin{1});
% result = fPostProcessing(qm, qmd, qmdd, 'none', robot, trajectory);

[feasible, ~, ~] = fIsMotorFeasible(robot, result, 'Test');

result.feasible = feasible;

% Evaluate cost function
%-Cost function: energy, 2NormVelocity, 2NormTorque, and 2NormElongation.
if strcmp(trajectory.objFun, 'energy')
    objFun = result.energy_dissipated;
elseif strcmp(trajectory.objFun, 'RMS_Velocity')
    objFun = norm(result.qmd)/sqrt(length(result.qmd));
elseif strcmp(trajectory.objFun, 'RMS_Torque')
    objFun = norm(result.tauM)/sqrt(length(result.qmd));
elseif strcmp(trajectory.objFun, 'RMS_Elongation')
    objFun = norm(result.elong)/sqrt(length(result.qmd));    
else
    error('Please select adequate value of objective function')
end

result.objFun = objFun;

end

