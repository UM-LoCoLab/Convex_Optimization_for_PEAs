function [ result ] = energy_consumption_linear_spring(robot, trajectory, k )
%ENERGY_CONSUMPTION_LINEAR_SPRING Summary of this function goes here
%   Detailed explanation goes here

%--Loading robot parameters
r = robot.r;
alpha1 = 1/k;

%-- Loading trajectory
ql = trajectory.ql;
qld = trajectory.qld;
qldd = trajectory.qldd;
torEla = -trajectory.torque;
torElad = -trajectory.torqued;
torEladd = -trajectory.torquedd;

qm = (ql - torEla*alpha1)*r;
qmd = (qld - torElad*alpha1)*r;
qmdd = (qldd - torEladd*alpha1)*r;

result = fPostProcessing(qm, qmd, qmdd, 'none', robot, trajectory);
result.k = k;

end

