function [optiRobot] = fPostProcessing(qm, qmd, qmdd, cvx_status, robot,...
    trajectoryA, varargin)
%-Varargin in the case that uncertain torque affects the motor

%% Post processing of the optimization results from CVX

%--Loading robot parameters
iM = robot.iM;
r = robot.r;
bm = robot.bm;
eta = robot.eta;

%--Loading trajectory A
time = trajectoryA.time;
ql = trajectoryA.ql;

%--Torque done by the spring on the load
torqueLoad   = -trajectoryA.torque;
torqueElastic   = trajectoryA.tau_P;
if numel(varargin) > 0
    torqueUncer     = varargin{1};
else
    torqueUncer     = zeros(size(torqueElastic));
end
%-- Post processing
%tauM = iM*qmdd + bm*qmd -((torqueLoad+torqueElastic)./(r*eta)) + torqueUncer;
tauM = trajectoryA.e_vec + trajectoryA.F_matr * trajectoryA.k_p;
tauM(end+1) = tauM(1);
tauM = tauM + torqueUncer;

%-- Elongation. Position of the motor after transmission
%elong = ql-(qm./r);
elong = ql;          % Check
minElongIndx = find(abs(elong) == min(abs(elong)));
elong       = elong - elong(minElongIndx);

%--CREATE POST PROCESSING FUNCTION - Walking
resultsCVX.qm = qm;
resultsCVX.qmd = qmd;
resultsCVX.qmdd = qmdd;
resultsCVX.elong = elong;
resultsCVX.tauM = tauM;
resultsCVX.solverStatus = cvx_status;
resultsCVX.torqueElastic = torqueElastic;

trajectory.ql = ql;
trajectory.qld = trajectoryA.qld;
trajectory.qldd = trajectoryA.qldd;
trajectory.time = time;

%--Getting results
optiRobot = fOptimizationPostprocessing(robot, trajectory, resultsCVX);




end