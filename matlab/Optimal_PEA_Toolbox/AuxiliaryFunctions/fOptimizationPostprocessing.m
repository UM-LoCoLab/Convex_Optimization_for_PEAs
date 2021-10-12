function [optiRobot] = fOptimizationPostprocessing(robot, trajectory, resultsOpt)
%---------------------------------------------
%[optiRobot] = fOptimizationPostprocessing(robot, trajectory, resultsOpt)
%---------------------------------------------
%Function to analize and organize the results from the optimization solvers

ql = trajectory.ql;
qld = trajectory.qld;
qldd = trajectory.qldd;
r = robot.r;
bm = robot.bm;
kt = robot.kt;
R = robot.R;
km = robot.km;
time = trajectory.time;
qm = resultsOpt.qm;
qmd = resultsOpt.qmd;
qmdd = resultsOpt.qmdd;
elong = resultsOpt.elong;
tauM = resultsOpt.tauM;
cvx_status = resultsOpt.solverStatus;
torqueElastic = resultsOpt.torqueElastic;

% Calculating new expressions from the optimization
qmAT = qm./r;
qmdAT = qmd./r;
taumAT = tauM.*r;
power = tauM.*qmd;

%----------------Checking individual terms
% inertia = trapz(time,iM*qmdd.*qmd);
% viscousF = trapz(time,bm*qmd.*qmd);
% externalLoad = trapz(time,(torqueElastic./r).*qmd);
% externalLoad2 = trapz(qm,(torqueElastic./r));
% energyInTheEE = trapz(elong.',torqueElastic);
% TotalEnergy = inertia+viscousF-externalLoad;

optiRobot.energy_tauMqmd = trapz(time,tauM.*qmd);
optiRobot.energy_windingHeat = trapz(time,tauM.^2)/km^2;
optiRobot.energy_viscousFriction = trapz(time,qmd.^2) *bm;
optiRobot.energy_total = optiRobot.energy_tauMqmd+ optiRobot.energy_windingHeat;
optiRobot.energy_dissipated = optiRobot.energy_windingHeat + ...
    optiRobot.energy_viscousFriction;
optiRobot.power = power;
optiRobot.peakPower = norm(power, Inf);

optiRobot.solverStatus = cvx_status;
optiRobot.qmAT = qmAT;
optiRobot.qmdAT = qmdAT;
optiRobot.taumAT = taumAT;
optiRobot.tauM = tauM;
optiRobot.qm = qm;
optiRobot.qmd = qmd;
optiRobot.qmdd = qmdd;
optiRobot.time = time;
optiRobot.ql = ql;
optiRobot.qld = qld;
optiRobot.qldd = qldd;
optiRobot.tauElastic = torqueElastic;
optiRobot.elong = elong;

end

