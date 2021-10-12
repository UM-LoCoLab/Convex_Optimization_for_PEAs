function [feasible, satisfiedConstraint, pST] = fIsMotorFeasible(robot, ...
    results, titleString)
% Loading motor parameters
Vin     = robot.voltage;
kt      = robot.kt;
R       = robot.R;
% Loading motor trajectories
qmd     = results.qmd;
tauM    = results.tauM;
% Increase 2% the right hand side to check constraints satisfaction
marign2MakeProbFeas = 1.02;
% Are motor constraints satisfied?
satisfiedConstraint.SpeedTorque = all(abs(tauM)<=...
    (kt/R*Vin-kt^2/R*abs(qmd))*marign2MakeProbFeas );
satisfiedConstraint.Torque      = all(abs(tauM)<= ...
    (robot.peakTorque)*marign2MakeProbFeas);
% Are all motor constraints satisfied?
feasible = all(...
    [satisfiedConstraint.SpeedTorque,satisfiedConstraint.Torque]);

vMax            = Vin/kt;
vAtMaxTorque    = -(robot.peakTorque - kt/R*Vin)*R/(kt^2);
maxTorque       = robot.peakTorque;

v = [-vMax 0; -vAtMaxTorque -maxTorque; vAtMaxTorque -maxTorque;...
    vMax 0; vAtMaxTorque maxTorque; -vAtMaxTorque maxTorque];
f = [1, 2, 3, 4, 5, 6];

% Plot speed torque
pST.title = sprintf('%s - Motor Speed-Torque',...
    titleString);
pST.f       = f;
pST.v       = v;

%Example of plotting using the patch approach
% figure, hold on, title(plotSpeedTorque.title)
% patch('Faces',f,'Vertices',v,'FaceColor','red','FaceAlpha',.1)
% plot(qmd, tauM);
% ylabel('Motor Torque [Nm]'), xlabel('Motor Velocity [rad/sec]')
% hold off
end