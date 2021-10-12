function [d, e] = robust_constraints(trajectory, robot)
%-- Loading constraints for each of the trajectories. Each trajectory is a
%column of the matrix of constraints

%--Loading robot parameters
iM = robot.iM;
r = robot.r;
bm = robot.bm;
eta = robot.eta;
vIn = robot.voltage;
R = robot.R;
kt = robot.kt;
% km = robot.km;
% time = trajectory.time;

mB = trajectory.mBar;
mUnc = trajectory.tauSUnc_m;
maxElong = robot.maxElong;

qldUnc = trajectory.qldUnc_a;
qlddUnc = trajectory.qlddUnc_a;
tauMax = robot.peakTorque;
qmdMax = robot.maxVelo;
etaUnc = robot.etaUnc;
tauUncLow = -robot.tauUnc;
% Uncertainty in compliance
comUnc = robot.comUnc;

ql = trajectory.ql;
qld = trajectory.qld;
qldd = trajectory.qldd;
torqNorma = trajectory.torque./mB;
torqNormad = trajectory.torqued./mB;
torqNormadd = trajectory.torquedd./mB;
torqueEla = trajectory.torque;
torqueElad = trajectory.torqued;
torqueEladd = trajectory.torquedd;

%-- Defining parameters
gamma1 = -(iM*torqueEladd*r + bm*torqueElad*r);
gamma2 = iM*qldd*r + bm*qld*r - torqueEla/(eta*r);
gamma1Norma = -(iM*torqNormadd*r + bm*torqNormad*r);

%-- Elongation Constraint
d1 = [torqNorma; -torqNorma];
e1 = ones(length(d1), 1)*maxElong/(mB+mUnc);

%-- Torque Constraint
d2 = [gamma1Norma; -gamma1Norma];
fParcial = [-iM*r*(qldd+qlddUnc)-bm*r*(qld+qldUnc)+tauMax+tauUncLow;...
    iM*r*(qldd-qlddUnc)+bm*r*(qld-qldUnc)+tauMax+tauUncLow];
mUncVec = (mB+sign(fParcial)*mUnc);
torqNormaPartial = [torqNorma; -torqNorma];
torqFactor = ((eta+sign(torqNormaPartial)*etaUnc)*r);
e2 = torqNormaPartial./torqFactor+fParcial./mUncVec;

%-- Torque and velocity constraints
%-- A
d3a = [gamma1Norma - torqNormad*kt^2*r/R];
e3aParcial = (vIn*kt/R - iM*r*(qldd+qlddUnc) - bm*r*(qld+qldUnc)+...
    +tauUncLow - kt^2*r/R.*(qld+qldUnc));
uncNume = ((eta+sign(torqNorma)*etaUnc)*r);
e3a = torqNorma./uncNume + e3aParcial./(mB + sign(e3aParcial)*mUnc);
%-- B
d3b = [gamma1Norma + kt^2*r/R.*torqNormad];
e3bParcial = (-iM*(qldd+qlddUnc)*r - bm*(qld+qldUnc)*r+...
    +tauUncLow + vIn*kt/R + kt^2*r/R.*(qld-qldUnc));
uncNume = ((eta+sign(torqNorma)*etaUnc)*r);
e3b = torqNorma./uncNume + e3bParcial./(mB + sign(e3bParcial)*mUnc);
%-- C
d3c = [-gamma1Norma - kt^2*r/R.*torqNormad];
e3cParcial = (+iM*(qldd-qlddUnc)*r + bm*(qld-qldUnc)*r+...
    +tauUncLow + vIn*kt/R - kt^2*r/R.*(qld+qldUnc));
uncNume = ((eta+sign(-torqNorma)*etaUnc)*r);
e3c = -torqNorma./uncNume + e3cParcial./(mB + sign(e3cParcial)*mUnc);
%-- D
d3d = [-gamma1Norma + kt^2*r/R.*torqNormad];
e3dParcial = (+iM*(qldd-qlddUnc)*r + bm*(qld-qldUnc)*r+...
    +tauUncLow + vIn*kt/R + kt^2*r/R.*(qld-qldUnc));
uncNume = ((eta+sign(-torqNorma)*etaUnc)*r);
e3d = -torqNorma./uncNume + e3dParcial./(mB + sign(e3dParcial)*mUnc);

d = [d1; d2; d3a; d3b; d3c; d3d];
% d = [d1; d2; d3a];
%-- Uncertainty due to stiffness implementation
% d = d + abs(d)*kUnc;

d1 = d + d*comUnc;
d2 = d - d*comUnc;
d = max(d1, d2);
% dup = zeros(length(d));
% for i = 1:length(d)
%     d1 = d + d*kUnc;
%     d2 = d - d*kUnc;
%     dup(i) = max(d1, d2);
% end

e = [e1; e2; e3a; e3b; e3c; e3d];
% e = [e1; e2; e3a];
end