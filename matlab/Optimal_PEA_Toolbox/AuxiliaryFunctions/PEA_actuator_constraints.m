function [Aineq, bineq, L] = PEA_actuator_constraints(...
    motorTransmission, trajectory, cost)

%%%%%%%%%%%%%%%%%%%%%%%%%% Load uncertainty values %%%%%%%%%%%%%%%%%%%%%%%%
if (isfield(trajectory, 'etaUnc'))
    etaUnc  = trajectory.etaUnc;
else
    etaUnc  = 0;
end
%Multiplicative uncertainty in spring torque
if (isfield(trajectory, 'tauSUnc_m'))
    tauSUnc_m  = trajectory.tauSUnc_m;
else
    tauSUnc_m  = 0; %Equivalent to 8.8kg/65.1kg
end

%Additive uncertainty in spring torque
if (isfield(trajectory, 'tauSUnc_a'))
    tauSUnc_a  = trajectory.tauSUnc_a;
else
    tauSUnc_a  = 0; 
end

if (isfield(motorTransmission, 'tauUnc'))
    tauUnc  = motorTransmission.tauUnc;
else
    tauUnc  = 0; 
end

%Additive and multiplicative uncertainty in velocity
if (isfield(trajectory, 'qldUnc_a'))
    qldUnc_a  = trajectory.qldUnc_a;
else
    qldUnc_a  = 0; 
end
if (isfield(trajectory, 'qldUnc_m'))
    qldUnc_m  = trajectory.qldUnc_m;
else
    qldUnc_m  = 0; 
end

%Additive and multiplicative uncertainty in acceleration
if (isfield(trajectory, 'qlddUnc_a'))
    qlddUnc_a  = trajectory.qlddUnc_a;
else
    qlddUnc_a  = 0; 
end
if (isfield(trajectory, 'qlddUnc_m'))
    qlddUnc_m  = trajectory.qlddUnc_m;
else
    qlddUnc_m  = 0; 
end

% Uncertainty compliance (1+-comp)
if (isfield(motorTransmission, 'comUncNonLinear'))
    comUncNonLinear  = motorTransmission.comUncNonLinear;
else
    comUncNonLinear  = 0; 
end

%-- Actuator parameters
kt      = motorTransmission.kt;
R       = motorTransmission.R;        %-Motor resistance [Ohms]
vIn     = motorTransmission.voltage;
eta     = motorTransmission.eta;
r       = motorTransmission.r;
bm      = motorTransmission.bm;
iM      = motorTransmission.iM;

%--Torque done by the spring on the load
torque      = -trajectory.torque(1:end-1);
n           = length(torque);
%Loading kinematics
qld         = trajectory.qld(1:end-1);
qldd        = trajectory.qldd(1:end-1);
%-Creating elongation matrix
L = create_elongation_matrix(trajectory);
% From the cost: qmd = a+B*alpha
a   = cost.a;
B   = cost.B;
% From the cost: tau_m = e + F*alpha
e   = cost.e;
F   = cost.F;

%-------------------------NOMINAL CONSTRAINTS -----------------------------
%-Maximum spring torque constraint
% Spring Torque = L*k_p
%Aineq_mElong         = [L; -L];
%[m_Aineq2, ~]        = size(Aineq_mElong);
%bineq_mElong         = ones(m_Aineq2, 1)*motorTransmission.maxElong;

%-Maximum torque constraint
Aineq_mTroque  = [F; -F];
bineq_mTorque  = [-e+ones(n, 1)*motorTransmission.peakTorque;...
    e+ones(n, 1)*motorTransmission.peakTorque];

%-Torque-Velocity Constraint - 03/26/19 Edgar's notebook.
Aineq_torVelo = [(F+kt^2/R*B);...
    (F-kt^2/R*B);...
    (-F+kt^2/R*B);...
    (-F-kt^2/R*B)];
bineq_torVelo = [kt/R*vIn*ones(n, 1)-e-kt^2/R*a;...
    kt/R*vIn*ones(n, 1)-e+kt^2/R*a;...
    kt/R*vIn*ones(n, 1)+e-kt^2/R*a;...
    kt/R*vIn*ones(n, 1)+e+kt^2/R*a];

%------------------------- ROBUST CONSTRAINTS -----------------------------
% Robust counterpart is all in the LHS in bineq

%-Robust against uncertain torque
bineq_mTorque  = bineq_mTorque - ones(2*n, 1)*abs(tauUnc);
bineq_torVelo  = bineq_torVelo - ones(4*n, 1)*abs(tauUnc);

%-Merging constratins
%Aineq = [Aineq_mElong; Aineq_mTroque; Aineq_torVelo];
%bineq = [bineq_mElong; bineq_mTorque; bineq_torVelo];
Aineq = [Aineq_mTroque; Aineq_torVelo];
bineq = [bineq_mTorque; bineq_torVelo];
%-Robust agains kinematics of the load

% unc_kinematics = [zeros(n,1);...
%     zeros(n,1);...
%     % Start of torque constraints
%     bm*r*(qldUnc_m*abs(qld)+ones(n, 1)*qldUnc_a)+...
%     iM*r*(qlddUnc_m*abs(qldd)+ones(n, 1)*qlddUnc_a);...
%     %
%     bm*r*(qldUnc_m*abs(qld)+ones(n, 1)*qldUnc_a)+...
%     iM*r*(qlddUnc_m*abs(qldd)+ones(n, 1)*qlddUnc_a);...
%     % Start of torque-velocity constraints
%     (bm*r+(kt^2)/R*r)*(qldUnc_m*abs(qld)+ones(n, 1)*qldUnc_a)+...
%     iM*r*(qlddUnc_m*abs(qldd)+ones(n, 1)*qlddUnc_a);...
%     %
%     (bm*r+(kt^2)/R*r)*(qldUnc_m*abs(qld)+ones(n, 1)*qldUnc_a)+...
%     iM*r*(qlddUnc_m*abs(qldd)+ones(n, 1)*qlddUnc_a);...
%     %
%     (bm*r+(kt^2)/R*r)*(qldUnc_m*abs(qld)+ones(n, 1)*qldUnc_a)+...
%     iM*r*(qlddUnc_m*abs(qldd)+ones(n, 1)*qlddUnc_a);...
%     %
%     (bm*r+(kt^2)/R*r)*(qldUnc_m*abs(qld)+ones(n, 1)*qldUnc_a)+...
%     iM*r*(qlddUnc_m*abs(qldd)+ones(n, 1)*qlddUnc_a);...
%     %
%     ];

unc_kinematics = [
    % Start of torque constraints
    bm*r*(qldUnc_m*abs(qld)+ones(n, 1)*qldUnc_a)+...
    iM*r*(qlddUnc_m*abs(qldd)+ones(n, 1)*qlddUnc_a);...
    %
    bm*r*(qldUnc_m*abs(qld)+ones(n, 1)*qldUnc_a)+...
    iM*r*(qlddUnc_m*abs(qldd)+ones(n, 1)*qlddUnc_a);...
    % Start of torque-velocity constraints
    (bm*r+(kt^2)/R*r)*(qldUnc_m*abs(qld)+ones(n, 1)*qldUnc_a)+...
    iM*r*(qlddUnc_m*abs(qldd)+ones(n, 1)*qlddUnc_a);...
    %
    (bm*r+(kt^2)/R*r)*(qldUnc_m*abs(qld)+ones(n, 1)*qldUnc_a)+...
    iM*r*(qlddUnc_m*abs(qldd)+ones(n, 1)*qlddUnc_a);...
    %
    (bm*r+(kt^2)/R*r)*(qldUnc_m*abs(qld)+ones(n, 1)*qldUnc_a)+...
    iM*r*(qlddUnc_m*abs(qldd)+ones(n, 1)*qlddUnc_a);...
    %
    (bm*r+(kt^2)/R*r)*(qldUnc_m*abs(qld)+ones(n, 1)*qldUnc_a)+...
    iM*r*(qlddUnc_m*abs(qldd)+ones(n, 1)*qlddUnc_a);...
    %
    ];

bineq = bineq - unc_kinematics;

%-Robust against kinetics of the load
% unc_kinetics = 1/(eta*(1-etaUnc)*r)*...
%     [zeros(n,1);...
%     zeros(n, 1);...
%     tauSUnc_m*abs(torque)+ones(n,1)*tauSUnc_a;...
%     tauSUnc_m*abs(torque)+ones(n,1)*tauSUnc_a;...
%     tauSUnc_m*abs(torque)+ones(n,1)*tauSUnc_a;...
%     tauSUnc_m*abs(torque)+ones(n,1)*tauSUnc_a;...
%     tauSUnc_m*abs(torque)+ones(n,1)*tauSUnc_a;...
%     tauSUnc_m*abs(torque)+ones(n,1)*tauSUnc_a;];

unc_kinetics = 1/(eta*(1-etaUnc)*r)*...
    [tauSUnc_m*abs(torque)+ones(n,1)*tauSUnc_a;...
    tauSUnc_m*abs(torque)+ones(n,1)*tauSUnc_a;...
    tauSUnc_m*abs(torque)+ones(n,1)*tauSUnc_a;...
    tauSUnc_m*abs(torque)+ones(n,1)*tauSUnc_a;...
    tauSUnc_m*abs(torque)+ones(n,1)*tauSUnc_a;...
    tauSUnc_m*abs(torque)+ones(n,1)*tauSUnc_a;];

bineq = bineq - unc_kinetics;

%-Robust against uncertain implementation // Uncertainty in compliance.
%AineqRobust = [Aineq_mElong; Aineq_mTroque; Aineq_torVelo];
AineqRobust = [Aineq_mTroque; Aineq_torVelo];
unc_compli  = abs(AineqRobust*ones(n, 1)*comUncNonLinear);
% unc_compli_A  = max((AineqRobust)*ones(n, 1)*comUnc, 0);
% unc_compli_B  = max(-(AineqRobust)*ones(n, 1)*comUnc, 0);
% unc_compli    = max(unc_compli_A, unc_compli_B)

bineq = bineq - unc_compli;

end

