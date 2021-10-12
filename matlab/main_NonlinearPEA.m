%% Energy-Optimal Design of Nonlinear Parallel Elasticity
%  Sicong Guo, Edgar Bolívar - Oct. 2021.
close all
clear
clc

% Nuke workspace
clc, clearvars -except results* runMonteCarlo

%------------------------------- Define PEA parameters --------------------
%task(ref. trajectories): walking running CubicSpring CubicLinearSpring
requirements.task     = 'walking';
requirements.joint    = 'hip';
requirements.wSpeed   = 'fast';
requirements.userMass = 69.1;         %69.1 Kg [Winter's dataset]; 65kg
%-Motor: {'EC45', 'EC30', 'ILM85x26', 'ILM70x8', 'ILM85x04', 'ActPackDephy'}
requirements.motor = 'ActPackDephy';
%requirements.MCTrials   = 10;    %Trials for monte carlo simulation 7k

%---------------------------------- PERFORM EXPERIMENT---------------------
requirements.objFun = 'energy'; % energy, RMS Torque, RMS Velocity, RMS Elongation

try
    addpath(genpath('Optimal_PEA_Toolbox'))
    Leg2 = performExperiments(requirements, false);
    requirements.motor = 'ActPackDephy';
    OSLm = performExperiments(requirements, false);
    rmpath(genpath('Optimal_PEA_Toolbox'))
catch e
    %Making sure the path is not modified from original state
    rmpath(genpath('Optimal_PEA_Toolbox'))
    %-Print the reason for the actual error
    e.getReport
end
Le2StiNom = uniquetol(Leg2.linearPEA.stiffness);
OSLStiNom = uniquetol(OSLm.linearPEA.stiffness);
maxDuaLe2 = Leg2.linearPEA.dualMax;
indDuaLe2 = Leg2.linearPEA.dualInd;
LHSLe2 = Leg2.linearPEA.LHSActCons;
RHSLe2 = Leg2.linearPEA.RHSActCons;
maxDuaOSL = OSLm.linearPEA.dualMax;
indDuaOSL = OSLm.linearPEA.dualInd;
LHSOSL = OSLm.linearPEA.LHSActCons;
RHSOSL = OSLm.linearPEA.RHSActCons;

% Select values of stiffness that are unique within 1 Nm/rad
indUni = @(x) [true; abs(x(2:end) - x(1)) > 1];

Le2StiNom = Le2StiNom(indUni(Le2StiNom));
OSLStiNom = OSLStiNom(indUni(OSLStiNom));

fprintf('Objective - %s\n', requirements.objFun)
fprintf('Leg2 spring stiffness - Nominal: %4.2f [Nm/rad]\n', Le2StiNom)
fprintf('Leg2 cost nonlinear (nominal)  : %4.2f []\n', Leg2.nonlinearPEA.cost)
fprintf('Leg2 cost linear (nominal)     : %4.2f []\n', Leg2.linearPEA.cost)
fprintf('Leg2 cost rigid                : %4.2f []\n', Leg2.robotRigid.cost)
fprintf('Leg2 dual variable: %3.2d  indx: %d \n', maxDuaLe2, indDuaLe2)
fprintf('Leg2 - LHS: %3.3f - RHS: %3.3f\n', LHSLe2, RHSLe2)
fprintf('OSL spring stiffness - Nominal : %4.2f [Nm/rad]\n', OSLStiNom)
fprintf('OSL cost nonlinear (nominal)   : %4.2f []\n', OSLm.nonlinearPEA.cost)
fprintf('OSL cost linear (nominal)      : %4.2f []\n', OSLm.linearPEA.cost)
fprintf('OSL cost rigid                 : %4.2f []\n', OSLm.robotRigid.cost)
fprintf('OSL dual variable: %3.2d   indx: %d \n', maxDuaOSL, indDuaOSL)
fprintf('OSL - LHS: %3.3f - RHS: %3.3f\n', LHSOSL, RHSOSL)

figure()
plot(OSLm.linearPEA.elong,OSLm.linearPEA.tauElastic,'DisplayName','Linear PEA')
hold on
grid on
p=plot(Leg2.nonlinearPEA.elong,Leg2.nonlinearPEA.tauElastic,'Marker','d','DisplayName','Nonlinear PEA');
p.MarkerIndices = 1:100:length(OSLm.linearPEA.tauElastic);
l1 = legend('show');
legend 'boxoff'
set(l1,'Interpreter','latex','location','best');
xlabel('$\delta_{p}$ [rad]','Interpreter','latex')
ylabel('$\tau_{p}$ [N$\cdot$m]','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gcf, 'PaperUnits', 'inches');
x_width=3.5;y_width=2;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(gcf,'kp_profile','epsc')

figure()
plot(OSLm.trajectory.time,(1/OSLm.PEAParam.r)*OSLm.trajectory.torque,'LineWidth',2,'DisplayName','Load')
hold on
grid on
p=plot(Leg2.trajectory.time,Leg2.linearPEA.tauM,'Marker','o','DisplayName','Linear PEA');
p.MarkerIndices = 1:100:length(Leg2.linearPEA.tauM);
p=plot(Leg2.trajectory.time,Leg2.nonlinearPEA.tauM,'Marker','d','DisplayName','Nonlinear PEA');
p.MarkerIndices = 1:100:length(Leg2.nonlinearPEA.tauM);
l1 = legend('show');
legend 'boxoff'
set(l1,'Interpreter','latex','FontSize',14);
xlabel('$t$ (s)','Interpreter','latex','FontSize',16)
ylabel('$\tau_{m}$ (N$\cdot$m)','Interpreter','latex','FontSize',16)
title([requirements.wSpeed ' ' requirements.joint],'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','FontSize',16)

%% Rigid Actuator
%rng(0,'twister');
pos_data_rms = rms(Leg2.trajectory.ql);
vel_data_rms = rms(Leg2.trajectory.qld);
acc_data_rms = rms(Leg2.trajectory.qldd);
trq_data_rms = rms(Leg2.trajectory.torque);
dtrq_data_rms = rms(Leg2.trajectory.torqued);
n = length(Leg2.trajectory.ql);
error = 0.2;
mc_no = 10000;
Em_new = zeros(mc_no,1);
e_l = -1;
e_u = 1;
dt = Leg2.trajectory.time(2)-Leg2.trajectory.time(1);
for i=1:mc_no
    pos_data_new = Leg2.trajectory.ql + pos_data_rms*error*((e_u-e_l).*rand(n,1) + e_l);
    vel_data_new = Leg2.trajectory.qld + vel_data_rms*error*((e_u-e_l).*rand(n,1) + e_l);
    acc_data_new = Leg2.trajectory.qldd + acc_data_rms*error*((e_u-e_l).*rand(n,1) + e_l);
    trq_data_new = Leg2.trajectory.torque + trq_data_rms*error*((e_u-e_l).*rand(n,1) + e_l);
    
    e_new = Leg2.PEAParam.iM * Leg2.PEAParam.r * acc_data_new(1:end-1) + Leg2.PEAParam.bm * Leg2.PEAParam.r * vel_data_new(1:end-1) + (1/Leg2.PEAParam.r) * trq_data_new(1:end-1);
    F_new = sparse((1/Leg2.PEAParam.r)*diag(pos_data_new(1:end-1)));
    tau_m_new = e_new;
    dq_m_new = Leg2.PEAParam.r * vel_data_new(1:end-1);
    dq_m_new(end+1) = dq_m_new(1);
    tau_m_new(end+1) = tau_m_new(1);
    Em_new(i) = trapz(Leg2.trajectory.time, (Leg2.PEAParam.km^(-2)).*tau_m_new.^2 + tau_m_new.*dq_m_new);
end
tau_m_0 = Leg2.PEAParam.iM * Leg2.PEAParam.r * Leg2.trajectory.qldd + Leg2.PEAParam.bm * Leg2.PEAParam.r * Leg2.trajectory.qld + (1/Leg2.PEAParam.r) * Leg2.trajectory.torque;
dq_m_0 = Leg2.PEAParam.r * Leg2.trajectory.qld;
Em_0 = trapz(Leg2.trajectory.time, (Leg2.PEAParam.km^(-2)).*tau_m_0.^2 + tau_m_0.*dq_m_0);
Em_trial = mean(Em_new);
Em_std = std(Em_new);
Em_act_per = (Em_trial-Em_0)*100/Em_0;
Em_std_act = (Em_trial+2*Em_std-Em_0)*100/Em_0 - Em_act_per;
disp(['Rigid: Em_0: ' num2str(Em_0) '; Em_trial: ' num2str(Em_trial) ' (' num2str(Em_act_per) '%); stdev: ' num2str(Em_std_act) '%'])

%% Linear PEA
Em_new = zeros(mc_no,1);
for i=1:mc_no
    pos_data_new = Leg2.trajectory.ql + pos_data_rms*error*((e_u-e_l).*rand(n,1) + e_l);
    vel_data_new = Leg2.trajectory.qld + vel_data_rms*error*((e_u-e_l).*rand(n,1) + e_l);
    acc_data_new = Leg2.trajectory.qldd + acc_data_rms*error*((e_u-e_l).*rand(n,1) + e_l);
    trq_data_new = Leg2.trajectory.torque + trq_data_rms*error*((e_u-e_l).*rand(n,1) + e_l);
    
    e_new = Leg2.PEAParam.iM * Leg2.PEAParam.r * acc_data_new(1:end-1) + Leg2.PEAParam.bm * Leg2.PEAParam.r * vel_data_new(1:end-1) + (1/Leg2.PEAParam.r) * trq_data_new(1:end-1);
    F_new = sparse((1/Leg2.PEAParam.r)*diag(pos_data_new(1:end-1)));
    tau_m_new = e_new + F_new * Leg2.linearPEA.stiffness(1:end-1);
    dq_m_new = Leg2.PEAParam.r * vel_data_new(1:end-1);
    dq_m_new(end+1) = dq_m_new(1);
    tau_m_new(end+1) = tau_m_new(1);
    Em_new(i) = trapz(Leg2.trajectory.time, (Leg2.PEAParam.km^(-2)).*tau_m_new.^2 + tau_m_new.*dq_m_new);
end
Em_per = (Leg2.linearPEA.energy_total-Em_0)*100/Em_0;
disp(['L: Em_0: ' num2str(Em_0) '; Em_nom: ' num2str(Leg2.linearPEA.energy_total) ' (' num2str(Em_per) '%)'])
Em_trial = mean(Em_new);
Em_std = std(Em_new);
Em_act_per = (Em_trial-Em_0)*100/Em_0;
Em_std_act = (Em_trial+2*Em_std-Em_0)*100/Em_0 - Em_act_per;
disp(['Em_trial: ' num2str(Em_trial) ' (' num2str(Em_act_per) '%); stdev: ' num2str(Em_std_act) '%'])

%% Nonlinear PEA
Em_new = zeros(mc_no,1);
for i=1:mc_no
    pos_data_new = Leg2.trajectory.ql + pos_data_rms*error*((e_u-e_l).*rand(n,1) + e_l);
    vel_data_new = Leg2.trajectory.qld + vel_data_rms*error*((e_u-e_l).*rand(n,1) + e_l);
    acc_data_new = Leg2.trajectory.qldd + acc_data_rms*error*((e_u-e_l).*rand(n,1) + e_l);
    trq_data_new = Leg2.trajectory.torque + trq_data_rms*error*((e_u-e_l).*rand(n,1) + e_l);
    
    e_new = Leg2.PEAParam.iM * Leg2.PEAParam.r * acc_data_new(1:end-1) + Leg2.PEAParam.bm * Leg2.PEAParam.r * vel_data_new(1:end-1) + (1/Leg2.PEAParam.r) * trq_data_new(1:end-1);
    F_new = sparse((1/Leg2.PEAParam.r)*diag(pos_data_new(1:end-1)));
    tau_m_new = e_new + F_new * Leg2.nonlinearPEA.stiffness(1:end-1);
    dq_m_new = Leg2.PEAParam.r * vel_data_new(1:end-1);
    dq_m_new(end+1) = dq_m_new(1);
    tau_m_new(end+1) = tau_m_new(1);
    Em_new(i) = trapz(Leg2.trajectory.time, (Leg2.PEAParam.km^(-2)).*tau_m_new.^2 + tau_m_new.*dq_m_new);
end
Em_per = (Leg2.nonlinearPEA.energy_total-Em_0)*100/Em_0;
disp(['NL: Em_0: ' num2str(Em_0) '; Em_nom: ' num2str(Leg2.nonlinearPEA.energy_total) ' (' num2str(Em_per) '%)'])
Em_trial = mean(Em_new);
Em_std = std(Em_new);
Em_act_per = (Em_trial-Em_0)*100/Em_0;
Em_std_act = (Em_trial+2*Em_std-Em_0)*100/Em_0 - Em_act_per;
disp(['Em_trial: ' num2str(Em_trial) ' (' num2str(Em_act_per) '%); stdev: ' num2str(Em_std_act) '%'])

figure()
plot(Leg2.trajectory.time,Leg2.trajectory.torque,'LineWidth',2,'DisplayName','ideal')
hold on
plot(Leg2.trajectory.time,trq_data_new,'DisplayName','true')
l1 = legend('show');
legend 'boxoff'
set(l1,'Interpreter','latex');
xlabel('$t$ (s)','Interpreter','latex','FontSize',16)
ylabel('$\tau_{l}$ (N$\cdot$m)','Interpreter','latex','FontSize',16)
title([requirements.wSpeed ' ' requirements.joint],'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','FontSize',16)

figure()
plot(Leg2.trajectory.time,Leg2.trajectory.ql,'LineWidth',2,'DisplayName','ideal')
hold on
grid on
plot(Leg2.trajectory.time,pos_data_new,'DisplayName','noised')
l1 = legend('show');
legend 'boxoff'
set(l1,'Interpreter','latex','location','best');
xlabel('$t$ (s)','Interpreter','latex')
ylabel('$q_{l}$ (m)','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gcf, 'PaperUnits', 'inches');
x_width=3.5;y_width=2;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(gcf,'ql-noised','epsc')

figure()
plot(OSLm.trajectory.time,(1/OSLm.PEAParam.r)*OSLm.trajectory.torque,'LineWidth',2,'DisplayName','Load')
hold on
grid on
p=plot(Leg2.trajectory.time,Leg2.linearPEA.tauM,'Marker','o','LineWidth',2,'DisplayName','Linear PEA (ideal)');
p.MarkerIndices = 1:100:length(Leg2.linearPEA.tauM);
p=plot(Leg2.trajectory.time,Leg2.nonlinearPEA.tauM,'Marker','d','LineWidth',2,'DisplayName','Nonlinear PEA (ideal)');
p.MarkerIndices = 1:100:length(Leg2.nonlinearPEA.tauM);
p=plot(Leg2.trajectory.time,tau_m_new,'DisplayName','Nonlinear PEA (noised)');
p.MarkerIndices = 1:100:length(tau_m_new);
l1 = legend('show');
legend 'boxoff'
set(l1,'Interpreter','latex','location','northeast');
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$\tau_{m}$ [N$\cdot$m]','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gcf, 'PaperUnits', 'inches');
x_width=3.5;y_width=2;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(gcf,'taum-profile','epsc')

%% Monte Carlo Simulation Convergence Test (if needed)
%mc_convergence_test;