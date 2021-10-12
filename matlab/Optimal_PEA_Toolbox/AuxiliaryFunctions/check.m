clc
clear
close all

%% Sanity check
% create tau_p
trajectory.time = 0:0.01:2;
trajectory.ql = sin(2*pi*trajectory.time);
trajectory.qld = 2*pi*cos(2*pi*trajectory.time);
L_trq = create_elongation_matrix(trajectory);
trq = L_trq * (2*ones(length(trajectory.time)-1,1));
[~,argidx]=sort(-trajectory.ql);

figure()
plot(trajectory.time,-2*trajectory.ql(argidx))
hold on
plot(trajectory.time(1:end-1),trq)