% Script to create a torque elastic trajectory that justifies the use of non linear elastic elements
clc, clear variables, close all

%----------- Creating reference motion
[ql, torque, time, robot] = localfunGenerateMotion;

% --------------Finding the maximum and minimum points of torque out of the EOM
syms x;
qld = diff(ql,x);
qldd = diff(qld,x);
qld3 = diff(qldd,x);
qld4 = diff(qld3,x);
torqued = diff(torque,x);
torquedd = diff(torqued,x); 

% period = time(end)-time(1);
% maxTorque = double(subs(ql,0))^3*robot.k3;
% timeMaxTorque = vpasolve(torque == maxTorque,x,[time(1) time(end)]);
% timeZeroTorque = vpasolve(torque == 0,x,[time(1) time(end)]);
% 
% %--------------Defining the sample time and creating corresponding time vector
% samTime = 1e-3;
% %time = [0:samTime:double(timeMaxTorque)];
% time = double(linspace(timeZeroTorque,timeZeroTorque+period,1000));

nSamples2 = length(time);
qlDisc = zeros(1,nSamples2);
qldDisc = zeros(1,nSamples2);
qlddDisc = zeros(1,nSamples2);
torqueDisc = zeros(1,nSamples2);
torquedDisc = zeros(1,nSamples2);
torqueddDisc = zeros(1,nSamples2);

tic
for i = 1:length(time)
    x = time(i);
    torqueDisc(i) = round(double(subs(torque)),10);
    torquedDisc(i) = round(double(subs(torqued)),10);
    torqueddDisc(i) = round(double(subs(torquedd)),10);
    qlDisc(i) = round(double(subs(ql)),10);
    qldDisc(i) = round(double(subs(qld)),10);
    qlddDisc(i) = round(double(subs(qldd)),10);
end
toc

save('CubicLinearSpring.mat','qlDisc','qldDisc','qlddDisc',...
    'torqueDisc', 'torquedDisc', 'torqueddDisc', 'time');

%% Local function library
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ql,torque,time, robot] = localfunGenerateMotion

k3 = 40;
k1 = 10;

robot.k3 = k3;

T = 1e-4;                   %Sample time
L = 5000;                   %Length of the signal for fourier
timeTest = (0:L-1)*T;       %Time vector
iniElong = pi/4;

[t,x] = ode45(@fLoadDynamicsOnlyMass,timeTest,[iniElong,0]);

ql = x(:,1);
qld = x(:,2);
% --------------------------------------- Arranging reference signal
pointCrossing = qld(2:end).*qld(1:end-1);
indxPointCrossing = find(pointCrossing<0);
indxEnd = indxPointCrossing(2);
time = t(1:indxEnd);
ql = ql(1:indxEnd);

[reconstructedSignal] = fFourierDecomposition(time,ql,3);
clear ql
ql = reconstructedSignal.signalAnalytic;
time = linspace(0,1/reconstructedSignal.freq(end),1000);
clear torque
torque = -ql^3*k3 - ql*k1;            %Remember that fplot does not accurately represent the function
% torque = -ql*k3;                    %Remember that fplot does not accurately represent the function
% figure, fplot(ql,torque,[0,time(end)/2],'MeshDensity',2), hold on, grid on
close all
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%