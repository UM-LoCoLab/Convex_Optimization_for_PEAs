%% Test drive of the thermal brushless ODE. 

%   The first state, x(1)  = Temp. Housing - Temp. Ambient
%   The second state, x(2) = Temp. Winding - Temp. Ambient
%   The input is the motor current squared (i^2)
%
%   For the definition of the dynamics check with Edgar. Here are the eqns.
%   dxdt = A*x + B*i^2 + Nl*i^2
%
%   Requires: thermalBrushlessODE, load_SEA_Parameters, loadDataFromWeb
%             download dataset folder for legacy purposed (check error from
%             (first run for details).
%
%   Oct. 2020 - Edgar Bolivar-Nieto
clc, clearvars

% Load experimental data
% w/ housing
data.imaging.housing.trial001.raw = dlmread('..\2_dataset_thermalModeling\1023_housing_3.txt');
data.voltage.housing.trial001.raw = load('..\2_dataset_thermalModeling\1023_housing_3.mat');

data.imaging.housing.trial001.time = data.imaging.housing.trial001.raw(:,1);
data.imaging.housing.trial001.time = data.imaging.housing.trial001.time/60; %convert to minutes
data.imaging.housing.trial001.time_raw =data.imaging.housing.trial001.time ;

data.imaging.housing.trial001.temp.winding_raw = (data.imaging.housing.trial001.raw(:,2)*2+data.imaging.housing.trial001.raw(:,3))/3; % adding wieghts
data.imaging.housing.trial001.temp.housing_raw = data.imaging.housing.trial001.raw(:,4);

data.time = data.imaging.housing.trial001.time(1:30:end);
data.housing = data.imaging.housing.trial001.temp.housing_raw(1:30:end);
data.winding = data.imaging.housing.trial001.temp.winding_raw(1:30:end);

% Defining motor parameters. Load parameters are irrelevant
requirements.motor      = 'ActPackDephy';
requirements.task       = 'running';
requirements.joint      = 'ankle';
requirements.userMass   = 70;         %69.1 Kg [Winter's dataset];
requirements.objFun     = 'Energy'; 

[motor, ~] = load_SEA_Parameters(requirements); % For legacy purposes

% Unwrap motor values
Rwa = motor.Rwa;
Rwh = motor.Rwh;     % [K/W]
Rha = motor.Rha;     % [K/W]
Cwa = motor.Cwa;     % [W*s/K]
Cha = motor.Cha;     % [W*s/K]
Ramb= motor.R;       % [Ohms] Resistance at room temperature
aCu = motor.alphaCu; % [1/100 percent per K]
im  = 8.0179;        % [Amps]
Ta  = data.housing(1);        % Ambient temperature in Celcius
Tmax = 100;
Tfinal = Tmax - Ta; % data.winding(end)-data.winding(1)
% 
% Rewrite values from Ung Hee's code
Rwh = 1.1842;
Rha = 3.5332;
Cwa = 36.0686;
Cha = 104.1564;
Twinding = Rwh*Cwa;
Tmotor   = Rha*Cha;

% Program the dynamics
A  = [-(Rha+Rwh)/(Rha*Rwh*Cha), 1/(Rwh*Cha);
     1/(Rwh*Cwa), -1/(Rwh*Cwa)];
 
% Check Edgar's notes for derivation. x(1) = Th-Ta, x(2) = Tw-Ta
dxdt = @(t,x) A*x + [0;im^2*Ramb*(1 + aCu*x(2))/Cwa];
% Initial condition and simulation time.
x0      = [0, data.winding(1) - data.housing(1)];
tspan   = [0 70*60];

%Simulate nonlinear!
[t, xnl]  = ode45(dxdt, tspan, x0); 

%Simulate linearized!
dxdt = @(t,x) A*x + [0;im^2*Ramb*(1 + aCu*Tfinal)/Cwa];
[tl, xl]  = ode45(dxdt, tspan, x0); 

% Interpolate XNL
xnlInt = interp1(t/60, xnl(:,2), data.time);
xlInt = interp1(tl/60, xl(:,2), data.time);


%Ploting
figure, hold on, title('Motor Winding Temperature')
plot(data.time, xnlInt+Ta, '-.')
plot(data.time, xlInt+Ta, '--')
plot(data.time, data.winding)
% plot(t/60, xnl(:,1)+Ta, '-.')
% plot(tl/60, xl(:,1)+Ta, '--')
% plot(data.time, data.housing)
% legend('W. Sim. N.Lin.', 'W. Sim. Lin.', 'W. Mea.','H. Sim. N.Lin.', 'H. Sim. Lin.', 'H. Mea.')
leg1 = legend('Simulated Nonlinear', 'Simulated Linear', 'Measured');
set(leg1,'Box','off')
xlabel('Time [minutes]')
ylabel(['Temperature [', char(176), 'C]'])


error = abs((xlInt+Ta) - (xnlInt+Ta));
percent = error./(xlInt+Ta)*100;

figure, plot(data.time, percent)




%Export files to pgf plots
%-- Export data for PGFPLOTS
%-- Step 1: Resample the trajectories. (Suggested 200 points per graph)
nPoints = 20;
timeInter = linspace(t(1), t(end), nPoints);
wNlin = interp1(t, xnl(:,2)+Ta, timeInter);
wlin  = interp1(tl, xl(:,2)+Ta, timeInter);
dataW = interp1(data.time, data.winding, timeInter./60);


%-- Step 2: Export .dat      
%-- Exporting elongation and torque
vec2Write = [(timeInter./60).', wNlin.', wlin.' dataW.'];
fid = fopen('FIG_temperature.dat','w');
fprintf(fid,'Time WindingNonlinear WindingLinear DataWinding\n');
fprintf(fid,'%f %f %f %f\n', vec2Write.');
fclose(fid);


%% Steady State Analysis

x2ss = data.winding(end) - data.winding(1);
x1ss = data.housing(end) - data.housing(1);
% From the dynamic equations cancelling the time derivatives
Ploss = (x2ss-x1ss)/Rwh + x2ss/Rwa;
Rss = Ramb*(1 + motor.alphaCu*x2ss);
imss  = sqrt(Ploss/Rss);

%% Output from Ung Hee's paper
% nohousing = 
% 
%   struct with fields:
% 
%                      R_WH: 0.5199
%                      R_WA: 4.2685e+03
%                      R_HA: 4.2433
%                      C_WA: 6.5473
%                      C_HA: 91.1028
%         winding_estimated: [69×1 double]
%              TH_estimated: [68×1 double]
%         time_steady_state: [1×100 double]
%                        rc: [1×1 struct]
%                  resistor: [1×1 struct]
%     time_constant_winding: 3.4038
%       time_constant_motor: 386.5748


% housing = 
% 
%   struct with fields:
% 
%                      R_WH: 1.1842
%                      R_WA: 9.1845e+03
%                      R_HA: 3.5332
%                      C_WA: 36.0686
%                      C_HA: 104.1564
%         winding_estimated: [71×1 double]
%              TH_estimated: [70×1 double]
%         time_steady_state: [1×100 double]
%                        rc: [1×1 struct]
%                  resistor: [1×1 struct]
%     time_constant_winding: 42.7112


