%% Script to generate periodic trajectories using a Fourier decomposition 
clc, clearvars, close all

human.task = 'walking';         %-- Gait {running, walking, stair_ascent}
human.joint = 'knee';           %-- Joint to analyze {hip, knee, ankle}
human.mass = 1;                 %-- Mass of the user [Kg]

% fast_cadence, slow_cadence, normal_cadence    
human.wSpeed = 'fast_cadence';
localFun_GenerateAndSavePeriodicTrajectories(human)


%% -- Local function library
function localFun_GenerateAndSavePeriodicTrajectories(human)
mass = human.mass;                  %-- Mass of the user [Kg]
joint = human.joint;             %-- Joint to analyze {hip, knee, ankle}
nPoints = 1000;


if strcmp(human.task, 'running')
    load('../1_InputData_BiomechanicDatasets/Datasets/dataset_Novacheck.mat')    
    %--[2] T. F. Novacheck, The biomechanics of running,
    %--Gait Posture, vol. 7, no. 1, Jan. 1998.
    ciclePeriod = 1.3213 /2;
    
    %-- Generate nPoints points to generate a grid
    %-- Converting to [rad] for consistency with Winter data
    ql = novacheck_running.(joint).position *pi /180;
    %-- Multiplying by user weight and (-1) to get Dorsiflexion positive
    %   torque as a positive number [Nm].
    torque = -novacheck_running.(joint).torque*mass;
    time = linspace(0, ciclePeriod, nPoints).';
    textFile = sprintf('%s_%s_%dkg', ...
    human.task, human.joint, human.mass);
elseif strcmp(human.task, 'walking')
    load('../1_InputData_BiomechanicDatasets/Datasets/dataset_Winter.mat')
    wSpeed = human.wSpeed;
    %-------------------------------WINTER'S DATA SCALING
    if strcmp(wSpeed,'normal_cadence')
        %--Cadence, Steps per minute. Winter page 12
        spm = 105;
        %--Sample time
        sT = 60/spm*2/1001;
    elseif strcmp(wSpeed,'fast_cadence')
        %Cadence, Steps per minute. Winter page 12
        spm = 123;
        %--Sample time
        sT = 60/spm*2/1001;
    elseif strcmp(wSpeed,'slow_cadence')
        %Cadence, Steps per minute. Winter page 12
        spm = 87;
        %--Sample time
        sT = 60/spm*2/1001;
    else
        error('Please select a walking speed, e.g., normal, fast');
    end
    %-- Angle of the ankle joint,Positive Angle (deg) => Dorsiflexion.
    %-- Negative => Plantarflexion
    ql = level_walking.(wSpeed).(joint).position*pi/180;
    %-- Torque of the ankle joint, Negative Torque => Dorsiflexion.
    %-- Positive => Plantarflexion
    torque = -level_walking.(wSpeed).(joint).torque*mass;    
    time = [0:sT:sT*(1000)]';  
    %-Name of the file
    textFile = sprintf('%s_%s_%s_%dkg', ...
    human.task, human.wSpeed, human.joint, human.mass);    
elseif strcmp(human.task, 'stair_ascent')
    load('../1_InputData_BiomechanicDatasets/Datasets/dataset_Riener.mat')
    %-- Generate nPoints points to generate a grid
    %-- Converting to [rad] for consistency with Winter's data
    ql = riener.stairAscent.(joint).position *pi /180;
    %-- Multiplying by user weight and (-1) to get Dorsiflexion positive
    %   torque as a positive number [Nm].
    torque = -riener.stairAscent.(joint).torque*mass;
    time = riener.stairAscent.(joint).time;
    textFile = sprintf('%s_%s_%dkg', ...
    human.task, human.joint, human.mass);
else
    error('Select adequate taks')
end

periodic_ql = fFourierDecomposition(time, ql, 10);
periodic_torque = fFourierDecomposition(time, torque, 10);

%-- Find the lowest frequency different from 0
freqs = periodic_ql.freq;
delFreq = min( freqs (freqs > 0 ) );

time = linspace(0, 1/delFreq, nPoints);
ql = periodic_ql.signalAnalytic;
torque = periodic_torque.signalAnalytic;

syms x

qld = diff(ql,x);
qldd = diff(qld,x);
torqueD = diff(torque,x);
torqueDD = diff(torqueD,x);

qlDisc = zeros(1,nPoints);
qldDisc = zeros(1,nPoints);
qlddDisc = zeros(1,nPoints);
torqueDisc = zeros(1,nPoints);
torquedDisc = zeros(1,nPoints);
torqueddDisc = zeros(1,nPoints);

tic
for i = 1:length(time)
    qlDisc(i) = round(double(subs(ql,time(i))),10);
    qldDisc(i) = round(double(subs(qld,time(i))),10);
    qlddDisc(i) = round(double(subs(qldd,time(i))),10);
    torqueDisc(i) = round(double(subs(torque,time(i))),10);
    torquedDisc(i) = round(double(subs(torqueD,time(i))),10);
    torqueddDisc(i) = round(double(subs(torqueDD,time(i))),10);
end
toc

save(textFile, 'qlDisc', 'qldDisc', 'qlddDisc', 'torqueDisc',...
    'torquedDisc', 'torqueddDisc', 'time');
end