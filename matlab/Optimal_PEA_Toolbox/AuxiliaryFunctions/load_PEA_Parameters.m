function [robot, trajectory] = load_PEA_Parameters(requirements)
%LINKPARAMETERS Summary of this function goes here
%   Detailed explanation goes here

loadDataFromWeb()

motor = requirements.motor;

%Efficiency of the transmission
%Remember the asymmetry in the efficiency
robot.eta = 1;
%Max elongation [rad]
robot.maxElong = 90*pi/180; % Elongation in degrees

% Load parameters based on motor and transmission requirements
if strcmp(motor,'EC45')
    robot.kt = 36.9/1000;             %Torque constant motor [Nm/A]
    robot.R = 0.608;            %Terminal resistance phase to phase [ohms]
    robot.iM = 181/1000/100^2;  %Inertia of the motor [Kg*m^2]
    robot.r = 283;              %Reduction ratio of the transmission
    robot.bm = robot.iM*2;      %Viscous damping of the motor [Nm*s]
    robot.km = robot.kt/sqrt(robot.R);      %Motor constant
    %-
    robot.voltage = 36;
    %--Constraints
    robot.maxVelo = robot.voltage/robot.kt;     
    robot.RatedTorque = 0.128;                  
    robot.peakTorque = robot.RatedTorque*5;   %Peak torque [Nm]
elseif strcmp(motor,'EC30')
    % 30V in leg 1, motor is rated for 24V
    robot.voltage = 30;                     %Voltage power supply [Volts]
    robot.kt = 13.6/1000;                   %Torque constant motor [Nm/A]
    robot.R = 0.102;             %Terminal resistance phase to phase [ohms]
    robot.iM = 33.3/1000/100^2;             %Inertia of the motor [Kg*m^2]
    %(%knee 360, ankle 720, 425)
    robot.r = 600;                          %Reduction ratio of the transmission
    robot.bm = robot.iM/2;                  %Viscous damping of the motor [Nm*s]
    robot.km = robot.kt/sqrt(robot.R);      %Motor constant
    %--Constraints
    robot.maxVelo = robot.voltage/robot.kt; %Max motor velocity 1500 RPM [rad/sec]
    robot.RatedTorque = 0.135;              %Peak torque [Nm]
    robot.peakTorque = robot.RatedTorque*2.5;   %Peak torque [Nm]
    robot.L = 16.3e-6;                          %Inductance [uH]
elseif strcmp(motor,'ILM85x26')
    robot.kt = 0.24;                              %Torque constant motor [Nm/A]
    robot.R = 323/1000;                           %Terminal resistance [ohms]
    robot.r = 24;                           %reduction ratio
%     rotorInerRD = 1.15*(1e-2)^2;            %Inertia of the rotor (RoboDrive) [Kg*m^2]
%     rotorInerToby = 13098.7/1000*(1e-3)^2;  %Inertia of the rotor (Toby) [Kg*m^2]
%     robot.iM = rotorInerRD+rotorInerToby;   %Inertia of the motor [Kg*m^2]
%     robot.bm = robot.iM/2;                  %Viscous friction mechanism

    % Using: "Design and Validation of a Powered Knee-Ankle Prosthesis with
    % High-Torque, Low-Impedance Actuators" T-RO. Elery et al.
    robot.iM = 0.0696/(robot.r^2);          %Inertia of the motor [Kg*m^2]    
    robot.bm = 0.4169/(robot.r^2);          %Viscous friction [N.m.s/rad]
    robot.km = robot.kt/sqrt(robot.R);      %Motor constant
    %--Constraints
    robot.maxVelo = 1500*2*pi/60;       %Max motor velocity 1500 RPM [rad/sec]
    robot.peakTorque = 8.3;             %Peak torque [Nm]
    robot.RatedTorque = 2.6;            %Peak torque [Nm]
    robot.L = 920e-6;                   %Inductance [uH]
    robot.voltage = 48;                 %Voltage power supply [Volts]
%     robot.voltage = 24;                 %Voltage power supply [Volts]
elseif strcmp(motor,'ILM70x18')
    robot.kt = 0.18;                              %Torque constant motor [Nm/A]
    robot.R = 655/1000;                           %Terminal resistance [ohms]
    rotorInerRD = 0.34*(1e-2)^2;            %Inertia of the rotor (RoboDrive) [Kg*m^2]
    rotorInerToby = 13098.7/1000*(1e-3)^2;  %Inertia of the rotor (Toby) [Kg*m^2]
    robot.iM = rotorInerRD+rotorInerToby;   %Inertia of the motor [Kg*m^2]
    robot.r = 38;
    robot.bm = robot.iM/2;                  %Drag torque
    robot.km = robot.kt/sqrt(robot.R);                  %Motor constant
    %--Constraints
    robot.maxVelo = 2100*2*pi/60;           %Max motor velocity 2100 RPM [rad/sec]
    robot.peakTorque = 4;                   %Peak torque [Nm]
    robot.RatedTorque = 1.25;               %Peak torque [Nm]
    robot.L = 1350e-6;                       %Inductance [uH]
    robot.voltage = 48;                     %Voltage power supply [Volts]
elseif strcmp(motor,'ILM85x04')
    robot.kt    = 0.04;                        %Torque constant motor [Nm/A]
    robot.R     = 138/1000;                    %Terminal resistance [ohms]
    rotorInerRD     = 0.28*(1e-2)^2;           %Inertia of the rotor (RoboDrive) [Kg*m^2]
    rotorInerToby   = 13098.7/1000*(1e-3)^2;   %Inertia of the rotor (Toby) [Kg*m^2]
    robot.iM    = rotorInerRD+rotorInerToby;   %Inertia of the motor [Kg*m^2]
    robot.r     = 160;
    robot.bm    = robot.iM/2;                  %Drag torque
    robot.km    = robot.kt/sqrt(robot.R);      %Motor constant
    %--Constraints
    robot.maxVelo       = 9000*2*pi/60;        %Max motor velocity 2100 RPM [rad/sec]
    robot.peakTorque    = 1.2;                 %Peak torque [Nm]
    robot.RatedTorque   = 0.43;                %Peak torque [Nm]
    robot.L             = 120e-6;              %Inductance [uH]
    robot.voltage       = 48;                  %Voltage power supply [Volts]
elseif strcmp(motor,'ActPackDephy')
    % Data from [1] U. H. Lee, C. Pan, and E. J. Rouse, “Empirical Characterization of 
    % a High-performance Exterior-rotor Type Brushless DC Motor and Drive,” 
    % in 2019 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS), 2019, pp. 8018–8025.
    robot.kt    = 0.14;                        %Torque constant motor [Nm/A]
%     robot.R     = 186/1000;                    %Terminal resistance [ohms]
    robot.R     = 279/1000;                    %Terminal resistance [ohms]
    robot.iM    = 0.12/1000;                   %Inertia of the motor [Kg*m^2]
    robot.r     = 50;                          %It is non constant
    robot.bm    = 0.16/1000;                   %Drag torque [Nm / (rad/s)]
    robot.km    = robot.kt/sqrt(robot.R);      %Motor constant
    robot.Rwa   = 3416.3;
    robot.Rwh   = 1.1;                         %Winding-Housing [K/W]
    robot.Rha   = 3.5;                         %Housing-Ambient [K/W]
    T_w         = 65;
    T_m         = 670.6;
    robot.Cwa   = T_w /robot.Rwa;              %Winding-Housing [W*s/K]
    robot.Cha   = T_m /robot.Rha;              %Housing-Ambient [W*s/K]
    %--Constraints
%     robot.voltage       = 36;                  %Voltage power supply [Volts]
    robot.voltage       = 48;                  %Voltage power supply [Volts]
    robot.maxVelo       = robot.voltage/robot.kt;   %Max motor velocity [rad/sec]
    robot.peakTorque    = 28.7*robot.kt;       %Peak torque [Nm]
    robot.RatedTorque   = 7.7*robot.kt;        %Max. continuous current [Nm]
    robot.L             = 138e-6;              %Inductance [uH]    
end
robot.noLoadSpeed = robot.voltage/robot.kt;
robot.stallTorque = robot.voltage*robot.kt/robot.R;
robot.alphaCu     = 0.00393;                  %[1/100 percent per K]

% Load the task requirements
task    = requirements.task;
joint   = requirements.joint;
mass    = requirements.userMass;

% Relative path to the reference trajectories
relPath = '1_InputData_BiomechanicDatasets/ReferenceTrajectories/';

if ( strcmp(task, 'walking') || strcmp(task, 'running')...
        || strcmp(task, 'stair_ascent') )    
    % Uncertainty values for a robust solution
    trajectory.mBar = mass;             %Nominal mass of user [kg]
    trajectory.mUnc = 8.8;              %Uncertainty in mass [+- kg]
    %-- Select walking speed from user or set to default
    if (isfield(requirements, 'wSpeed') && strcmp(task, 'walking'))
        wSpeed = requirements.wSpeed;
        trajectory.wSpeed   = wSpeed;
        nameFile = sprintf('%s_%s_%s_1kg', task, wSpeed, joint);
        pathFile = fullfile(relPath, nameFile);
        eval( append('load ', pathFile) );
    else
        nameFile = sprintf('%s_%s_1kg', task, joint);
        pathFile = fullfile(relPath, nameFile);
        eval( append('load ', pathFile) );
    end
    trajectory.time = time.';
    trajectory.ql = qlDisc.';
    trajectory.qld = qldDisc.';
    trajectory.qldd = qlddDisc.';
    trajectory.torque = torqueDisc.'*mass;
    trajectory.torqued = torquedDisc.'*mass;
    trajectory.torquedd = torqueddDisc.'*mass;
elseif ( strcmp(task, 'CubicSpring') )
    nameFile = 'CubicSpring5k';
    pathFile = fullfile(relPath, nameFile);
    eval( append('load ', pathFile) );
    trajectory.time = time.';
    trajectory.ql = qlDisc.';
    trajectory.qld = qldDisc.';
    trajectory.qldd = qlddDisc.';
    trajectory.torque = torqueDisc.';
    trajectory.torqued = torquedDisc.';
    trajectory.torquedd = torqueddDisc.';
elseif ( strcmp(task, 'CubicLinearSpring') )
    nameFile = 'CubicSpring5';
    pathFile = fullfile(relPath, nameFile);
    eval( append('load ', pathFile) );
    trajectory.time = time.';
    trajectory.ql = qlDisc.';
    trajectory.qld = qldDisc.';
    trajectory.qldd = qlddDisc.';
    trajectory.torque = torqueDisc.';
    trajectory.torqued = torquedDisc.';
    trajectory.torquedd = torqueddDisc.';
else
    error('Provide a proper reference task')
end

trajectory.joint    = joint;
trajectory.mass     = mass;
trajectory.task     = task;
trajectory.objFun   = requirements.objFun;

end


