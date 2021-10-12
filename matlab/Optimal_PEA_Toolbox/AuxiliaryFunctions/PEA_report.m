function PEA_report(customTrajectory, motor, robot, ...
    names, varargin)

try
    %open the file for writing
    file_ID = fopen('PEA_report.txt', 'wt');
    
    %motor data
    fprintf(file_ID, 'Trajectory Details:');
    fprintf(file_ID, '\n\t-Task: \t%s', customTrajectory.task);
    fprintf(file_ID, '\n\t-Joint: %s', customTrajectory.joint);
    if not (strcmpi(customTrajectory.task, 'CubicSpring'))
        fprintf(file_ID, '\n\t-Subject Mass: \t%d kg', customTrajectory.userMass);
    end
    fprintf(file_ID, '\n\nMotor Details:');
    fprintf(file_ID, '\n\t-Motor Name: \t\t\t\t%s', motor);
    fprintf(file_ID, '\n\t-Motor Torque Constant: \t%.4f Nm/A', robot.kt);
    fprintf(file_ID, '\n\t-Motor Terminal Resistance: %.4f Ohms', robot.R);
    fprintf(file_ID, '\n\t-Motor Inertia: \t\t\t%.2f mg*m^2', robot.iM*1000);    
    fprintf(file_ID, '\n\t-Viscous Friction: \t\t\t%.2f mNm*s', robot.bm*1000);
    fprintf(file_ID, '\n\t-Motor Constant: \t\t\t%.4f Nm/W^(1/2)', robot.km);

    fprintf(file_ID, '\n\nTransmission Details:');
    fprintf(file_ID, '\n\t-Reduction Ratio: \t\t\t%d', robot.r);
    fprintf(file_ID, '\n\t-Efficiency: \t\t\t\t%d%%', robot.eta*100);
    
    
    fprintf(file_ID, '\n\n================================================================================');
    fprintf(file_ID, '\nMotor Energy Consumption and Peak Power');
    fprintf(file_ID, '\n================================================================================');
    
    fprintf(file_ID, '\n\n');
    
    %SEA vs. rigid actuator chart
    for i = 1:length(names)
        if i == 1
            fprintf(file_ID, '\n\t\t\t\t\t\t\t%s\t', names{i});
        else
            fprintf(file_ID, '|\t %s \t', names{i});
        end        
    end
    fprintf(file_ID, '\n--------------------------------------------------------------------------------');
    
    for i = 1:length(varargin)
        if i == 1
            fprintf(file_ID, '\nTotal Energy Consumption: \t%.2f J\t\t', ...
                varargin{i}.energy_total);
        else
            fprintf(file_ID, '|\t%.2f J\t\t', varargin{i}.energy_total);
        end
    end
    
    for i = 1:length(varargin)
        if i == 1
            fprintf(file_ID, '\nEnergy Dissipated: \t\t\t%.2f J\t\t', ...
                varargin{i}.energy_dissipated);
        else
            fprintf(file_ID, '|\t%.2f J\t\t', varargin{i}.energy_dissipated);
        end
    end
    
    for i = 1:length(varargin)
        if i == 1
            fprintf(file_ID, '\nPeak Power: \t\t\t\t%.2f W\t', ...
                varargin{i}.peakPower);
        else
            fprintf(file_ID, '|\t%.2f W\t', varargin{i}.peakPower);
        end
    end
    
%         fprintf(file_ID, '\nEnergy Dissipated: \t\t\t%.2f J\t\t\t\t\t|\t%.2f J\t\t|\t%.2f J', ...
%         optimalSEAs.energy_dissipated, optimalLinear.energy_dissipated, robotRigid.energy_dissipated);
%     fprintf(file_ID, '\nPeak Power: \t\t\t\t%.2f W\t\t\t\t|\t%.2f W\t|\t%.2f W', ...
%         optimalSEAs.peakPower, optimalLinear.peakPower, robotRigid.peakPower);
    
    %close the file
    fclose(file_ID);
catch
    warning('Error: I/O Malfunction');
end




end

