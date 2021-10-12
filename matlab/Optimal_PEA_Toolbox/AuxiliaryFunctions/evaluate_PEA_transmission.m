function [PEA, PEATrans] = evaluate_PEA_transmission(...
    robot, trajectory, springKind)
%This function increases the reduction ratio until the system is infeasible

PEA             = optimal_SEA_Spring(robot, trajectory, springKind);
feasible        = PEA.feasible;
i = 1;
while(feasible) 
    PEATrans.Energy(i)  = PEA.energy_dissipated;
    PEATrans.r(i)       = robot.r;
    i                   = i+1;
    robot.r             = robot.r+1;
    PEA                 = optimal_SEA_Spring(robot, trajectory,springKind);
    feasible            = PEA.feasible;    
end

robot.r         = robot.r-1;
PEA             = optimal_SEA_Spring(robot, trajectory,springKind);
PEA.r           = robot.r;

end

