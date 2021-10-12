function [robotRigid, rigidTrans] = evaluate_transmission(...
    robot, trajectory)
%This function increases the reduction ratio until the system is infeasible

%-Computing rigid case (nominal conditions)
robotRigid      = fCostRigidCase(robot, trajectory);
feasible        = robotRigid.feasible;
i = 1;
while(feasible) 
    rigidTrans.Energy(i)    = robotRigid.energy_dissipated;
    rigidTrans.r(i)         = robot.r;
    i                       = i+1;
    robot.r                 = robot.r+1;
    robotRigid              = fCostRigidCase(robot, trajectory);
    feasible                = robotRigid.feasible;    
end

robot.r         = robot.r-1;
robotRigid      = fCostRigidCase(robot, trajectory);
robotRigid.r    = robot.r;

end

