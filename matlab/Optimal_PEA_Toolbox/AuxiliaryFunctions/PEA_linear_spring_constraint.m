function [Aineq, bineq] = PEA_linear_spring_constraint(trajectory)
%CONSTRAINTS_COMPLICANCE_TRAJECTORIES Summary of this function goes here
%   Detailed explanation goes here

% Length of the periodic motion
n = length(trajectory.torque) - 1;

% Create derivative matrix
D = spdiags([-ones(n,1), ones(n,1)], [-1, 0], n, n);
% diff_1 = x_1 - x_end (taking advantage of periodicity)
D(1, end) = -1;

%-Compliance is constant during the trajectory D\alpha == 0
Aineq_pos               = D;
[m_Aineq, ~]            = size(Aineq_pos);
bineq_pos               = zeros(m_Aineq, 1);
Aineq_neg               = -D;
[m_Aineq, ~]            = size(Aineq_neg);
bineq_neg               = zeros(m_Aineq, 1);
 
Aineq = [Aineq_pos; Aineq_neg];
bineq = [bineq_pos; bineq_neg];

end

