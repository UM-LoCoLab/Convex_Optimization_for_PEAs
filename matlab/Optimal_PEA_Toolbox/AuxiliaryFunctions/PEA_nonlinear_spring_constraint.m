function [Aineq, bineq] = PEA_nonlinear_spring_constraint(trajectory)
%CONSTRAINTS_COMPLICANCE_TRAJECTORIES Summary of this function goes here
%   Detailed explanation goes here

trq_mat = create_elongation_matrix(trajectory);

% Creating integral constraint int{taus*alpha*deltaTaus}
% torqueElastic = -trajectory.torque;
% deltaTauS = [0;torqueElastic(2:end)-torqueElastic(1:end-1)];
% intVector = [torqueElastic(2:end-1).*deltaTauS(2:end-1)+...
%     torqueElastic(2:end-1).*deltaTauS(3:end)]./2;
% intVector = [torqueElastic(1)*deltaTauS(2)/2; intVector;...
%     torqueElastic(end)*deltaTauS(end)/2];

%torqueElastic       = -trajectory.torque(1:end-1);
elongpElastic       = trajectory.ql(1:end-1);
deltaP           = diff(elongpElastic);
integralVector      = [deltaP; 0] + [0; deltaP];
integralVector(1)   = integralVector(1)/2;
integralVector(end) = integralVector(end)/2;
integralVector      = -integralVector.*elongpElastic;
% integralVector = [];

%-Elongation is monotonically increasing
Aineq               = -(trq_mat(2:end, :) - trq_mat(1:end-1, :));
%-Adding energy equal to zero as two complementary inequalities
Aineq               = sparse([Aineq; integralVector.'; -integralVector.']);
[m_Aineq, ~]        = size(Aineq);
bineq               = zeros(m_Aineq, 1);

end

