function elo_mat = create_elongation_matrix(trajectory)
%CONSTRAINTS_COMPLICANCE_TRAJECTORIES Summary of this function goes here
%   Detailed explanation goes here
%   Mod: This function actually generates tau_p vector from delta_p = q_l
%   vector

time    = trajectory.time;
deltaT  = time(2) - time(1);

%--Torque done by the spring on the load
%torque      = -trajectory.torque(1:end-1);
%torqued     = -trajectory.torqued(1:end-1);
delta_p      = -trajectory.ql(1:end-1);
delta_pd     = -trajectory.qld(1:end-1);

%-Index when elongation is the closest to zero
minElongIndx = find(abs(delta_p) == min(abs(delta_p)));
delta_p       = delta_p - delta_p(minElongIndx);

%-Ordering matrix with respect to order of elongation
n   = length(delta_p);
[~, indx] = sort(delta_p);
i = (1:length(indx));
ordMatr = sparse(i, indx, ones(length(indx), 1), n, n);

%-Creating torque matrix using Trapezoidal integration.
elo_mat         = zeros(n, n);
elo_mat(2,1)    = delta_pd(1);
elo_mat(2,2)    = delta_pd(2);
for i = 3:length(delta_pd)
    elo_mat(i,:) = elo_mat(i-1,:);
    elo_mat(i,i) = delta_pd(i);
end

%Initial value always divided by 2
elo_mat(:,1)    = elo_mat(:,1)./2;
%Diagonal value always divided by 2
for i = 1:length(delta_pd)
    elo_mat(i,i) = elo_mat(i,i)/2;
end

% %-Substracting torque at index "minElongIndx"
% del_vector  = elo_mat(minElongIndx, :);
% for i = 1:length(delta_pd)
%     elo_mat(i,:) = elo_mat(i, :) - del_vector;
% end

elo_mat     = elo_mat*deltaT;

%-Sorting elongation with respect to torque
elo_mat     = ordMatr*elo_mat;

end

