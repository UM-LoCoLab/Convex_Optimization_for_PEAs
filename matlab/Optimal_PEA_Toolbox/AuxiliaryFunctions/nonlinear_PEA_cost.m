function [costParameters] = nonlinear_PEA_cost(robot, trajectory)

iM      = robot.iM;
r       = robot.r;
km      = robot.km;
bm      = robot.bm;
% bm = 0;
eta     = robot.eta;
time    = trajectory.time;
deltaT  = time(2) - time(1);
% tauUnc  = robot.tauUnc;

%-Constant matrix that computes numerical derivative as a matrix operation
D1 = fnumerical_derivative_matrix_FixedSampleRate(deltaT, ...
            length(time)-1, 1);

%- Loading trajectories
% ql      = trajectoryA.ql(1:end-1);
ql      = trajectory.ql(1:end-1);
qld     = trajectory.qld(1:end-1);
qldd    = trajectory.qldd(1:end-1);
%--Torque done by the spring on the load
torque      = -trajectory.torque(1:end-1);     % = tau_L = -tau_ImportedData
torqued     = -trajectory.torqued(1:end-1);
torquedd    = -trajectory.torquedd(1:end-1);

%-- From 03/14/19 Edgar's lab book
a = qld*r;  
B = -spdiags(zeros(size(r*torqued)), 0, length(torqued), length(torqued));
c = qldd*r;
D = sparse(zeros(size(diag(r*torquedd))));
e = iM*c + bm*a - torque/(eta*r);   % (?)
F = sparse((1/r)*diag(ql));
G = deltaT*(((F.')*F)/(km^2)+bm*((B.')*B));
h = deltaT*(2*e.'*F/(km^2)+2*bm*a.'*B);
w = deltaT*((e.')*e/(km^2)+bm*(a.')*(a)-torque.'*qld/eta);

%-Define the affine-quadratic representation as input of "quad_form" (cvx)
[m, n]  = size(F);
diagQ   = ones(n,1)*deltaT/(km^2);
Q_Fcvx  = spdiags(diagQ, 0, m, n);

[m, n]  = size(B);
diagQ   = ones(n,1)*deltaT*bm;
Q_Bcvx  = spdiags(diagQ, 0, m, n);

costParameters.a   = a;
costParameters.B   = B;
costParameters.c   = c;
costParameters.D   = D;
costParameters.e   = e;
costParameters.F   = F;
costParameters.G   = G;
costParameters.h   = h;
costParameters.w   = w;
costParameters.Q_Fcvx   = Q_Fcvx;
costParameters.Q_Bcvx   = Q_Bcvx;
end