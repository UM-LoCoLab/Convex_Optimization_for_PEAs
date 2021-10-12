function [optiRobot] = optimal_PEA_Spring(motorTransmission, ...
    trajectory, spring)
% optimal_PEA_Spring  Returns optimal nonlinear PEA using nominal and
% uncertain values of the trajectory and robot parameters.
%
%   OptimalSEA = optimal_PEA_Spring(motorTransmission, trajectory, spring)
%
%   -- INPUTS:
%           motorTransmission: [struct]
%               Parameters of the electric motor and robot system.
%           trajectory: [struct] Desired load trajectory.
%
%   -- OUTPUTS:
%           optiRobot   [struct] Optimal PEA parameters and energy
%                                measurements

%--Energy optimization
[cost] = nonlinear_PEA_cost(motorTransmission, trajectory);

% DEFINE ACTUATOR CONSTRAINTS AND CONSTRUCTION CONSTRAINTS

%--Constraint to satisfy actuator constraints
[Aineq_const, bineq_const, L] = PEA_actuator_constraints(...
    motorTransmission, trajectory, cost);

%-Define the kind of spring
if (strcmp(spring, 'nonlinear'))
    [AineqSp, bineqSp] = PEA_nonlinear_spring_constraint(trajectory);      % In process
elseif (strcmp(spring, 'linear'))
    [AineqSp, bineqSp] = PEA_linear_spring_constraint(trajectory);         %
else
    error('Please enter an adequate value for the spring.')
end

%- Programming plan:
%------------------------------------- INTERVENTION OF THE CODE STARTS HERE
% clear Aineq bineq
% Aineq = AineqSp;
% bineq = bineqSp;
%--- Impose RMS constraint
% Ball centered at e_nom with radious equivalnet to 10% change
n = length(cost.e);
radSqr = ( abs(cost.e)*1 ).^(2);
P = spdiags(radSqr, 0, n,n);
Pinv = spdiags(radSqr.^(-1), 0, n,n);
Psqr = spdiags(radSqr.^(1/2), 0, n,n);

if not(strcmp(trajectory.objFun, 'energy') || ...
        strcmp(trajectory.objFun, 'RMS Velocity') || ...
        strcmp(trajectory.objFun, 'RMS Torque') ||...
        strcmp(trajectory.objFun, 'RMS Elongation'))
    error('Select adequate string for the cost function.')
end

%     %----------------------------------------------------------------------
%     % Solving using Gurobi directly
%     model.Q     = cost.G;
%     model.obj   = cost.h;
%     model.A     = Aineq;
%     model.rhs   = bineq;
%     model.sense = '<';
%     gurobi_write(model, 'qp.lp');
%     params.outputflag = 1;
%     results     = gurobi(model, params);
%     cvx_status  = 'Solved';
%     kp      = results.x;
%     %----------------------------------------------------------------------
%     F = cost.F;
%     e = cost.e;
%     tauRms = motorTransmission.RatedTorque;
%     n = length(cost.h);
%     
%     cvx_begin SDP
%     cvx_solver mosek    %-Mosek, Gurobi, SDPT3, sedumi
% %     variable kp(n, 1)
%     variables kp(n, 1) lam gam
%     minimize(quad_form( kp, cost.G ) + cost.h*kp + cost.w)
%     subject to
%     Aineq*kp <= bineq;
% %     quad_form(kp, F.'*F) + 2*e.'*F*kp + e.'*e <= tauRms^2*n;
%     quad_form(kp, F.'*F) + (-gam) <= tauRms^2*n;    
%     lam  >= 0;
%     [(lam*Pinv-speye(n,n)), (-F*kp-lam*Pinv*e);...
%       (-F*kp-lam*Pinv*e).', (lam*(e.'*Pinv*e-1)-gam)] >= 0;
%     cvx_end
%     costCVX = quad_form( kp, cost.G ) + cost.h*kp + cost.w;
%     
%     %-- Checking if the solution of cvx is useful
%     if not(strcmp(cvx_status, 'Solved'))
%         warning(['CVX status:', cvx_status])
%     end   
    
n = length(cost.h);
F = cost.F;
e = cost.e;
tauRms = motorTransmission.RatedTorque;
oneOverSqrtN = 1/sqrt(n);


%---------------------------------------------FORMULATION CVX PROGRAM
cvx_begin
    cvx_solver Mosek            %-Mosek, Gurobi, SDPT3, sedumi
    cvx_precision default       %-best, default
    dual variable rmsTorCon;
    dual variable actCon;
    dual variable sprCon;

    variables kp(n, 1) lam gam
    expression costCVX        
    % Select cost function according to input
    if(strcmp(trajectory.objFun, 'RMS Velocity'))
        expression qmdOpti
        qmdOpti  = cost.a + cost.B*kp;
        costCVX = norm(qmdOpti)*oneOverSqrtN;
    elseif(strcmp(trajectory.objFun, 'RMS Torque'))
        expression tauMOpti
        tauMOpti = cost.e + cost.F*kp;
        costCVX = norm(tauMOpti)*oneOverSqrtN;
    elseif(strcmp(trajectory.objFun, 'energy'))
        costCVX = quad_form( kp, cost.G ) + cost.h*kp + cost.w;
    elseif  (strcmp(trajectory.objFun, 'RMS Elongation'))
        costCVX = norm(L*kp)*oneOverSqrtN;
    end
    minimize(costCVX)
    subject to
        sprCon : AineqSp*kp <= bineqSp;
        actCon : Aineq_const*kp <= bineq_const;
%         rmsTorCon : quad_form(kp, F.'*F) + 2*e.'*F*kp + e.'*e...
%             <= tauRms^2*n;
%         quad_form(kp, F.'*F) + (-gam) <= tauRms^2*n;
%         lam  >= 0;
%         [(lam*Pinv-speye(n,n)), (-F*kp-lam*Pinv*e);...
%             (-F*kp-lam*Pinv*e).', (lam*(e.'*Pinv*e-1)-gam)] >= 0;
cvx_end
%--------------------------------------------------------------------------

%-- Checking if the solution of cvx can be used
if not(strcmp(cvx_status, 'Solved'))
    warning(['CVX status:', cvx_status])
end
% Re-evaluating the cost (just in case)
if(strcmp(trajectory.objFun, 'RMS Velocity'))
    costCVX = norm(qmdOpti)*oneOverSqrtN;
elseif(strcmp(trajectory.objFun, 'RMS Torque'))
    costCVX = norm(tauMOpti)*oneOverSqrtN;
elseif(strcmp(trajectory.objFun, 'energy'))
    costCVX = quad_form( kp, cost.G ) + cost.h*kp + cost.w;
elseif  (strcmp(trajectory.objFun, 'RMS Elongation'))
    costCVX = norm(L*kp)*oneOverSqrtN;
end

%- Loading parameters
r       = motorTransmission.r;
time    = trajectory.time;
% deltaT  = time(2) - time(1);

%- Loading trajectories
ql      = trajectory.ql(1:end-1);
qld     = trajectory.qld(1:end-1);
% qld      = trajectory.qld(1:end-1);
%--Torque done by the spring on the load
torque      = -trajectory.torque(1:end-1);
torqued     = -trajectory.torqued(1:end-1);

%-Index when elongation is the closest to zero
minTorqueIndx = find(abs(-ql) == min(abs(-ql)));

%-Defining spring torque and setting torque equal to zero when elongation is
% close to zero.
tauP = cumtrapz(time(1:end-1), -qld.*kp);
tauP = tauP - tauP(minTorqueIndx);
trajectory.tau_P = tauP;
trajectory.tau_P(end+1) = tauP(1);
trajectory.e_vec = cost.e;
trajectory.F_matr = cost.F;
trajectory.k_p = kp;

%qm   = (ql - tauP)*r;
qm   = r*ql;
qmd  = cost.a+cost.B*kp;
qmdd = cost.c+cost.D*kp;

qm(end + 1)     = qm(1);
qmd(end + 1)    = qmd(1);
qmdd(end + 1)   = qmdd(1);

%-- Post processing the results
optiRobot = fPostProcessing(...
    qm, qmd, qmdd, cvx_status, motorTransmission, trajectory);
optiRobot.stiffness = kp;
%Periodic
optiRobot.stiffness(end+1) = optiRobot.stiffness(1);

%-- Checking constraints
[feasible, violatedConstraint, SpeedTorquePlot] = fIsMotorFeasible(...
    motorTransmission, optiRobot, sprintf('%s PEA', spring));
optiRobot.feasible = feasible;
optiRobot.violatedConstraint = violatedConstraint;
optiRobot.plot = SpeedTorquePlot;

% Testing robust RMS
dualVariables = [actCon; rmsTorCon];
[maxDua, indDua] = max(dualVariables);
if (indDua == length(dualVariables))
    warning('The RMS torque constraint is active')
else
    LHS = Aineq_const(indDua,:)*kp;
    RHS = bineq_const(indDua);
end

optiRobot.cosMat = cost;
optiRobot.Psqr = Psqr;
optiRobot.P = P;
optiRobot.Pinv = Pinv;
optiRobot.dualMax = maxDua;
optiRobot.dualInd = indDua;
optiRobot.LHSActCons = LHS;
optiRobot.RHSActCons = RHS;

%- Troubleshooting cost
optiRobot.cost = costCVX;

end