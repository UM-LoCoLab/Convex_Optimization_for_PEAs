function [dxdt] = fLoadDynamicsOnlyMass(t,x)
%FLOADDYNAMICS Summary of this function goes here
%   Detailed explanation goes here

iL = 0.125;                 %Kg*m^2

u = torqueElastic(t,x);

dxdt(1,1) = x(2);
dxdt(2,1) = -1/iL*u;

end

function u = torqueElastic(t,x)

k3 = 40;
k1 = 10;

% u = k3*x(1)^3;
u = k3*x(1)^3 + k1*x(1);
end