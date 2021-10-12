function [D] = fnumerical_derivative_matrix(refTime, derivative)
%% [D] = fnumerical_derivative_matrix(refTime)
%
% Function that calculates the matrix 'D' such that numerical
% differentiation can be represented as 'qd = D*q'.
%
% It handles constant and variable sample time. The vector 'refTime'
% inticates the vector of time for the vector 'q'.
%
% The function assumes that the vector 'q' representes periodic motion.
% This implies that q(refTime(1)) == q(refTime(end)) and
% period = refTime(end) - refTime(1).

period = refTime(end)-refTime(1);
%--Removing last element since it won't be required to calculate qd(end).
%--qd(end) = qd(1) from the periodic motion assumption.
refTime = refTime(1: end-1);
refTimeA = refTime;
refTimeB = [refTimeA(end)-period; refTimeA(1: end-1)];
%--Delta at time i is equal to ti-t(i-1) based on Stanford's hand out
deltTime = refTimeA - refTimeB;
%--Select the vector 'b' to accomodate for the proper order of the derivative
if derivative == 1
    b = [0;1;0];
elseif derivative == 2
    b = [0;0;1];
else
    error('Please insert a correct order for the derivative')
end
%--'a' is the vector of coefficients for the numerical differentiation.
a = zeros(length(deltTime), 3);
for i = 1: length(deltTime)
    if i == length(deltTime)
        A = [1 1 1;...
            -1 0 deltTime(1)/deltTime(i);...
            1/2 0 1/2*(deltTime(1)/deltTime(i))^2];
    else
        A = [1 1 1;...
            -1 0 deltTime(i+1)/deltTime(i);...
            1/2 0 1/2*(deltTime(i+1)/deltTime(i))^2];
    end
    a(i,:) = A\b;
end
%--Scaling the coefficients dividing by the sample time
for i = 1: length(a)
    if derivative == 1
        a(i,:) = a(i,:)./ (deltTime(i));
    elseif derivative == 2
        a(i,:) = a(i,:)./ (deltTime(i)^2);
    else
        error('Please insert a correct order for the derivative')
    end
end
%-- Exploiding sparsity pattern
D = spdiags(a(2:end, 1), -1, length(a), length(a)) + ...
    spdiags(a(:,2), 0, length(a), length(a)) + ...
    spdiags([0;a(:,3)], 1, length(a), length(a));
D(1,end) = a(1,1);
D(end,1) = a(end,3);
D = sparse(D);
end