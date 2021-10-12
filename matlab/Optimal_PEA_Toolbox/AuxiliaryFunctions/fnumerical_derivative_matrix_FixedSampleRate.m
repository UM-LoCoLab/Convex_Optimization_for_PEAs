function [D] = fnumerical_derivative_matrix_FixedSampleRate(delTime, n, ...
    derivative)
%% Function to calculate the derivative as a matrix multiplication
%-- It is assumed that the trajectory is periodic, i.e., ql(1) == ql(end)

delT = delTime;

%-- using 3 coefficients. Third order approximation
if derivative == 1
    e = ones(n,1);
    D = spdiags([-e e], [-1,1], n, n);
    %-- Exploiting periodicity
    D(end, 1) = 1;
    D(1, end) = -1;
    D = D./(2*delT);    
elseif derivative == 2
    e = ones(n,1);
    D = spdiags([e -2*e e], -1:1, n, n);
    %-- Exploiting periodicity
    D(end, 1) = 1;
    D(1, end) = 1;
    D = D./(delT^2);
else
    error('Please insert a correct order for the derivative')
end


% %-- using 5 coefficients. Fifth order approximation
% 
% 
% %-- Taylor table. Please refer to 08/06/18 of Edgar's lab notebook for
% %derivation. Basic idea: Tailor expansion about f(i) of f(i-2), f(i-1), f(i+1),
% %f(i+2). Expressing the result as [a0, a1, a2, a3, a4]*[f(i-2); f(i-1); f(i);
% % f(i+1); f(i+2)]
% 
% A = [1, 1, 1, 1, 1;...
%     -2 -1 0 1 2;...
%     2 sym(1/2) 0 sym(1/2) 2;...
%     -sym(8/6) -sym(1/6) 0 sym(1/6) sym(8/6);...
%     sym(16/24) sym(1/24) 0 sym(1/24) sym(16/24)];
% 
% if derivative == 1
%     %-- b selects the desired order of derivative. 
%     % [(f(i); f(i)'; f(i)''; f(i)'''; f(i)'''']
%     b = [0; 1; 0; 0; 0];
%     coefA = A\b;
%     coefA = double(coefA);
%     e = ones(n,1);
%     D = spdiags([e*coefA(1) e*coefA(2) e*coefA(3) e*coefA(4) e*coefA(5)],...
%         -2:2, n, n);
%     %-- Exploiting periodicity    
%     D(1, end) = coefA(2);
%     D(1, end-1) = coefA(1);
%     D(2, end) = coefA(1);
%     D(end-1, 1) = coefA(5);
%     D(end, 1) = coefA(4);
%     D(end, 2) = coefA(5);
%     D = D./(delT);    
% elseif derivative == 2
%     %-- b selects the desired order of derivative. 
%     % [(f(i); f(i)'; f(i)''; f(i)'''; f(i)'''']
%     b = [0; 0; 1; 0; 0];
%     coefA = A\b;
%     coefA = double(coefA);
%     e = ones(n,1);
%     D = spdiags([e*coefA(1) e*coefA(2) e*coefA(3) e*coefA(4) e*coefA(5)],...
%         -2:2, n, n);
%     %-- Exploiting periodicity    
%     D(1, end) = coefA(2);
%     D(1, end-1) = coefA(1);
%     D(2, end) = coefA(1);
%     D(end-1, 1) = coefA(5);
%     D(end, 1) = coefA(4);
%     D(end, 2) = coefA(5);
%     D = D./(delT^2);
% else
%     error('Please insert a correct order for the derivative')
% end


end