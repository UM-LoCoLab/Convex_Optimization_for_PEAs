function [x] = sampleFromSphere(dimSp, radSp)
%SAMPLEFROMSPHERE returns a point in R^n that is within the n ball of
%a given radious.

n = ceil(rand(1,1)*dimSp);
r = sign(rand(1,1)-0.5)*radSp;
try
    x = sparse(n, 1, r, dimSp, 1);
catch e
    %Making sure the path is not modified from original state
    rmpath(genpath('Optimal_SEA_Toolbox'))
    %-Print the reason for the actual error
    e.getReport
end

% %   The definition of the sample space is in 
% %   https://en.wikipedia.org/wiki/N-sphere#cite_note-2
% if dimSp < 3
%     error('sampleFromSphere supports only dimensions higher or equal than 3')
% end
% 
% 
% theta = zeros(dimSp-1, 1);
% theta(1:dimSp-2) = rand(dimSp-2, 1)*pi;
% theta(dimSp-1) = rand(1, 1)*2*pi;
% % radSp = rand(1,1)*radSp;
% radSp = radSp;
% 
% x = zeros(dimSp, 1);
% for i = 1:dimSp
%     if i == 1
%         x(1) = radSp*cos(theta(i));
%     else
%         mul = 1;
%         for j = 1:i - 1
%             mul = mul*sin(theta(j));
%         end
%         if i == dimSp
%             x(i) = radSp*mul;  
%         else
%             x(i) = radSp*mul*cos(theta(i));  
%         end
%     end
% end



end

