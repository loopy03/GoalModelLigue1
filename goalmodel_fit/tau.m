function [ taufactor ] = tau( goals1 , goals2 , lambda1 , lambda2 , rho )
% Compute the Dixon-Coles adjustment to the log-likelihood.

taufactor  = ones(length(goals1),1) ;
 
if rho ~= 0
    y0 = goals2 == 0 ;
    y1 = goals2 == 1 ;
    x0 = goals1 == 0 ;
    x1 = goals1 == 1 ;
    
    idx00 = x0 & y0 ;
    taufactor(idx00) = 1 - lambda1(idx00) .* lambda2(idx00) * rho ;
    idx01 = x0 & y1 ;
    taufactor(idx01) = 1 + lambda1(idx01) * rho ;
    idx10 = x1 & y0 ;
    taufactor(idx10) = 1 + lambda2(idx10) * rho ;
    idx11 = x1 & y1 ;
    taufactor(idx11) = 1 - rho ;
end
end

