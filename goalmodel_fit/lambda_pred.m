function [ out ] = lambda_pred(plist , team1 , team2 , x1 , x2)
% Compute the expected value in goalmodel.
%
% Compute the expected value for the poisson and negative binomial
% from a list of paramteres and data. This function is used many
% places, like when making predictions and in the negloklik function.


% Team 1 log-expected goals
eta1 = plist.intercept + plist.attack(team1) - plist.defense(team2) ;

if isfield(plist,'hfa')
    eta1 = eta1 + plist.hfa ;
end

% Team 2 log-expected goals
eta2 = plist.intercept + plist.attack(team2) - plist.defense(team1) ;

% Psychological underestimation factor (Rue & Salvesen 1997)
if isfield(plist,'gamma')
    deltas = delta_ab(plist.attack(team1),plist.defense(team1),plist.attack(team2), plist.defense(team2)) ;
    gamma_delta = plist.gamma*deltas ; 
    eta1 = eta1 - gamma_delta ;
    eta2 = eta2 + gamma_delta ;
end

% Additional covariates.
%%% TO DO

% Link function
lambda1 = exp(eta1) ;
lambda2 = exp(eta2) ;

out = struct('expg1',lambda1,'expg2',lambda2) ;
end

