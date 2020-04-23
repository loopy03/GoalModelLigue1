function out = lambdaPred(plist , team1 , team2)
% LAMBDAPRED Compute the expected value in goalmodel
% Compute the expected value for the poisson and negative binomial
% from a list of paramteres and data. This function is used many
% places, like when making predictions and in the negloklik function.
%
% Inputs:
%   plist  structure of parameters
%      plist.pTeam  table of attack and defence paramters
%      plist.hfa  home field advantage (optional)
%      plist.gamma  psychological underestimation factor 
%                  (Rue & Salvesen 1997) (optional)
%      plist.rho  Dixon-Coles adjustment paramters (Dixon & Coles 1996)
%                  (optional)
%   team1  n by 1 cell array
%   team2  n by 1 cell array
%   
% Output:
%   out  structure of expected goals
%      out.expg1  expected goals of team1
%      out.expg2  expected goals of team2
%

% Team 1 log-expected goals
eta1 = plist.intercept + plist.pTeam(team1,:).attack - plist.pTeam(team2,:).defence ;

if isfield(plist,'hfa')
    eta1 = eta1 + plist.hfa ;
end

% Team 2 log-expected goals
eta2 = plist.intercept + plist.pTeam(team2,:).attack - plist.pTeam(team1,:).defence ;

% Psychological underestimation factor (Rue & Salvesen 1997)
if isfield(plist,'gamma')
    deltas = delta_ab(plist.pTeam(team1,:).attack,plist.pTeam(team1,:).defence,plist.pTeam(team2,:).attack,plist.pTeam(team2,:).defence) ;
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

