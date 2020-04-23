function [ prob ] = predictResult( modelFit , team1 , team2)
% PREDICTRESULT Compute the expected probabilities team1 vs. team2
% according to modelFit goal model
%   
% Predict goals
dgoals = predictGoals(modelFit,team1,team2) ;

prob = zeros(length(team1),3) ;
for ii=1:length(team1)
    prob(ii,1) = sum(sum(tril(dgoals{ii},-1))) ;
    prob(ii,2) = trace(dgoals{ii}) ;
    prob(ii,3) = sum(sum(triu(dgoals{ii},1))) ;
end
end

