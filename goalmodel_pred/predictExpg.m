function [ out ] = predictExpg( modelFit , team1 , team2 )
% PREDICTEXPG Compute expected goals for team1 against team2 according to
% modelFit goal model
%
% Inputs :
%  modelFit  goal model object
%  team1  n by 1 cell of strings containing team1 names
%  team2  n by 1 cell of strings containing team2 names
%
% Output :
%  out structure
%    out.team1
%    out.team2
%    out.expg1  expected goals for team1
%    out.expg2  expected goals for team2
%
ee = lambdaPred(modelFit.parameters,team1,team2) ;
out = struct('team1',team1,'team2',team2,'expg1',num2cell(ee.expg1),'expg2',num2cell(ee.expg2)) ;

end

