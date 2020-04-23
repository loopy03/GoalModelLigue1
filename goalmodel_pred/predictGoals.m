function [ res ] = predictGoals( modelFit , team1 , team2 , varargin )
% PREDICTGOALS Summary of this function goes here
% 
% Inputs :
%  modelFit  goal model object
%  team1  n by 1 cell of strings containing team1 names
%  team2  n by 1 cell of strings containing team2 names
%  Optional (Name,Value) :
%    'MaxGoal'  maximum goal difference accept for probabilites computation
%
%  Output :
%    res  length(team1) by 1 cell of probabilies matrices
%
% defaults for optional parameters
lwrx = 25;

% parse optional input parameters
v = 1;
while v < numel(varargin)
  switch lower(varargin{v})
  case 'maxgoal'
    assert(v+1<=numel(varargin));
    v = v+1;
    lwrx = varargin{v};
  otherwise
    error('Unsupported parameter: %s',varargin{v});
  end
  v = v+1;
end

% Predict the expected goals
out = predictExpg(modelFit,team1,team2) ;
expg = [out.expg1; out.expg2]' ;

% find the upper limit of where to evaluate the probability function
upperProb = 0.999 ;

% Compute maxGoal
if ismember(lower(modelFit.model),{'poisson','gaussian','cmp'})
    maxGoal = poissinv(upperProb,modelFit.maxGoal) ;
elseif strcmpi(modelFit.model,'negbin')
    r = 1/modelFit.parameters.dispersion ;
    p = r/(r+modelFit.maxGoal) ;  
    maxGoal = nbininv(upperProb,r,p) ;
end
maxGoal = max(lwrx,maxGoal) ;

res = cell(length(team1),1) ;

% Compute the probabilities matrices containing the probability of each
% goal difference
for ii=1:length(team1)
    if ismember(lower(modelFit.model),{'poisson','gaussian'})
        resTmp = poisspdf(0:maxGoal,expg(ii,1))' * poisspdf(0:maxGoal,expg(ii,2)) ;
    elseif strcmpi(modelFit.model,'negbin')
        p1 = r/(r+expg(ii,1)) ;  
        p2 = r/(r+expg(ii,2)) ;
        resTmp = nbinpdf(0:maxGoal,r,p1)' * nbinpdf(0:maxGoal,r,p2) ;
    elseif strcmpi(modelFit.model,'cmp')
        % TO DO
        ll1 = lambdaCMP(expg(ii,1),modelFit.parameters.dispersion,'fast') ;
        ll2 = lambdaCMP(expg(ii,2),modelFit.parameters.dispersion,'fast') ;
        
        resTmp = cmppdf(0:maxGoal,ll1,modelFit.parameters.dispersion)' .* cmppdf(0:maxGoal,ll2,modelFit.parameters.dispersion) ;
    end
    
    % Dixon-Coles adjustement
    if isfield(modelFit.parameters,'rho')
        scalingVec = tau( [0 1 0 1] , [0 0 1 1] , expg(ii,1)*ones(4,1) , expg(ii,2)*ones(4,1) , modelFit.parameters.rho ) ;
        scalingMat = [scalingVec(1) scalingVec(3); scalingVec(2) scalingVec(4)] ;
        resTmp(1:2,1:2) = resTmp(1:2,1:2) .* scalingMat ;
    end
        
    % Normalize
    resTmp = resTmp ./ sum(sum(resTmp)) ;
    
    res{ii} = resTmp ;
end
end

