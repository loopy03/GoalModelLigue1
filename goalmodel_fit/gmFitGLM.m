function [ out ] = gmFitGLM( goals1 , goals2 , team1 , team2 , model , varargin )
% GMFITGLM Fit models using the built-in glm.fit() function.
% This is an internal function that in some cases provides the final parameter
% estimates. When model = 'negbin' or 'cmp', or when there are fixed
% coefficients, this function provides only starting values for later
% optimization. Therefore the output of this function (coefficeints, loglik)
% does not neccecarily reflect the input.
%
% Inputs:
%   goals1  n by 1 double of goals scored by team1
%   goals2  n by 1 double of goals scored by team2
%   team1  n by 1 cell array
%   team2  n by 1 cell array
%   model  string containing the goal model resolution method ('poisson',
%          'negbin', 'cmp' or 'gaussian')
%   Optional (Name,Value) :
%        'hfa' logical (default : true) indicated whether a home field 
%        advantage is considered for team1
%        'weights' empty (default) or n by 1 double that determine the 
%        influence of each match on the final parameter estimates.  
% 
% Output:
%   out structure of fitglm results 
%        out.pList  structure of parameters (pTeam, intercept, hfa, sigma)
%        out.loglikelihood  loglikelihood of the model distribution
%        out.AIC  Akaike information criterion
%        out.nParEst  number of estimated parameters
%
% Defaults for optional parameters
hfa = true ;
weights = [] ;

% Parse optional input parameters
v = 1;
while v < numel(varargin)
  switch lower(varargin{v})
  case 'hfa'
    assert(v+1<=numel(varargin));
    v = v+1;
    hfa = varargin{v};
  case 'weights'
    assert(v+1<=numel(varargin));
    v = v+1;
    weights = varargin{v};
  otherwise
    error('Unsupported optional input parameter: %s',varargin{v});
  end
  v = v+1;
end

% Prepare some
allTeams = unique([team1;team2]) ;
nTeams = length(allTeams) ;
team1Stacked = [team1;team2] ;
team2Stacked = [team2;team1] ;

% Response vector
yy = [goals1;goals2] ;

% Attack and defence matrices
xmatA = zeros(2*length(goals2),nTeams-1) ;
xmatD = xmatA ;

for ii = 2:nTeams
    t1Idx = strcmp(allTeams(ii),team1Stacked) ;
    t2Idx = strcmp(allTeams(ii),team2Stacked) ;
    xmatA(:,ii-1) = t1Idx ;
    xmatD(:,ii-1) = -t2Idx ;
end

% Sum-to-zero constraint for the first team
xmatA(strcmp(allTeams(1),team1Stacked),:) = -1 ;
xmatD(strcmp(allTeams(1),team2Stacked),:) = -1 ;
colNames = allTeams(2:end) ;

% Intecept 
intercept = ones(2*length(goals2),1) ;

% Add home field advantage
if hfa
    xmatHfa = [ones(length(goals2),1);zeros(length(goals2),1)] ;
    xmat = [intercept xmatHfa xmatA xmatD] ;
    varNames = ['intercept';'hfa';strcat('attack_',colNames);strcat('defence_',colNames);'goals'] ;
else
    xmat = [intercept xmatA xmatD] ;
    varNames = ['intercept';strcat('attack_',colNames);strcat('defence_',colNames);'goals'] ;
end

% Add additional covariates
% TO DO !

% Preparation to fitglm
if ismember(lower(model),{'poisson','negbin','cmp'})
    glmDistribution = 'poisson' ;
elseif strcmpi(model,'gaussian')
    glmDistribution = 'normal' ;
end

varNames = matlab.lang.makeValidName(varNames);
lastwarn('Nothing to see here', 'this:is:not:a:warning');
if isempty(weights)
    glmRes = fitglm(xmat,yy,'Distribution',glmDistribution,'Link','log','Intercept',false,'VarNames',varNames);
else
    glmRes = fitglm(xmat,yy,'Distribution',glmDistribution,'Link','log','Intercept',false,'VarNames',varNames','Weights',[weights(:);weights(:)]);
end

% Check for warning
[lastMsg, lastID] = lastwarn;
warningWasThrown = ~strcmp(lastID, 'this:is:not:a:warning');
if warningWasThrown
    converged = false ;
    warningMsg = lastMsg ;
else
    converged = true ;
    warningMsg = '' ;
end

% A few parsing
idxA = cellfun('isempty',strfind(glmRes.Coefficients.Properties.RowNames,'attack_')) ;
attackParams = glmRes.Coefficients.Estimate(~idxA) ;
attackParams = [-sum(attackParams);attackParams] ;

idxD = cellfun('isempty',strfind(glmRes.Coefficients.Properties.RowNames,'defence_')) ;
defenceParams = glmRes.Coefficients.Estimate(~idxD) ;
defenceParams = [sum(defenceParams);defenceParams] ;

pTeam = table(attackParams,defenceParams,'VariableNames',{'attack','defence'},'RowNames',allTeams) ;
pList = struct('pTeam',pTeam,'intercept',glmRes.Coefficients.Estimate('intercept')) ;

if strcmpi(model,'gaussian')
    pList.sigma = std(glmRes.Residuals.Raw) ;
end

if hfa
    pList.hfa = glmRes.Coefficients.Estimate('hfa') ;
end

loglik = glmRes.LogLikelihood ;
AIC = glmRes.ModelCriterion.AIC ;

out = struct('pList',pList,'loglikelihood',loglik,'AIC',AIC,...
        'nParEst',size(xmat,2),'converged',converged,'warningMsg',warningMsg) ;
end

