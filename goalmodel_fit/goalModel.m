function [ out ] = goalModel( goals1 , goals2 , team1 , team2 , varargin )
% GOALMODEL Fitting models for team1 (scoring goals1) vs. team2 
%    (scoring goals2) 
%    This function is used to fit models of goals scored in sports 
%    competitions. At a minimum this function estimates 'attack' and 
%    'defence' ratings for all teams, but other covariates can be 
%    included, as well as other adjustments. The underlying statistical 
%    model can be either a Poisson, Negative Binomial, or Gaussian model.
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
%          advantage is considered for team1
%        'dc'  logical (default : false). If TRUE an adjustment for low 
%          scoring goals is included in the model.
%        'rs'  logical (default : false). If TRUE an adjustment for teams 
%          to over and under estimate the opponent.
%        'weights' empty (default) or n by 1 double that determine the 
%          influence of each match on the final parameter estimates.
%        'model' String indicating whether the goals follow a 'poisson' 
%          model (default), a Negative Binomial ('negbin') or a Gaussian 
%          ('gaussian) model.
% 
% Output:
%   out structure of results 
%        out.parameters  structure of parameters (pTeam, intercept, hfa, etc.)
%        out.loglikelihood  loglikelihood of the model distribution
%        out.AIC  Akaike information criterion
%        out.nParEst  number of estimated parameters
%        out.rSquared  rSquared estimation
%
warningMsg = '' ;

% Defaults for optional parameters
hfa = true ; dc = false ; rs = false ;
model = 'poisson' ;
weights = [] ;

% Parse optional input parameters
v = 1;
while v < numel(varargin)
  switch lower(varargin{v})
  case 'hfa'
    assert(v+1<=numel(varargin));
    v = v+1;
    hfa = varargin{v};
  case 'dc'
    assert(v+1<=numel(varargin));
    v = v+1;
    dc = varargin{v};
  case 'rs'
    assert(v+1<=numel(varargin));
    v = v+1;
    rs = varargin{v};
  case 'weights'
    assert(v+1<=numel(varargin));
    v = v+1;
    weights = varargin{v};
  case 'model'
    assert(v+1<=numel(varargin));
    v = v+1;
    model = varargin{v};  
  otherwise
    error('Unsupported optional input parameter: %s',varargin{v});
  end
  v = v+1;
end

% Stop if length conditions are not respect
if length(goals1)~=length(goals2)
    error('Error in goalModel.m : goals1 and goals2 must be the same length')
end
if length(goals2)~=length(goals1)
    error('Error in goalModel.m : goals and team must be the same length')
end
if length(team1)~=length(team2)
    error('Error in goalModel.m : team1 and team2 must be the same length')
end
if length(goals1)<1
    error('Error in goalModel.m : goals and team must not be empty')
end
if isnumeric(goals1)==0
    error('Error in goalModel.m : goals1 must be numerical values')
end
if isnumeric(goals2)==0
    error('Error in goalModel.m : goals2 must be numerical values')
end
if ismember(lower(model),{'poisson','negbin','gaussian','cmp'})==0
    error('Error in goalModel.m : unknown model')
end
  
% Check if data are suitable for DC adjustement
if dc
    if isempty(find(goals1<=1 & goals2<=1,1))
        error('Error in goalModel.m : Dixon-Coles adjustment is not applicable when there are no instances both teams scoring 1 goal or less.')
    end
    if strcmpi(model,'gaussian')
        error('Error in goalModel.m : Dixon-Coles adjustment does not work with a Gaussian model.')
    end    
end

% TO DO : check for connected data

% Check if the weights vector is suitable
if isempty(weights)==0
    if all(isnan(weights))~=0
        error('Error in goalModel.m : weights contains NaN values')
    end
    if isnumeric(weights)==0
        error('Error in goalModel.m : weights must be numerical values')
    end
    if length(goals1)~=length(weights)
        error('Error in goalModel.m : weights, team and goals must be the same length')
    end
    if all(weights>0)==0
        error('Error in goalModel.m : weights values must be positive')
    end
    if all(weights==0)
        error('Error in goalModel.m : weights values are all equal to 0')
    end
end

% TO DO : chack for additional covariates

% Male sure that team are at cell type
team1 = cell(team1) ;
team2 = cell(team2) ;

% Some useful quantities.
allTeams = unique([team1;team2]) ;
nTeams = length(allTeams) ;
nGames = length(goals1) ;

% If it is sufficient to fit the model with glm.fit().
mDefault = ismember(lower(model),{'poisson','gaussian'}) & dc== false & rs==false ;

% Which fitter method wad used
if mDefault
    fitter = 'fitglm' ;
else
    fitter = 'gm' ;
end

% Fit a model with fitglm. Some model classes are compatible with this function.
% If not, the results from fitting this simples model is used as starting values.

% TODO: If all attack and defence and hfa and intercept are fixed, this is not
% needed for fitting the rest of the models. This step should be skiped
% to save computing time.

% TODO: Maybe even negbin models can be fitted with this approach, with
% the glm.nb() function from the MASS pacakge.

% compute the fitglm function
gmFitGlmRes = gmFitGLM(goals1,goals2,team1,team2,model,'hfa',hfa,'weights',weights) ;

if mDefault
    if gmFitGlmRes.converged==0
        wTmp = 'figlm.m did not converge. Parameters estimates are unreliable.' ;
        warningMsg = [warningMsg,' / ',wTmp] ;
        warning(wTmp)
    end
    
    % TO DO : boudary info
    
    pList = gmFitGlmRes.pList ;
    logLikelihood = gmFitGlmRes.loglikelihood ;
    nParEst = gmFitGlmRes.nParEst ;
    AIC = gmFitGlmRes.AIC ;
    converged = gmFitGlmRes.converged ;
else
    % Inital values from gmFitGlmRes
    pListInit = gmFitGlmRes.pList ;
    pListInit.pTeam = table2array(pListInit.pTeam(2:end,:)) ;
    
    if dc
        pListInit.rho = 0.01 ;
    end
    if rs
        pListInit.gamma = 0.0 ;
    end
    if strcmpi(model,'negbin')
        % on log scale suring estimation
        % start must be close to zero (almost Poisson)
        pListInit.dispersion = -10 ;
    elseif strcmpi(model,'gaussian')
        % on log scale suring estimation
        pListInit.sigma = log(pListInit.sigma) ;
    elseif strcmpi(model,'cmp')
        % on log scale during estimation.
        % start with a value exp(0)=1, which is the same as Poisson.
        pListInit.dispersion = 0 ;
    end
    
    % Transformation structure -> double array
    parInitsCell = struct2cell(pListInit) ;
    a=cell2mat(parInitsCell(1)) ;
    parInitsCell{1} = a(:) ; 
    parInits = cell2mat(parInitsCell)' ;

    % Optimization
    optimFun = @(params) (negLogLik(params,goals1,goals2,team1,team2,hfa,model,pListInit,weights)) ;

     if strcmpi(model,'cmp')
        % Utilisation de fmincon pour pouvoir ajouter des boundaries car
        % sinon, calcul de fminunc plante...
        options = optimoptions(@fmincon,'Algorithm','interior-point','MaxIterations',1000,'MaxFunctionEvaluations',1000,'StepTolerance',10e-3,...
                             'CheckGradients',true);
        lb = -Inf*ones(1,length(parInits)) ;
        ub = Inf*ones(1,length(parInits)) ;
        lb(end) = log(0.7) ; % Fixe la dispersion à être au dessus de 0.7 
        ub(end) = log(10) ; % Fixe la dispersion à être en dessous de 10 
        [optimRes,optimValue,exitflag,info] = fmincon(optimFun,parInits,[],[],[],[],lb,ub,[],options) ;
    
    else
        options = optimoptions(@fminunc,'Algorithm','quasi-newton','StepTolerance',10e-3);
        [optimRes,optimValue,exitflag,info] = fminunc(optimFun,parInits,options) ;
    end
    
    converged = true ;
    if exitflag == 0
        converged = false ;
        wmTmp = 'fminunc.m did not converged. Parameters estimates may be unreliable.' ;
        warningMsg = [warningMsg,' / ',wmTmp] ;
        warning(wmTmp)
        warning(info.message)
    end
    
    % Transformation double array -> structure
    allTeams = unique([team1;team2]) ;
    parCell{1} = optimRes(1:numel(pListInit.pTeam))' ;
    parCell = [parCell;num2cell(optimRes(numel(pListInit.pTeam)+1:end)')] ;
    a=cell2mat(parCell(1)) ;
    parCell{1} = [a(1:end/2) a(end/2+1:end)] ;
    pList = cell2struct(parCell,fields(pListInit)) ;
    
    % Calculate missing attack and deffence parameters
    pList.pTeam = [-sum(pList.pTeam); pList.pTeam] ;

    pList.pTeam = array2table(pList.pTeam,'RowNames',allTeams,'VariableNames',{'attack','defence'}) ;
         
    logLikelihood = -optimValue ;
    nParEst = length(optimRes) ;
    AIC = 2*nParEst - 2*logLikelihood ;
    
    % Rescale dispersion
    if isfield(pList,'dispersion')
        pList.dispersion = exp(pList.dispersion) ;
    end
    
    % Rescale sigma
    if isfield(pList,'sigma')
        pList.sigma = exp(pList.sigma) ;
    end
end

% Compute R squared
allGoals = [goals1;goals2] ;
if isempty(weights)
    meanGoals = mean(allGoals) ;
else
    meanGoals = sum([weights(:);weights(:)].*allGoals)/(2*sum(weights)) ;
end

% Deviances neede for R squared
if strcmpi(model,'poisson')
    if isempty(weights)
        logLikelihoodSatured = sum(log(poisspdf(allGoals,allGoals))) ;
        logLikelihoodNull = sum(log(poisspdf(allGoals,meanGoals))) ;
    else
        logLikelihoodSatured = sum(log(poisspdf(allGoals,allGoals)).*[weights(:);weights(:)]) ;
        logLikelihoodNull = sum(log(poisspdf(allGoals,meanGoals)).*[weights(:);weights(:)]) ;
    end
elseif strcmpi(model,'negbin')
    % TO DO : must be checked carrefully
    if isempty(weights) 
%         dispersion0Tmp = mle(allGoals,'pdf',@(allGoals,r)nbinpdf(allGoals,dispersion0Tmp,dispersion0Tmp/(dispersion0Tmp+meanGoals)),'start',1) ;
%         logLikelihoodSatured = sum(log(nbinpdf(allGoals,10^6,10^6/(allGoals+10^6))),'start',1) ;
%         logLikelihoodNull = sum(log(nbinpdf(allGoals,dispersion0Tmp,pdispersion0Tmp/(dispersion0Tmp+meanGoals)))) ;
          logLikelihoodSatured = NaN ;
          logLikelihoodNull = NaN ;
    else
%         dispersion0Tmp = mle(allGoals,'pdf',@(allGoals,r)nbinpdf(allGoals,dispersion0Tmp,dispersion0Tmp/(dispersion0Tmp+meanGoals)).^[weights(:);weights(:)]) ;
%         logLikelihoodSatured = sum(log(nbinpdf(allGoals,10^6,10^6/(allGoals+10^6))).*[weights(:);weights(:)]) ;
%         logLikelihoodNull = sum(log(nbinpdf(allGoals,dispersion0Tmp,pdispersion0Tmp/(dispersion0Tmp+meanGoals))).*[weights(:);weights(:)]) ;
          logLikelihoodSatured = NaN ;
          logLikelihoodNull = NaN ;
    end
elseif strcmpi(model,'gaussian')
    if isempty(weights) 
        sigma0Tmp = std(allGoals) ;
        logLikelihoodSatured = sum(log(normpdf(allGoals,allGoals,sigma0Tmp))) ;
        logLikelihoodNull = sum(log(normpdf(allGoals,meanGoals,sigma0Tmp))) ;
    else
        sigma0Tmp = sqrt(sum([weights(:);weights(:)] .* (allGoals - meanGoals).^2)) ; 
        logLikelihoodSatured = sum(log(normpdf(allGoals,allGoals,sigma0Tmp)).*[weights(:);weights(:)]) ;
        logLikelihoodNull = sum(log(normpdf(allGoals,meanGoals,sigma0Tmp)).*[weights(:);weights(:)]) ;
    end
elseif strcmpi(model,'cmp')
    if isempty(weights) 
        logLikelihoodSatured = NaN ;
        logLikelihoodNull = NaN ;
    else
        logLikelihoodSatured = NaN ;
        logLikelihoodNull = NaN ;
    end
end

deviance = 2 * (logLikelihoodSatured - logLikelihood) ;
devianceNull = 2 * (logLikelihoodSatured - logLikelihoodNull) ;

rSquared = 1 - (deviance / devianceNull) ;

% Max goal
maxGoal = max(allGoals) ;

out = struct('parameters',pList,'loglikelihood',logLikelihood,'nParEst',nParEst,...
        'AIC',AIC,'rSquared',rSquared,'allTeams',{allTeams},'nGames',nGames,...
        'model',model,'converged',converged,'maxGoal',maxGoal,...
        'fitter',fitter,'warningMsg',warningMsg) ;

end

