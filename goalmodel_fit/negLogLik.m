function negLL = negLogLik(params , goals1 , goals2 , team1 , team2 , hfa , model , pListSkeleton , weights)
% NEGLOGLIK Summary of this function goes here
%   Detailed explanation goes here

% Transformation double array -> structure
allTeams = unique([team1;team2]) ;
parCell{1} = params(1:numel(pListSkeleton.pTeam))' ;
parCell = [parCell;num2cell(params(numel(pListSkeleton.pTeam)+1:end)')] ;
a=cell2mat(parCell(1)) ;
parCell{1} = [a(1:end/2) a(end/2+1:end)] ;
pList = cell2struct(parCell,fields(pListSkeleton)) ;

% Add sum to zero constraint on parameters
pList.pTeam = [-sum(pList.pTeam); pList.pTeam] ;
pList.pTeam = array2table(pList.pTeam,'RowNames',allTeams,'VariableNames',{'attack','defence'}) ;

% Expectd goals (Poisson & nbin parameters)
expg = lambdaPred(pList , team1 , team2) ;

% Log-likelihood computation
if strcmpi(model,'poisson')
    logLik1 = log(poisspdf(goals1,expg.expg1)) ;
    logLik2 = log(poisspdf(goals2,expg.expg2)) ;
elseif strcmpi(model,'negbin')
    dispTmp = 1/exp(pList.dispersion) ;
    p1 = dispTmp ./ ( dispTmp + expg.expg1 ) ;
    p2 = dispTmp ./ ( dispTmp + expg.expg2 ) ;
    logLik1 = log(nbinpdf(goals1,dispTmp,p1)) ;
    logLik2 = log(nbinpdf(goals2,dispTmp,p2)) ;
elseif strcmpi(model,'gaussian')
    logLik1 = log(normpdf(goals1,expg.expg1,exp(pList.sigma))) ;
    logLik2 = log(normpdf(goals2,expg.expg2,exp(pList.sigma))) ;
elseif strcmpi(model,'cmp')
    expLogUpsilon = exp(pList.dispersion) ;   
    logLik1 = log(cmppdf(goals1,lambdaCMP(expg.expg1,expLogUpsilon,'fast'),expLogUpsilon)) ;
    logLik2 = log(cmppdf(goals2,lambdaCMP(expg.expg2,expLogUpsilon,'fast'),expLogUpsilon)) ;   
end

logLikTerms = logLik1 + logLik2 ;

% Dixon-Coles adjustement
if isfield(pList,'rho')
    dcAdj = tau(goals1,goals2,expg.expg1,expg.expg2,pList.rho) ;
    logLikTerms = logLikTerms + log(dcAdj) ;
end

if isempty(weights)
    logLik = sum(logLikTerms) ;
else
    logLik = sum(logLikTerms.*weights(:)) ;
end

negLL = -logLik ;

end

