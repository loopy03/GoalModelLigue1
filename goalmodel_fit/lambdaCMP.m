function [ res ] = lambdaCMP( varargin )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
mu = varargin{1} ;
upsilon = varargin{2} ;
if nargin>=3
    method = varargin{3} ;
else
    method = 'sum' ;
end
if nargin==4
    error = varargin{4} ;
else
    error = 0.01 ;
end

lambdas = lambdaApprox(mu,upsilon) ;

if strcmpi(method,'fast')
    res = lambdas ;
elseif strcmpi(method,'sum')
    ul = max([50*ones(length(lambdas),1),ceil((lambdas./error).^(1./upsilon))],[],2) ;
    
    for ii=1:length(lambdas)
        lambdaInt = max([0.001*ones(2,1) , [-5;5] + lambdas(ii)*ones(2,1)],[],2) ;
            
        % Upper extending
        fun = @(lambda) lambdaCond(lambda,mu(ii),upsilon,ul(ii)) ;
        while fun(lambdaInt(2))*fun(lambdaInt(1))>=0
            lambdaInt(2) = lambdaInt(2) + 0.1*(lambdaInt(2)-lambdaInt(1)) ;
        end
       
        % Solve
        rootRes = fzero(fun,lambdaInt) ;
        lambdas(ii) = rootRes ;
    end
    res = lambdas ;
end

end

