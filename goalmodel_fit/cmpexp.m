function [ ee ] = cmpexp( varargin )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% Deal with inputs
lambda = varargin{1} ;
upsilon = varargin{2} ;
if nargin ==4
    error = varargin{4} ;
else
    error = 0.01 ;
end
if nargin>=3
    method = varargin{3} ;
else
    method = 'sum' ;
end

if strcmpi(method,'fast')
    ee = lambda.^(1./upsilon) - ((upsilon - 1) ./ (2*upsilon)) ;
elseif strcmpi(method,'sum')
    % upper limit for computing the normalizing constant
    ul = ceil((lambda ./ error).^(1./upsilon)) ;
    
    % normalizing constant
    cc = CMPNormalizingConstant(lambda,upsilon,error) ;
    
    % expectation
    ee = zeros(length(lambda),1) ;
    for ii=1:length(lambda)
        ee(ii) = sum(exp((([0:ul] * log(lambda(ii)))) - (upsilon * log(factorial([0:ul]))) - log(cc(ii))) .* [0:ul]) ;
    end
end

end

