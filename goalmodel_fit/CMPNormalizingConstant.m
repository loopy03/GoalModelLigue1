function [ cc ] = CMPNormalizingConstant( varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Deal with inputs
lambda = varargin{1} ;
upsilon = varargin{2} ;
if nargin ==3
    error = varargin{3} ;
else
    error = 0.01 ;
end

% Upper limit of the sum

ul = max([50*ones(length(lambda),1),ceil((lambda./error).^(1./upsilon))],[],2) ;
ul = 50*ones(length(lambda),1) ;

% Compute the normalizing constant
cc = zeros(length(lambda),1) ;
for ii=1:length(lambda)
    cc(ii) = sum(exp(([0:ul(ii)] * log(lambda(ii))) - (upsilon * log(factorial([0:ul(ii)]))))) ; 
end

end

