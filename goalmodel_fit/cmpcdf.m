function [ pp ] = cmpcdf(varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% Deal with inputs
x = varargin{1} ;
lambda = varargin{2} ;
upsilon = varargin{3} ;
if nargin ==5
    error = varargin{5} ;
else
    error = 0.01 ;
end
if nargin>=4
    lower = varargin{4} ;
else
    lower = true ;
end

cc = CMPNormalizingConstant(lambda,upsilon,error) ;

pp = zeros(length(lambda),1) ;
for ii=1:length(lambda)
    pp(ii) = sum(exp( ([0:x(ii)] * log(lambda(ii))) - (upsilon * log(factorial([0:x(ii)]))) - log(cc(ii)) )) ; 
end

if lower==false
    pp = 1 - pp ;
end
end

