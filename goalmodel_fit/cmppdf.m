function [ prob ] = cmppdf(varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Deal with inputs
x = varargin{1} ;
lambda = varargin{2} ;
upsilon = varargin{3} ;
if nargin ==4
    error = varargin{4} ;
else
    error = 0.01 ;
end

cc = CMPNormalizingConstant(lambda,upsilon,error) ;

logProb = (x.*log(lambda)) - (upsilon.*log(factorial(x))) - log(cc) ;

prob = exp(logProb) ;

end

