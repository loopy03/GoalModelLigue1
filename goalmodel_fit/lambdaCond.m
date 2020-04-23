function [ cond ] = lambdaCond( varargin )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
lambda = varargin{1} ;
mu = varargin{2} ;
upsilon = varargin{3} ;
if nargin==4
    ul = varargin{4} ;
else
    ul = 100 ;
end
cond = sum(exp((([0:ul] * log(lambda))) - (upsilon * log(factorial([0:ul])))) .* ([0:ul] - mu)) ;
end

