function [ approx ] = lambdaApprox( mu , upsilon )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
approx = (mu + ((upsilon - 1) ./ (2.*upsilon))) .^ upsilon ;
end

