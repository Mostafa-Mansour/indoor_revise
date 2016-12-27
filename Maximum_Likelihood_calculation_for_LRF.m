% function [ Maxliklehood ] = Maximum_Likelihood_calculation_for_LRF( Y,S,R )

function [ Maxliklehood ] = Maximum_Likelihood_calculation_for_LRF( Y,sx,R )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%% Y is the vector of our measurements 
%%% R is the variance of the measurements error
%%% H is the measurement matrix
%%% X is the state vector 

m=length(Y); % number of our meaurements 

% sum = 0; % for calculating the summution in the exponential part   
% mul= 1; % for calculating the multiplaction part o the equation (denominator)

% for i=1:m
%     
%     sum= sum + (((Y(i)- (y/cosd(k+ i*phi)))^2) /r);
%     
%     mul = mul * ((2*pi)^(m/2)) * r^m ;
%     
%     
% end

% sum = ((Y-S)^2) /R ;
sum = (Y-sx)'*(R^-1)*(Y-sx);
Maxliklehood = (exp(-0.5 * sum)); 

end

