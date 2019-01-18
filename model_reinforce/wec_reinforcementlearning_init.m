%% Wave energy converter simulation
% Created: ECE 530, 7 Nov 2017
% Sub optimal PTO damping = 1.1364e+06
clc
close all
clear
format compact



%% WEC parameters

A = 3.3e5;  % added mass kg
B = 9.7e4;  % damping N/(m/s)
K = 8.8e5;  % hydrostatic stiffness (bouyancy) N/m
m = 9.7e4;  % dry mass kg

% Fastest period ~ 1 second -> 2*pi rad/sec
wp = (2*pi/1)*1;  % pole for derivative approximation

zb_0 = 0;   % initial cond
dzb_0 = 0;  % initial cond



%% Load water surface elevation time series
load eta

t = [0:length(eta.values)-1]*eta.sampleTime;
t = t';
zw_timeseries = [t eta.values];

% Gently phase in the waves to avoid abrupt initial derivatives
mask1 = [0:0.01:1]';
remainingLength = length(zw_timeseries(:,2)) - length(mask1);
mask = [mask1 ; ones(remainingLength,1)];
zw_timeseries(:,2) = zw_timeseries(:,2) .* mask;



%% PTO
% Can define Bpto here



%% Clean up
disp('Great job!  You have done it again!!')
