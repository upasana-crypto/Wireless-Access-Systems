clear all;
close all;
clc;

%% multi-path channel
L = 10;
lamda = 1/5;

% draw random coefficients from Gaussian distribution
h_coeff = [1+1i,2+3i];

% construct h
h = zeros(1,length(h_coeff));
for ind_h = 1:length(h_coeff)
    h(ind_h)= h(ind_h)+h_coeff(ind_h)*exp(-lamda*(ind_h-1));    
end
P_h = power_calculation(h);

P_temp_h = power_calculation(temp_h);
%% 

P_h_wanted = 1;

% construct h_scaled
h_sc = zeros(size(h));
for ind_h_sc=1: length(h)
    h_sc(ind_h_sc) = sqrt(P_h_wanted/P_h) * h(ind_h);
end

P_h_sc = power_calculation(h_sc);

