% Script to test k-means code using coresets and plotting results
% NOTE: Edit input variables where indicated. Do not edit the rest of the
% code.
%Diego Lorenzo-Casabuena Gonzalez, 2017. Prof. Dan Feldman.
%%

clear;
close all;
clc;

%%%%%%%%%%%%%%%%%%%% Input variables (Edit) %%%%%%%%%%%%%%%%%%%%

mode = 2; % 1 for k-median, 2 for k-means problem
k = 3; % Number of j-subspaces to solve for
j = 1; % Type of subspace (0 = points, 1 = lines, ...) Note: Coreset code only available for j = 0,1
d = 2; % Dimension of input point set 
num_points = 20; % Number of points to generate for test input point set
max_value = 100; % Point set has range [0,max_value]
will_plot = true; % Plot results? Note: Can only plot 2D and 3D results for j=1

% Generate point set for the test 
P = rand(d,num_points)*max_value;
P = round(P); % Uncomment to run code on integer set of points

use_coresets = false; % If true, the problem is sped up with coresets

%%%%%%%%%%%%%%%%%%%% Code Test (Do not edit) %%%%%%%%%%%%%%%%%%%%

% Get results
if use_coresets
    [set,error] = k_means_by_coreset_reduction(P,k,j,mode);
else
    [set,error] = k_means(P,1,k,j,mode);
end

type = ['b*','g*','m*','y*','b+','g+','m+','y+'];


%Plot results (Postponed)
if will_plot && k < 9 && j == 1
    
    if d == 3 
        plot3(P(1,:),P(2,:),P(3,:), 'rx'); % Plot points
        hold on;
        for i = 1:k     % Plot subspaces (solution)
            plot3(final_set(1,:,i),final_set(2,:,i),final_set(3,:,i), type(i));   
        end
    end

    if d == 2
        plot(P(1,:),P(2,:), 'rx'); % Plot points
        hold on;
        for i = 1:k     % Plot subspaces (solution);
            plot(set(1,:,i),set(2,:,i), type(i));
        end
    else
        msg = 'Can only plot 2-D and 3-D results';
        err(msg);
    end     
end