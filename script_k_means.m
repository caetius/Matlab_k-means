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
num_points = 200; % Number of points to generate for test input point set
max_value = 100; % Point set has range [0,max_value]
will_plot = true; % Plot results? Note: Can only plot 2D and 3D results for j=1

% Generate point set for the test 
P = rand(d,num_points)*max_value;
P = round(P); % Uncomment to run code on integer set of points

%%%%%%%%%%%%%%%%%%%% Code Test (Do not edit) %%%%%%%%%%%%%%%%%%%%

% Get results
[set,error] = k_means_by_coreset_reduction(P,k,j,mode);

% Plot results (Postponed)
% if will_plot && k < 8 && j == 1
% 
%     type = ['b*','g*','y*','m*','c*','w*','k*'];
%     temp = -max_value:1:max_value;
%     
%     if d == 3 
%         plot3(P(1,:),P(2,:),P(3,:), 'rx'); % Plot points
%         hold on;
%         for i = 1:k     % Plot subspaces (solution)
%             [~, final_set, ~, ] = projectPointsOntoSubspace(rand(d,1), set(:,:,i));
%             plot3(final_set(1,:)*temp,final_set(2,:)*temp,final_set(3,:)*temp, type(i));   
%         end
%     end
% 
%     if d == 2
%         plot(P(1,:),P(2,:), 'rx'); % Plot points
%         hold on;
%         for i = 1:k     % Plot subspaces (solution)
%             [~, final_set, ~] = projectPointsOntoSubspace(rand(d,1), set(:,:,i));
%             plot(final_set(1,:)*temp,final_set(2,:)*temp, type(i));
%         end
%     else
%         msg = 'Can only plot 2- and 3-D results';
%         err(msg);
%     end     
% end