% k-means_plot Matlab code.
% Diego Lorenzo-Casabuena Gonzalez, 2017. Prof. Dan Feldman
%%

% Plot results for k-line-means in 3D.
function [] = k_means_plot(P,k,j,mode)

[final_set,lowest_cost] = k_means(P,k,j,mode);

plot3(P(1,:),P(2,:),P(3,:), 'bx')
hold on;
for i = 1:size(final_set,3)
    plot3(final_set(1,:,i),final_set(2,:,i),final_set(3,:,i))
end

end