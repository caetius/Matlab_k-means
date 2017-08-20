% k-means_plot Matlab code.
% Diego Lorenzo-Casabuena Gonzalez, 2017. Prof. Dan Feldman
%%

% Plot results for k-line-means in 3D.
function [] = k_means_plot(P,k,j,mode)

[final_set,lowest_cost] = k_means(P,1,k,j,mode);

type = ['b-','g-','r-','b+','g+','r+']

temp = -10:1:10;
plot3(P(1,:),P(2,:),P(3,:), 'rx');
hold on;
for i = 1:k
    plot3(final_set(1,:,i)*temp,final_set(2,:,i)*temp,final_set(3,:,i)*temp, type(i))
end

end