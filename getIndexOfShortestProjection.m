% Find index of and projection of set of points to closest Euclidean subspace
% Input: Set of points in R^d and 3D array of sets of (j+1)-d-dimensional points defining a j-subspace 
% Output: Indices of subspace closest to each point and the corresponding
% projection.
%Diego Lorenzo-Casabuena Gonzalez, 2017. Prof. Dan Feldman
%%

function [index, projection] = getIndexOfShortestProjection(point, subspaces)

num_sub = size(subspaces,3);
current_score = intmax;

for i = 1:num_sub
   
    [proj, ~, ~] = projectPointsOntoSubspace(point, subspaces(:,:,i));
    dist = abs(point-proj);
    if dist < current_score
       current_score = dist;
       index = i;
       projection = proj;
    end
    
end

end