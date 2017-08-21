% Find index of shortest projection of a point to an indexed set of Euclidean subspace
% Input: A point in R^d and 3D array of sets of (j+1)-d-dimensional points defining a j-subspace 
% Output: Indices of subspace closest to each point and the corresponding
% projection.
%Diego Lorenzo-Casabuena Gonzalez, 2017. Prof. Dan Feldman
%%

function [index, projection] = getIndexOfShortestProjection(point, subspaces,mode)

num_sub = size(subspaces,3);
current_score = intmax;

for i = 1:num_sub
   
    [proj, ~, ~,dist] = projectPointsOntoSubspace(point, subspaces(:,:,i));
    dist = findDistanceScore(dist, mode, 1);
    if dist < current_score
       current_score = dist;
       index = i;
       projection = proj;
    end
    
end

end