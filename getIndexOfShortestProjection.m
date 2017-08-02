function [index, projection] = getIndexOfShortestProjection(point, subspaces)

num_sub = size(subspaces,3);
current_score = intmax;

for i = 1:num_sub
   
    [proj, basis, subspace_t] = projectPointsOntoSubspace(point, subspaces(:,:,i));
    dist = abs(point-proj);
    if dist < current_score
       current_score = dist;
       index = i;
       projection = proj;
    end
    
end

end