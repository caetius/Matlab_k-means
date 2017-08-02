
% Algorithm to find the projection of matrix set of points on a subspace 
% Input: Set of points and vectors defining subspace
% Output: Projection on subspace specified for all points

function [proj, basis, subspace_t] = projectPointsOntoSubspace(points, subspace)

    translation = subspace(:,1);   
    subspace_t = bsxfun(@minus,subspace,translation) % Make j vectors from j+1 points
    subspace_t = subspace_t(:,2:size(subspace_t,2))
    
   [basis,R] = qr(subspace_t,0) % Orthogonalise columns to form basis with QR
    proj = basis*basis'*points

end