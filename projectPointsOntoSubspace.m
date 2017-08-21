% k-means-by-coreset-reduction Matlab code.
% Algorithm to find the projection of matrix set of points on a subspace 

% Methodology: Given j+1 points defining a subspace in R^j, we shift all
% points and the subspace points to the orgin, QR the subspace points to
% get j orthonormal vectors defining the subspace and we take the
% projection of the translated points to the newly defined subspace.
% Finally, we take the distance between the translated points and the (also
% translated) projection.

% Input: Set of points and vectors defining subspace
% Output: Projection on subspace specified for all points
% Diego Lorenzo-Casabuena Gonzalez, 2017. Prof. Dan Feldman
%%

function [proj, basis, subspace_t, dist] = projectPointsOntoSubspace(points, subspace)

    translation = subspace(:,1);   
    subspace_t = bsxfun(@minus,subspace,translation); % Make j vectors from j+1 points
    subspace_t = subspace_t(:,2:size(subspace_t,2));
    points_t = bsxfun(@minus,points,translation);
    
    [basis,R] = qr(subspace_t,0); % Orthogonalise columns to form basis with QR
    proj = basis*basis'*points_t;
    dist = points-proj;
   
end