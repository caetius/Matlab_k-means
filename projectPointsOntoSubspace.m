% k-means-by-coreset-reduction Matlab code.
% Algorithm to find the projection of matrix set of points on a subspace 
% Input: Set of points and vectors defining subspace
% Output: Projection on subspace specified for all points
% Diego Lorenzo-Casabuena Gonzalez, 2017. Prof. Dan Feldman
%%

function [proj, basis, subspace_t] = projectPointsOntoSubspace(points, subspace)

    translation = subspace(:,1);   
    subspace_t = bsxfun(@minus,subspace,translation); % Make j vectors from j+1 points
    subspace_t = subspace_t(:,2:size(subspace_t,2));
    %points_t = bsxfun(@minus,points,translation);  % Shift points
    
   [basis,R] = qr(subspace_t,0); % Orthogonalise columns to form basis with QR 
   proj = basis*basis'*points;

   % Shift points and basis back from origin
  % basis = bsxfun(@plus,basis,translation);
   %proj = bsxfun(@plus,proj,translation);
 %  points = bsxfun(@plus,points,translation); % NOTE: Not necessary 
   
end