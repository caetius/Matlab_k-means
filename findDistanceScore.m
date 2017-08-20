% Find vector of distances between two sets of points
% Input: Two sets of points (mxn) and a mode (1=median, 2=mean)
% Output: Vector of weighted distances (1xn)
% Diego Lorenzo-Casabuena Gonzalez, 2017. Prof. Dan Feldman
%%
function [scores] = findDistanceScore(points1, points2, mode, weights)

if size(points1) ~= size(points2)
   msg = 'Dimensions dont match';
   err(msg);
end

scores = ((sum(abs(points1-points2).^mode)).^(1/mode)).*weights;

end