% Find vector of distances between two sets of points
% Input: Set of distances (mxn), a mode (1=median, 2=mean), and weights for
% each score.
% Output: Vector of weighted distances (1xn)
% Diego Lorenzo-Casabuena Gonzalez, 2017. Prof. Dan Feldman
%%
function [scores] = findDistanceScore(dist, mode, weights)

scores = ((sum(abs(dist).^mode)).^(1/mode)).*weights;

end