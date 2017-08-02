% k-means Matlab code.
% Diego Lorenzo-Casabuena Gonzalez, 2017. Prof. Dan Feldman
%%

function [final_set,lowest_cost] = k_means(P,W,k,j,mode)

% First, check the input weights vector. If program is used without coresets, weights 
% are set to 1 by default. Else, custom weights are used.
[P_rows,P_cols] = size(P);
if size(W) ~= P_cols
    W = ones(1,P_cols);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Find all subsets of size k of the points set, P.
%       Uses basis vectors and an indexing from nchoosek
%%%%%%%%%%%%%%%%%%%%%%%%%%

choice_j = reshape(nchoosek(1:P_cols,j+1)',[1,nchoosek(P_cols,j+1)*size(nchoosek(1:P_cols,j+1),2)]); % Get all subsets of length j as a vector of indices.
id_P = eye(P_cols); % Basis vectors
j_matrix = id_P(:,choice_j);  %Construct matrix from basis vectors
set_of_jsubspaces = P*j_matrix; % Get all subspaces of size j using matrix of basis vectors
j_set_cols = size(set_of_jsubspaces,2); 

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Find all k combinations of the given j-subspaces in (1) to yield index vector.
%       Choose a chunk of k elements in choice_k to get each combinations of subspaces  
%%%%%%%%%%%%%%%%%%%%%%%%%%

num_j_subspaces = j_set_cols/(j+1); % Total number of subspaces possible from n points (the subspaces are gemerated above)
choice_k = reshape(nchoosek(1:num_j_subspaces,k)',[1,nchoosek(num_j_subspaces,k)*size(nchoosek(1:num_j_subspaces,k)',1)]); % Get all combinations of subspaces of length k (indices).
num_k_sets = size(choice_k,2) / k; % Number of candidate solutions


%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) For each element in sol_space evaluate objective function, store
%   results in a tensor
%%%%%%%%%%%%%%%%%%%%%%%%%%

norms_matrix = zeros(num_j_subspaces,P_cols);
for i = 1:num_j_subspaces
 
    first_col = (j+1)*(i-1)+1; % Get next index of first element in subspace (each subspace is separated by j cols)
    subspace = set_of_jsubspaces(:,first_col:first_col+j); % Find current subspace
    [proj, basis, subspace_t] = projectPointsOntoSubspace(P, subspace)
    norms_matrix(i,:) = findDistanceScore(P, proj, mode, W);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4) Find sets that minimize objective function given matrix of distance
% norms for all subspaces and points
%%%%%%%%%%%%%%%%%%%%%%%%%%

best_cols = zeros(1,P_cols);
best_set = 0;
lowest_cost = intmax;
for i = 1:num_k_sets     % For every set of subspaces
    indices = choice_k(k*(i-1)+1:k*i);    % Retrieve indices of subspaces in each set
    for col = 1:P_cols     % For every subspace in indices, find matrix with min distances to all points in P and store set and norm if
        best_cols(col) = min(norms_matrix(indices,col));    % smallest so far
    end
    current_cost = norm(best_cols); 
    if current_cost < lowest_cost
        lowest_cost = current_cost;
        best_set = i;
    end  
end
    

% Get set of j-subspaces (and the subspaces themselves) at index best_set
% Subspaces are given as basis of vectors
final_set = zeros(P_rows,j+1,k);
j_index = choice_k(:,k*(best_set-1)+1:k*best_set);
for i = 1:k
    size(set_of_jsubspaces(:,j_index(i):j_index(i)+j))
    size(final_set)
    final_set(:,:,i) = set_of_jsubspaces(:,j_index(i):j_index(i)+j);
    %my_set = set_of_jsubspaces(:,j_index(i):j_index(i)+j);
    %translation = my_set(:,1);
    %my_set_t = bsxfun(@minus,my_set,translation); % Make j vectors from j+1 points
    %[Q,R] = qr(my_set_t(:,2:size(my_set_t,2)),0);
    %final_set(:,:,i) = bsxfun(@plus,Q,translation);
    
end
% end fn
end
