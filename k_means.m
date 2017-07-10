% k-means Matlab code.
% Diego Lorenzo-Casabuena Gonzalez, 2017. Prof. Dan Feldman
%%

function [final_set,lowest_cost] = k_means(P,k,j,mode)

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Find all subsets of size k of the points set, P.
%       Uses basis vectors and an indexing from nchoosek
%%%%%%%%%%%%%%%%%%%%%%%%%%
[P_rows,P_cols] = size(P);
choice_j = reshape(nchoosek(1:P_cols,j+1)',[1,nchoosek(P_cols,j+1)*size(nchoosek(1:P_cols,j+1),2)]); % Get all subsets of length j as a vector of indices.
col_dim_j = (j+1)*nchoosek(P_cols,(j+1));
id_P = eye(P_cols); % Basis vectors
j_matrix = zeros(P_cols,col_dim_j);
for i = 1:col_dim_j
    j_matrix(:,i) = id_P(:,choice_j(i));    %Construct matrix from basis vectors. TODO: - Find a way to do this without for loops
end
set_of_jsubspaces = P*j_matrix; % Get all subspaces of size j using matrix of basis vectors
j_set_cols = size(set_of_jsubspaces,2); 

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Find all k combinations of the given j-subspaces in (1) -> to yield index vector.
%       Choose a chunk of k elements in choice_k to get each combinations of subspaces  
%%%%%%%%%%%%%%%%%%%%%%%%%%
num_j_subspaces = j_set_cols/(j+1); % Total number of subspaces possible from n points (the subspaces are gemerated above)
choice_k = reshape(nchoosek(1:num_j_subspaces,k)',[1,nchoosek(num_j_subspaces,k)*size(nchoosek(1:num_j_subspaces,k)',1)]); % Get all combinations of subspaces of length k (indices).
num_k_sets = size(choice_k,2) / k; % Number of candidate solutions


%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) For each element in sol_space evaluate objective function, store
%   results in a tensor
%%%%%%%%%%%%%%%%%%%%%%%%%%
%id_kj = eye(j_set_cols);
dist_tensor = zeros(P_rows,P_cols,num_j_subspaces); % Store disances between all points in P and each j-th subspace. Each set of distances is a matrix nxd
norms_matrix = zeros(P_cols,num_j_subspaces);
for i = 1:num_j_subspaces
    first_col = (j+1)*(i-1)+1; % Get next index of first element in subspace (each subspace is separated by j cols)
    subspace = set_of_jsubspaces(:,first_col:first_col+j); % Find current subspace
    translation = sum(subspace,2)/size(subspace,2);
    subspace_t = bsxfun(@minus,subspace,translation); % translate subspace so that it passes through the origin
    P_t = bsxfun(@minus,P,translation);
    % The idea here is to orthogonalise the columns of the subspace in
    % order to find the projection 
    [subspace_t,R] = qr(subspace_t,0); % Orthogonalise columns to form basis with QR
    dist_tensor(:,:,i) = (subspace_t*inv(subspace_t'*subspace_t)*subspace_t')*P_t - P_t; % Compute distances between all points and their projection in the ith j-subspace
    %% Find norms of all matrices 
    for j_ = 1:P_cols
        norms_matrix(i,j_) = norm(dist_tensor(:,j_,i).^mode);
    end
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
final_set = zeros(P_rows,j+1,k);
j_index = choice_k(:,k*(best_set-1)+1:k*best_set);
for i = 1:k
    final_set(:,:,i) = set_of_jsubspaces(:,j_index(i):j_index(i)+j);
end
% end fn
end

