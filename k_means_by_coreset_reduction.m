% k-means-by-coreset-reduction Matlab code.
% Diego Lorenzo-Casabuena Gonzalez, 2017. Prof. Dan Feldman
%%


% Compute coreset and solution to k-line-means/median problem.
function [set,error] = k_means_by_coreset_reduction(P,k,j,mode,epsilon)

    if j > 1
       msg = 'This code is only valid for k-point- and k-line-means. (j=0,1)';
       error(msg);
    end

    coreset = P;
    core_dims = size(coreset);
    sample_size = 20;
    last = core_dims(2);
    middle = last - sample_size;
    weights = ones(1,last); % Subsets of the weights are stored at these indices on every iteration
    W = 0;
    while(last > sample_size)
       
       % 1) Compute L^* using k_means for a small subset of P
       tree_subset = coreset(:,middle:last);
       [L,e] = k_means(tree_subset,weights(middle:last),k,j,mode);
       
       % 2) Find projection ::: TODO: - How to do this in multiple lines. 
       % For now assume one line
       translation = sum(tree_subset,2)/size(tree_subset,2);
       tree_subset_t =  bsxfun(@minus,tree_subset,translation);
       proj = zeros(size(tree_subset),k);
       basis = qr(bsxfun(@minus,L,translation),0);
       compare_matrix = 
       for i = 1:core_dim(2)
            for j_ = 1:k
                
                points_by_subspace(:,:,j_) = horzcat(projection(:,:,i),tree_subset_t(i));
            
              
       end
       (basis*inv(basis'*basis)*basis')*tree_subset_t;
       % 3) Compute importance for all p in tree_subset (imp is 1xsample_size vector)
       imp = importance_of_projection(P_) + importance_of_points(tree_subset,L);
       
       % 4) Sample k/epsilon points from distribution imp and assign
       % weights
       indices = round(tree_subset*rand(1,k/epsilon));
       sample = tree_subset(:,indices);
       
       % 5) Update coreset
       last = last - (sample_size-size(sample,2));
       coreset(:,middle:last) = sample;
       weights(middle:last) = 1./imp(indices); % Set found weights
       middle = last - sample_size;
       W = weights(middle:last);
    end

    [set,error] = k_means(coreset,weights,k,j,mode);

end


% Total importance given by importance_of_projection + importance_of_points
% Work with indices of original matrix

function [P_imp_proj] = importance_of_projection(P_)
    P_imp_proj = ones(1,size(P_,2))./size(P_,2); % Initial importance
    while size(P_,2) > 0
        center = sum(P_,2)/size(P_,2);      % Center of projected points
        dist = bsxfun(@minus,P_,center);    % Find distances to center
        norms = zeros(1,size(dist,2));      
        for i=0:size(dist,2)                
            norms(i) = norm(dist(:,i));     % Get distance norms
        end
        small_ind = sort(norms,2,'ascend'); % Sort norms
        % Must get index of larger array here: maybe from previous
        % iteration and setting an initial permutation of size n.
        l_ind = small_ind; % TODO: - Read above
        P_ = P_(small_ind);                 % Set P_ to n/2 points in P_ closest to center
        for i = 1:size(P_,2)                
             P_imp_proj(l_ind);                  % Update importance of remaining elements
        end
        
    end

end

function [P_imp] = importance_of_points(P,L)
    %% could take this from k_means code

end

