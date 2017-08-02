% k-means-by-coreset-reduction Matlab code.
% Diego Lorenzo-Casabuena Gonzalez, 2017. Prof. Dan Feldman
%%


% Compute coreset and solution to k-line-means/median problem.
function [set,error] = k_means_by_coreset_reduction(P,k,j,mode)

    if j ~= 1
       msg = 'This code is only valid for k-point- and k-line-means. (j=0,1)';
       error(msg);
    end

    coreset = P;    % Starts off containing all coreset values, and by the end of the algorithm the first sample_size values are the final coreset
    input_dims = size(coreset); % Size of input points
    subproblem_size = 10;    % Algorithm divided into subproblems of this size (approx. binary tree)
    sample_size = round(log(subproblem_size));    % Final size of coreset. Also size of output at each subproblem
    last = input_dims(2);
    middle = last - subproblem_size + 1;
    weights = ones(1,last); % Subsets of the weights are stored at these indices on every iteration
    W = 0;
    
    while(last > subproblem_size)
       
       % 1) Compute L^* using k_means for a small subset of P
       tree_subset = coreset(:,middle:last);
       [L,e] = k_means(tree_subset,weights(middle:last),k,j,mode);
       
       % 2) Find projection for all points and assign each point to its
       % closest subspace (First column is always zeroes)
       points_by_subspace = zeros(input_dims(1),1,size(L,3));
       proj_by_subspace = zeros(input_dims(1),1,size(L,3));
       
       assignment_indices = zeros(1,size(tree_subset,2));
       
       % Iterate through every point and assign it to its closest subspace
       for i = 1:size(tree_subset,2)
           [index, proj] = getIndexOfShortestProjection(tree_subset(:,i), L);
           assignment_indices(i) = index;
           if proj_by_subspace(:,:,index) == zeros(input_dims(1),1)
               proj_by_subspace(:,:,index) = proj;
               points_by_subspace(:,:,index) = tree_subset(:,i);
           else
               horzcat(proj_by_subspace(:,:,index),proj)
               proj_by_subspace(:,:,index) = horzcat(proj_by_subspace(:,:,index),proj);
               points_by_subspace(:,:,index) = horzcat(points_by_subspace(:,:,index),tree_subset(:,i));
           end
       end

       % 3) Compute importance for all points in each subspace (imp is 1xsample_size vector)
       for i = 1:size(L,3)
            imp = [imp importance_of_projection(point_projections, mode) + importance_of_points(tree_subset,point_projections)];
       end
       
       % 4) Sample k/epsilon points from distribution imp and assign
       % weights
       indices = round(size(tree_subset,2)*rand(1,sample_size));
       points_final = reshape(points_by_subspace,size(points_by_subspace,1),size(points_by_subspace,2)*size(points_by_subspace,3));
       sample = points_final(:,indices);
       
       % 5) Update coreset
       last = last - (subproblem_size-size(sample,2)+1);
       coreset(:,middle:last) = sample;
       weights(middle:last) = ones(middle:last)./imp(indices); % Set found weights
       middle = last - subproblem_size;
       W = weights(middle:last);
       
       disp('Completed: ');
       disp((input_dims(2)-last)/input_dims(2)); % Display progress
       
    end

    [set,error] = k_means(coreset(1:subproblem_size),weights,k,j,mode);

end


% Importance of projection: Recursive algorithm assigns importance based on
% each projection's distance to the mean point of projection. Uses indices
% of initial array

function [P_imp_proj] = importance_of_projection(proj,mode)
    P_imp_proj = ones(1,size(proj,2))./size(proj,2);    % Initial importance
    center = sum(proj,2)/size(proj,2);                  % Mean projection
    dist = (sum(abs(bsxfun(@minus,proj,center)).^mode)).^(1/mode); % Distance of each projection to the mean
    [sorted_list,indices] = sort(dist,'descend');       % Sort distances and assign importance based on indices of sort operation
    
    while size(indices) > 0
        indices = indices(size(indices)/2:size(indices));   % Take |indices|/2 closest points to center
        P_imp_proj(indices) = ones(1,size(indices))./size(indices);    % Reassign weights
    end
end

% Calculate importance of points based on distance to projection / total
% projection distance
function [P_imp] = importance_of_points(points,proj,mode)
    P_imp = zeros(core_dims(1),size(points,2));
    total_proj = sum(norm(sum(proj-points,2),mode));    % Sum of distances from all points to their projections
        
    for i = 1:size(points)
        P_imp(i) = sum(norm(proj(i)-points(i),mode))/total_proj;    % Importance = distance from point to projection / Sum of distances from all points to their projections
    end
end

