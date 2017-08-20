% k-means-by-coreset-reduction Matlab code.
% Diego Lorenzo-Casabuena Gonzalez, 2017. Prof. Dan Feldman
%%


% Compute coreset and solution to k-line-means/median problem.
function [set,error] = k_means_by_coreset_reduction(P,k,j,mode)

    if j > 1 || j < 0
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
    
    tStart = tic;   % Start timer
    
    while(middle > 1)
       
       % 1) Compute L^* using k_means for a small subset of P
       point_subset = coreset(:,middle:last);
       [L,e] = k_means(point_subset,weights(middle:last),k,j,mode);
       
       % 2) Find projection for all points and assign each point to its
       % closest subspace (First column is always zeroes)
      
        data_points_proj_ind = zeros(2*input_dims(1)+1,size(point_subset,2));
        
       % Iterate through every point and assign it to its closest subspace
       for i = 1:size(point_subset,2)
           [index, proj] = getIndexOfShortestProjection(point_subset(:,i), L); 
           data_points_proj_ind(1:input_dims(1),i) = point_subset(:,i);
           data_points_proj_ind(input_dims(1)+1:2*input_dims(1),i) = proj;
           data_points_proj_ind(2*input_dims(1)+1,i) = index;
       end
      
       % 3) Compute importance for all points in each subspace (imp is 1xnum cols of point_subset vector)
       imp = zeros(1,1);
       points_final = zeros(input_dims(1),1);
       for i = 1:size(L,3)
           subspace_points = data_points_proj_ind(:,ismember(data_points_proj_ind(2*input_dims(1)+1,:),i));
           imp = [imp (importance_of_projection(subspace_points(input_dims(1)+1:2*input_dims(1),:), mode) + importance_of_points(subspace_points(1:input_dims(1),:),subspace_points(input_dims(1)+1:2*input_dims(1),:),mode))];
           points_final = [points_final subspace_points(1:input_dims(1),:)];
       end
       imp = imp(2:size(imp,2));  % Remove first index containing initial zero value
       points_final = points_final(:,2:size(points_final,2)); % Same for points_final
       
       % 4) Sample k/epsilon points from distribution imp and assign
       % weights
       indices = randperm(subproblem_size,sample_size);   % Selection currently not independent
       sample = points_final(:,indices);
       
       % 5) Update coreset
       last = middle + sample_size - 1;
       coreset(:,middle:last) = sample;
       weights(middle:last) = ones(1,(last+1)-middle)./imp(indices); % Set found weights
       middle = last - subproblem_size + 1;
       %Termination Case:
       if middle < 1
           middle = 1;
       end
       W = weights(middle:last);
       
       % Display percentage of progress
       perc = 100*(input_dims(2)-last)/input_dims(2);
       str = ['Completed: ', num2str(perc), '%'];
       clc;
       disp(str); % Display progress
       
    end

    [set,error] = k_means(coreset(:,1:last),weights,k,j,mode);
    clc;
    disp('Completed: 100.00%')
    tElapsed = toc(tStart); % Stop timer
    str = ['Time elapsed: ', num2str(tElapsed)];
    disp(str); %Display time taken to solve problem
end


% Importance of projection: Recursive algorithm assigns importance based on
% each projection's distance to the mean point of projection. Uses indices
% of initial array

function [P_imp_proj] = importance_of_projection(proj,mode)
    P_imp_proj = ones(1,size(proj,2))./size(proj,2);    % Initial importance
    center = sum(proj,2)/size(proj,2);                  % Mean projection
    dist = (sum(abs(bsxfun(@minus,proj,center)).^mode)).^(1/mode); % Distance of each projection to the mean
    [~,indices] = sort(dist,'descend');       % Sort distances and assign importance based on indices of sort operation
    
    while size(indices,2) > 1
        round((size(indices,2)+1)/2);
        indices = indices(round((size(indices,2)+1)/2):size(indices,2));   % Take |indices|/2 closest points to center
        P_imp_proj(indices) = ones(1,size(indices,2))./size(indices,2);    % Reassign weights
    end
end

% Calculate importance of points based on distance to projection / total
% projection distance
function [P_imp] = importance_of_points(points,proj,mode)
    P_imp = zeros(1,size(points,2));
    total_proj = sum(norm(sum(proj-points,2),mode));    % Sum of distances from all points to their projections
    for i = 1:size(points,2)
        P_imp(i) = sum(norm(proj(i)-points(i),mode))/total_proj;    % Importance = distance from point to projection / Sum of distances from all points to their projections
    end
end

