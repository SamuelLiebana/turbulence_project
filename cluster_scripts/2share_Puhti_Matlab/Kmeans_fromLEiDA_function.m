function mink = Puhti_leida_function


%% 2 - Cluster the Leading Eigenvectors
    %cd /users/tuulari/
    load('V1_all_new.mat')
    %disp('Clustering the eigenvectors into')
    % V1_all is a matrix containing all the eigenvectors:
    % Columns: N_areas are brain areas (variables)
    % Rows: (Tmax-2)*n_Sessions are all time points (independent observations)
    
    % USER: Set maximum/minimum number of clusters
    % There is no fixed number of FC states that the brain can display
    % Keep the range small for the first trials
    % Extend the range depending on the hypothesis of each work
    maxk=20;
    mink=2;
    rangeK=mink:maxk;
    
    % Set the parameters for Kmeans clustering
    Kmeans_results=cell(size(rangeK));
    
    for k=1:length(rangeK)
        % disp(['- ' num2str(rangeK(k)) ' FC states'])
        % Distance can be cosine, cityblock, default one is sqeuclidean
        [IDX, C, SUMD, D]=kmeans(V1_all,rangeK(k),'Distance','cosine','Replicates',200,'MaxIter',500,'Display','final','Options',statset('UseParallel',1));
        [~, ind_sort]=sort(hist(IDX,1:rangeK(k)),'descend');
        [~,idx_sort]=sort(ind_sort,'ascend');
        Kmeans_results{k}.IDX=idx_sort(IDX);   % Cluster time course - numeric column vectors
        Kmeans_results{k}.C=C(ind_sort,:);       % Cluster centroids (FC patterns)
        Kmeans_results{k}.SUMD=SUMD(ind_sort); % Within-cluster sums of point-to-centroid distances
        Kmeans_results{k}.D=D(:,ind_sort);       % Distance from each point to every centroid 
    end
    
    save Kmeans_results.mat
end 
