function [D, D_all, facets, goodIdx] = genDepartureSet(m, n, scalar, seed)
% GENDEPARTURESET_ULTRAFAST
% If n <= 8, use full projection and real convex hull.
% If n > 8, use random top-norm selection without computing convex hull.

    if exist('seed','var') && ~isempty(seed)
        rng(seed);
    end

    fprintf('[Info] Running genDepartureSet_ultrafast with n = %d\n', n);

    % Step 1: Generate D
    D = randi(scalar, m, n);
    fprintf('[Step 1] Generated D of size [%d x %d]\n', size(D,1), size(D,2));

    if n <= 8
        %% Full projection version
        fprintf('[Mode] Full projection (n <= 8)\n');

        % Step 2: Generate all subset projections
        axesIdx = 1:n;
        combList = cell(n,1);
        for i = 1:n
            combList{i} = nchoosek(axesIdx, i);
        end
        rowsPer = sum(cellfun(@(C) size(C,1), combList));

        P_cell = cell(m,1);
        parfor i = 1:m
            local_points = zeros(rowsPer, n);
            local_cnt = 0;
            for dim = 1:n
                combos = combList{dim};
                for t = 1:size(combos,1)
                    local_cnt = local_cnt + 1;
                    vec = zeros(1,n);
                    idx = combos(t,:);
                    vec(idx) = D(i, idx);
                    local_points(local_cnt, :) = vec;
                end
            end
            P_cell{i} = local_points;
        end

        P_all = [vertcat(P_cell{:}); zeros(1, n)];
        D_all = [D; P_all];

        % Step 3: Convex hull
        fprintf('[Step 3] Computing convex hull...\n');
        t_start = tic;
        facets = convhulln(D_all);
        t_elapsed = toc(t_start);
        numFacets = size(facets,1);
        fprintf('[Info] Convex hull produced %d facets in %.2f seconds\n', numFacets, t_elapsed);

        % Step 4: Identify good facets
        valid = false(numFacets,1);
        parfor k = 1:numFacets
            verts_k = D_all(facets(k,:), :);
            zeroCounts = sum(verts_k == 0, 2);
            if sum(zeroCounts > 0) <= 1
                valid(k) = true;
            end
        end
        goodIdx = find(valid);
        fprintf('[Step 4] Found %d good facets\n', numel(goodIdx));

        if isempty(goodIdx)
            error('No facet satisfies the zero-component condition');
        end

    else
        %% Ultrafast northeast random sampling
        fprintf('[Mode] Ultrafast random sampling (n > 8)\n');

        % Step 2: Generate one-zero projections
        P_list = [];
        for i = 1:m
            for j = 1:n
                vec = D(i, :);
                vec(j) = 0;
                P_list = [P_list; vec];
            end
        end

        D_all = [D; P_list];
        fprintf('[Step 2] Generated D_all with projections, size [%d x %d]\n', size(D_all,1), size(D_all,2));

        % Step 3: Select top 10% northeast points
        norms = sum(D_all, 2);
        [~, sorted_idx] = sort(norms, 'descend');
        top_idx = sorted_idx(1:ceil(0.1 * numel(sorted_idx)));
        D_top = D_all(top_idx, :);
        fprintf('[Step 3] Selected top 10%% points, size [%d x %d]\n', size(D_top,1), size(D_top,2));

        % Step 4: Randomly sample pseudo-facets
        numFacets = 100; % target number of facets
        facets = zeros(numFacets, n);
        for f = 1:numFacets
            facets(f,:) = randsample(size(D_top,1), n);
        end
        goodIdx = (1:numFacets)';

        fprintf('[Step 4] Randomly sampled %d pseudo-facets from D_top\n', numFacets);
    end
end
