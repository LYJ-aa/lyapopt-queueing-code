function lambda = genLambda(D, D_top, facets, goodIdx, epsilon, seed, scenario)
% GENLAMBDA_FAST
% Generate a lambda vector based on the selected northeast facets.

    if exist('seed','var') && ~isempty(seed)
        rng(seed);
    else
        rng('shuffle');
    end

    n = size(facets, 2);

    switch scenario
        case 'b_boundary'
            % Sample uniformly from a random good facet
            fidx = goodIdx(randi(numel(goodIdx)));
            verts_idx = facets(fidx, :);
            verts = D_top(verts_idx, :);
            u = rand(n, 1); 
            u = u / sum(u);
            lambda = (u' * verts)';

        case 'b_epsilon'
            % Sample from facet and apply small perturbation
            fidx = goodIdx(randi(numel(goodIdx)));
            verts_idx = facets(fidx, :);
            verts = D_top(verts_idx, :);
            u = rand(n, 1);
            u = u / sum(u);
            lambda_b = (u' * verts)';
            epsVec = epsilon * (rand(n,1) + ones(n,1));
            lambda = max(lambda_b .* (1 - epsVec), 0);

        otherwise
            error('Unknown scenario: %s', scenario);
    end
end

