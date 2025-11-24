%% Parameters
n = 5;
m = n*10;
T = 1000;
runs = 100;
num_lambdas = 2000;
scenario = 'b_boundary';
sigma_A = 1 * ones(n,1);
sigma_D = 1 * ones(n,1);
epsilon = 0.1;
scalar = 10;
alpha = 0.05;   % for final proportion CI
eps = 1e-6;

%% Logs (total metric)
better_log              = false(num_lambdas, 1);   % Lyap <= Max (total)
moderate_gain_log_total = false(num_lambdas, 1);   % >=0.1
big_gain_log_total      = false(num_lambdas, 1);   % >=0.2
bigg_gain_log_total     = false(num_lambdas, 1);   % >=0.3
Big_gain_log_total      = false(num_lambdas, 1);   % >=0.4
large_gain_log_total    = false(num_lambdas, 1);   % >=0.5

%% Generate departure set
[D, D_top, facets, goodIdx] = genDepartureSet(m, n, scalar, 42);

%% (Optional) Parallel setup
if isempty(gcp('nocreate')), parpool; end

%% Main simulation
for lambda_trial = 1:num_lambdas

    lambda = genLambda(D, D_top, facets, goodIdx, epsilon, [], scenario);

    final_max_total  = zeros(runs, 1);
    final_lyap_total = zeros(runs, 1);

    for run = 1:runs
        Q_max  = zeros(n,1);
        Q_lyap = zeros(n,1);

        for t = 1:T

            % MaxWeight decision
            [~, idx_max] = max(D * Q_max);
            d_max = D(idx_max, :)';

            % Lyapunov decision
            lyap_costs = sum((max(Q_lyap' - D, 0)).^2, 2);
            [~, idx_lyap] = min(lyap_costs);
            d_lyap = D(idx_lyap, :)';

            %% arrivals

            A_curr = zeros(n,1);
            for i = 1:n
                mu = lambda(i);
                v  = sigma_A(i)^2;

                if mu <= 0
                    A_curr(i) = 0;
                    continue;
                end

                if v < 0
                    error('Negative variance detected for arrivals at i = %d.', i);
                elseif v == 0
                    A_curr(i) = mu;
                else
                    % Gamma
                    k     = mu^2 / v;   % shape
                    theta = v / mu;     % scale
                    A_curr(i) = gamrnd(k, theta);
                end
            end

            % --- MaxWeight departures ---
            for i = 1:n
                if d_max(i) ~= 0
                    mu = d_max(i);
                    v  = sigma_D(i)^2;

                    if mu <= 0
                        d_max(i) = 0;
                        continue;
                    end

                    if v < 0
                        error('Negative variance detected for MaxWeight departure at i = %d.', i);
                    elseif v == 0
                        d_max(i) = mu;
                    else
                        % Gamma
                        k     = mu^2 / v;   % shape
                        theta = v / mu;     % scale
                        d_max(i) = gamrnd(k, theta);
                    end
                end
            end

            % --- LyapOpt departures ---
            for i = 1:n
                if d_lyap(i) ~= 0
                    mu = d_lyap(i);
                    v  = sigma_D(i)^2;

                    if mu <= 0
                        d_lyap(i) = 0;
                        continue;
                    end

                    if v < 0
                        error('Negative variance detected for LyapOpt departure at i = %d.', i);
                    elseif v == 0
                        d_lyap(i) = mu;
                    else
                        % Gamma
                        k     = mu^2 / v;   % shape
                        theta = v / mu;     % scale
                        d_lyap(i) = gamrnd(k, theta);
                    end
                end
            end


            % Queue update
            Q_max  = max(Q_max  - d_max, 0) + A_curr;
            Q_lyap = max(Q_lyap - d_lyap, 0) + A_curr;
        end

        % Only final-time total metric
        final_max_total(run)  = sum(Q_max);
        final_lyap_total(run) = sum(Q_lyap);
    end

    % Average over runs
    avg_max_total  = mean(final_max_total);
    avg_lyap_total = mean(final_lyap_total);

    % Lyap better
    if avg_lyap_total <= avg_max_total
        better_log(lambda_trial) = true;
    end

    % Relative improvement
    rel_imp = (avg_max_total - avg_lyap_total) / avg_max_total;
    if rel_imp >= 0.1
        moderate_gain_log_total(lambda_trial) = true;
    end
    if rel_imp >= 0.2
        big_gain_log_total(lambda_trial) = true;
    end
    if rel_imp >= 0.3
        bigg_gain_log_total(lambda_trial) = true;
    end
    if rel_imp >= 0.4
        Big_gain_log_total(lambda_trial) = true;
    end
    if rel_imp >= 0.5
        large_gain_log_total(lambda_trial) = true;
    end
end

%% Save summary only
save('lambda_eval_summary.mat', ...
    'better_log', ...
    'moderate_gain_log_total', ...
    'big_gain_log_total', ...
    'bigg_gain_log_total', ...
    'large_gain_log_total');

%% Print summary with 95% CI
fprintf('---\n');
N = num_lambdas;

print_with_error('Lyapunov better than MaxWeight:',   sum(better_log), N);
print_with_error('Rel. improvement total ≥ 0.1:',     sum(moderate_gain_log_total), N);
print_with_error('Rel. improvement total ≥ 0.2:',     sum(big_gain_log_total), N);
print_with_error('Rel. improvement total ≥ 0.3:',     sum(bigg_gain_log_total), N);
print_with_error('Rel. improvement total ≥ 0.4:',     sum(Big_gain_log_total), N);
print_with_error('Rel. improvement total ≥ 0.5:',     sum(large_gain_log_total), N);


%% ---- helper function ----
function print_with_error(label, count, N)
if N <= 0
    fprintf('%-40s N is zero.\n', label);
    return;
end
p = count / N;
z = 1.96;
se = sqrt(p*(1-p)/N);
hw = z * se;
fprintf('%-40s %.2f%% ± %.2f%% (95%% CI)\n', ...
    label, 100*p, 100*hw);
end
