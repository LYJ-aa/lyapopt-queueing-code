% Extended evaluation: Performance at time T across multiple b values
clear; clc;

%% Parameters
n = 2;              % Number of queues
m = 3;              % Number of feasible service schedules
T = 1000;
runs = 1;
scenario = 'b_boundary';
scalar = 1;

%% b values setup
b_values = scalar * (2:140);
num_bs = numel(b_values);
result_max = zeros(1, num_bs);
result_lyap = zeros(1, num_bs);

for bi = 1:num_bs
    b = b_values(bi);
    % Define departure set and lambda for this b
    lambda = scalar * [1, (b-1)/b];
    D = scalar * [b, 0; 0, 1; 1, (b-1)/b];

    % Preallocate metrics for this b
    Metric_max_total = zeros(runs, T);
    Metric_lyap_total = zeros(runs, T);

    for run = 1:runs
        Q_max = zeros(n,1);
        Q_lyap = zeros(n,1);

        for t = 1:T
            % MaxWeight decision
            [~, idx_max] = max(D * Q_max);
            d_max = D(idx_max, :)';
            % Lyapunov decision
            lyap_costs = sum((max(Q_lyap' - D, 0)).^2, 2);
            [~, idx_lyap] = min(lyap_costs);
            d_lyap = D(idx_lyap, :)';

            % Generate arrivals
            A_curr = lambda';
      

            % Update queues
            Q_max = max(Q_max - d_max, 0) + A_curr;
            Q_lyap = max(Q_lyap - d_lyap, 0) + A_curr;

            % Record total queue lengths
            Metric_max_total(run, t) = sum(Q_max);
            Metric_lyap_total(run, t) = sum(Q_lyap);
        end
    end

    % Compute mean at final time T
    mean_max_total = mean(Metric_max_total, 1);
    mean_lyap_total = mean(Metric_lyap_total, 1);
    result_max(bi) = mean_max_total(end);
    result_lyap(bi) = mean_lyap_total(end);
end

%% Plot results

fig = figure('Units','inches', 'Position',[1,1,6,4]);
set(fig, 'PaperUnits','inches', 'PaperSize',[6,4], 'PaperPosition',[0,0,6,4]);
plot(b_values/sqrt(2), result_max, 'b-+', 'DisplayName', 'MaxWeight');
hold on;
plot(b_values/sqrt(2), result_lyap, 'r-<', 'DisplayName', 'LyapOpt');
xlabel('B');
ylabel(sprintf('Total queue length', T));
legend('Location', 'northwest');
title('Total queue length across different B');


set(fig,'color','w');
set(findall(fig, '-property', 'FontSize'), 'FontSize', 10);
export_fig(fig,'performance_vs_b.pdf', '-pdf', '-nocrop');


