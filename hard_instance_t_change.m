

clear; clc;

%% Parameters
n = 2;              % Number of queues
m = 3;              % Number of feasible service schedules
T = 2000;
runs = 1;
sigma_A = 0 * ones(n,1);   % Arrival variance for each queue
sigma_D = 0 * ones(n,1);   % Service variance for each queue
scalar = 1;
alpha = 0.05;  % 95% confidence
b = scalar*10*sqrt(2);
%% Generate D and facets
D = scalar*[b,0;0,1;1,(b-1)/b];


%% Generate one lambda
lambda = scalar*[1,(b-1)/b];

%% Preallocate
Metric_max_sq = zeros(runs, T);
Metric_lyap_sq = zeros(runs, T);
Metric_max_total = zeros(runs, T);
Metric_lyap_total = zeros(runs, T);


%% Run simulations
for run = 1:runs
    Q_max = zeros(n,1);
    Q_lyap = zeros(n,1);

    for t = 1:T
        [~, idx_max] = max(D * Q_max);
        d_max = D(idx_max, :)';

        lyap_costs = sum((max(Q_lyap' - D, 0)).^2, 2);
        [~, idx_lyap] = min(lyap_costs);
        d_lyap = D(idx_lyap, :)';

        tol = 1e-8;   % Tolerance used to determine whether v = mu

        %% arrivals
        if all(sigma_A == 0)
            % deterministic arrivals
            A_curr = lambda';
        else
            % randomized arrivals
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
        end
       

        %% departures
        if all(sigma_D == 0)
            % deterministic departures
            d_max  = d_max;
            d_lyap = d_lyap;
        else
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
        end



        Q_max = max(Q_max - d_max, 0) + A_curr;
        Q_lyap = max(Q_lyap - d_lyap, 0) + A_curr;


        Metric_max_total(run, t)  = sum(Q_max);
        Metric_lyap_total(run, t) = sum(Q_lyap);
    end
end

%% Post-processing
df = runs - 1;
tcrit = tinv(1 - alpha/2, df);
x = 1:T;

mean_max_total = mean(Metric_max_total,1);
mean_lyap_total = mean(Metric_lyap_total,1);
se_max_total = std(Metric_max_total,0,1)/sqrt(runs);
se_lyap_total = std(Metric_lyap_total,0,1)/sqrt(runs);
hw_max_total = tcrit * se_max_total;
hw_lyap_total = tcrit * se_lyap_total;


%% Plot Total Queue Length with 90% CI
fig2 = figure('Units','inches', 'Position',[1,1,6,4]);
set(fig2, 'PaperUnits','inches', 'PaperSize',[6,4], 'PaperPosition',[0,0,6,4]);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 10);

ax = gca;
hold(ax, 'on');   


% 计算上下界
max_upper = mean_max_total + hw_max_total;
max_lower = mean_max_total - hw_max_total;
lyap_upper = mean_lyap_total + hw_lyap_total;
lyap_lower = mean_lyap_total - hw_lyap_total;

% 阴影带 (MaxWeight)
fill([x fliplr(x)], [max_upper fliplr(max_lower)], ...
     [0 0.4470 0.7410], 'EdgeColor','none','FaceAlpha',0.3,'HandleVisibility','off');

% 阴影带 (LyapOpt)
fill([x fliplr(x)], [lyap_upper fliplr(lyap_lower)], ...
     [0.8500 0.3250 0.0980], 'EdgeColor','none','FaceAlpha',0.3,'HandleVisibility','off');

% 均值曲线
plot(x, mean_max_total, 'b-', 'LineWidth', 1.5);
plot(x, mean_lyap_total, 'r-', 'LineWidth', 1.5);

xlabel('Time slot T');
ylabel('Total queue length');
title('Total queue length across different T');

legend('MaxWeight','LyapOpt', ...
       'Location','northwest');

ax = gca;
set(ax,'Box','on');      % 保留外框
set(ax,'Color','w');     % 坐标区白色
set(fig2,'Color','w');   % 整个 figure 白色

export_fig(fig2,'total_b20_vs_large_T.pdf','-pdf','-nocrop');  






