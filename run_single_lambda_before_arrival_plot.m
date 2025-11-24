


% clear; clc;

%% Parameters
n = 8;
m = n*10;
T = 100;
runs = 10000;
scenario = 'b_epsilon';
sigma_A = 1 * ones(n,1);
sigma_D = 1 * ones(n,1);
epsilon = 0.1;
scalar = 10;
alpha = 0.05;  % 95% confidence
%% Generate D and facets
 [D,D_top, facets, goodIdx] = genDepartureSet(m, n, scalar, 42);
 

%% Generate one lambda
 lambda = genLambda(D, D_top, facets, goodIdx, epsilon, [], scenario);

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
        % A_curr
        % d_lyap

        Q_max = max(Q_max - d_max, 0) + A_curr;
        Q_lyap = max(Q_lyap - d_lyap, 0) + A_curr;

        Metric_max_ba(run, t)  = sum(Q_max-A_curr); %before arrivals
        Metric_lyap_ba(run, t) = sum(Q_lyap-A_curr);
        Metric_max_total(run, t)  = sum(Q_max);
        Metric_lyap_total(run, t) = sum(Q_lyap);
    end
end

%% Post-processing
df = runs - 1;
tcrit = tinv(1 - alpha/2, df);
x = 1:T;

mean_max_ba = mean(Metric_max_ba,1);
mean_lyap_ba = mean(Metric_lyap_ba,1);
se_max_ba = std(Metric_max_ba,0,1)/sqrt(runs);
se_lyap_ba = std(Metric_lyap_ba,0,1)/sqrt(runs);
hw_max_ba = tcrit * se_max_ba;
hw_lyap_ba = tcrit * se_lyap_ba;

mean_max_total = mean(Metric_max_total,1);
mean_lyap_total = mean(Metric_lyap_total,1);
se_max_total = std(Metric_max_total,0,1)/sqrt(runs);
se_lyap_total = std(Metric_lyap_total,0,1)/sqrt(runs);
hw_max_total = tcrit * se_max_total;
hw_lyap_total = tcrit * se_lyap_total;




%% Plot total queue length before arrivals


fig1 = figure('Units','inches','Position',[1,1,6,4]);
set(fig1,'Color','w','InvertHardcopy','off');
axis tight
set(gca, 'Position', [0.13 0.15 0.75 0.75])
hold on;

% Draw shaded confidence intervals without legend entries
h1 = fill([x fliplr(x)], [mean_max_ba+hw_max_ba fliplr(mean_max_ba-hw_max_ba)], [0 0.4470 0.7410], 'EdgeColor','none', 'FaceAlpha',0.3);
h1.HandleVisibility = 'off';
h2 = fill([x fliplr(x)], [mean_lyap_ba+hw_lyap_ba fliplr(mean_lyap_ba-hw_lyap_ba)], [0.8500 0.3250 0.0980], 'EdgeColor','none', 'FaceAlpha',0.3);
h2.HandleVisibility = 'off';

% Plot performance lines with DisplayName for legend
plot(x, mean_max_ba, 'b-', 'LineWidth', 1.5, 'DisplayName', 'MaxWeight');
plot(x, mean_lyap_ba, 'r-', 'LineWidth', 1.5, 'DisplayName', 'LyapOpt');

xlabel('Time slot T');
ylabel('Total queue lengths before arrivals');
title('Total queue lengths before arrivals performance');
legend('show', 'Location', 'southeast');

% Show all box borders
box on;

set(findall(fig1, '-property', 'FontSize'), 'FontSize', 10);
export_fig(fig1, 'sq_n8.pdf', '-pdf','-nocrop');

%% Plot Total Queue Length
fig2 = figure('Units','inches','Position',[1,1,6,4]);
set(fig2,'Color','w','InvertHardcopy','off');
axis tight
set(gca, 'Position', [0.13 0.15 0.75 0.75])
hold on;

% Draw shaded confidence intervals without legend entries
h3 = fill([x fliplr(x)], [mean_max_total+hw_max_total fliplr(mean_max_total-hw_max_total)], [0 0.4470 0.7410], 'EdgeColor','none', 'FaceAlpha',0.3);
h3.HandleVisibility = 'off';
h4 = fill([x fliplr(x)], [mean_lyap_total+hw_lyap_total fliplr(mean_lyap_total-hw_lyap_total)], [0.8500 0.3250 0.0980], 'EdgeColor','none', 'FaceAlpha',0.3);
h4.HandleVisibility = 'off';

% Plot performance lines with DisplayName for legend
plot(x, mean_max_total, 'b-', 'LineWidth', 1.5, 'DisplayName', 'MaxWeight');
plot(x, mean_lyap_total, 'r-', 'LineWidth', 1.5, 'DisplayName', 'LyapOpt');

xlabel('Time slot T');
ylabel('Total queue length');
title('Total queue length performance');
legend('show', 'Location', 'southeast');

% Show all box borders
box on;

set(findall(fig2, '-property', 'FontSize'), 'FontSize', 10);
export_fig(fig2, 'total_n8.pdf', '-pdf','-nocrop');








