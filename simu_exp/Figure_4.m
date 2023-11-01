% This script reproduces the Figure 4 of our paper.

%% choose parameters

k = 20; % the sparsity level
A = 1; % the signal noise ratio
n = 750;% the sample size
p = 100;% the number of features
c = 0.5; % the feature correlation

option = struct;
option.q = 0.2;
option.beta = 'cv_beta';
option.frac = 1/5;
option.tests = 200; % the number of simulation instances
option.normalize = true;
option.lambda = 10.^[0: -0.01: -6];

% set random seed
option.rng = 1;

% settings for nu
expo = [0: 0.2: 2];
option.nu = 10.^expo;

num_nu = length(option.nu);

%% generate the linear transformation

D = zeros(5 * p, p);

for i = 1:p
    for  j = 1: 5
        k = mod(i+j-1, p) + 1;
        D(5 * (i-1) + j, i) = 1;
        D(5 * (i-1) + j, k) = -1;
    end
end

beta_true = zeros(p, 1);
for i = 1: k
    beta_true(i, 1) = A;
    if rem(i, 3) == 1
        beta_true(i, 1) = 0;
    end
end



%% start simulation experiments

% choose the respective D for each example
fprintf('Running simulation experiments.\n');
simu_data = split_knockoffs.private.simu_unit(n, p, D, A, c, k, option);
rawvalue = simu_data.raw;
save(sprintf('%s/result/temp',pwd));

%% save results
save(sprintf('%s/result/Figure_4', pwd));

%% start making plots

t_value = tinv([0.1, 0.9], option.tests - 1);
lower_bound = t_value(1);
upper_bound = t_value(2);

ax = [0,2,0,1];
expo = [0:0.2:2];
xlab = [0:0.2:2];
ylab = [0:0.2:1];

%% calculate Ws

fdr_sk = zeros(option.tests, num_nu, 2);
power_sk = zeros(option.tests, num_nu, 2);

gamma_true = D * beta_true;
dataframe = rawvalue;
for i = 1: option.tests
    result_each_time = dataframe{i};
    for j = 1: num_nu
        stats = result_each_time{j};
        Z = stats.Z;
        t_Z = stats.t_Z;
        r = stats.r;
        t_r = stats.t_r;
        Z_tilde = t_Z .* (real(r) == real(t_r));
        W = Z .* sign(Z - t_Z);
        %W = max(Z, t_Z) .* sign(Z - t_Z);
        %W = Z .* sign(Z - Z_tilde);
        %W = max(Z, Z_tilde) .* sign(Z - Z_tilde);
        S = knockoffs.select(W, option.q, 'knockoff');
        [fdr_sk(i, j , 1), power_sk(i, j , 1)] ...
            = split_knockoffs.private.simu_eval(gamma_true, S);
        S = knockoffs.select(W, option.q, 'knockoff+');
        [fdr_sk(i, j , 2), power_sk(i, j , 2)] ...
            = split_knockoffs.private.simu_eval(gamma_true, S);
    end
end

%% make plots for FDR and Power


x = expo;
fdr = fdr_sk;
fdr = reshape(fdr, [option.tests, num_nu, 2]);
fdr_mid = reshape(mean(fdr, 1), [num_nu, 2]);

power = power_sk;
power = reshape(power, [option.tests, num_nu, 2]);
power_mid = reshape(mean(power, 1), [num_nu, 2]);

sd_fdr = reshape(std(fdr, 0, 1), [num_nu, 2]);
fdr_top = fdr_mid + sd_fdr * upper_bound;
fdr_bot = fdr_mid + sd_fdr * lower_bound;

sd_power = reshape(std(power, 0, 1), [num_nu, 2]);
power_top = power_mid + sd_power * upper_bound;
power_bot = power_mid + sd_power * lower_bound;

x2 = [x, fliplr(x)];

inBetween_fdr = [fdr_top(:, 1)', fliplr(fdr_bot(:, 1)')];
inBetween_power = [power_top(:, 1)', fliplr(power_bot(:, 1)')];
inBetween_fdr_2 = [fdr_top(:, 2)', fliplr(fdr_bot(:, 2)')];
inBetween_power_2 = [power_top(:, 2)', fliplr(power_bot(:, 2)')];

fig = figure();
hold on;
set(fig, 'DefaultTextInterpreter', 'latex');
grid on
plot(x, fdr_mid(:, 1), 'r');
plot(x, power_mid(:, 1), '-.r');
plot(x, fdr_mid(:, 2), 'b');
plot(x, power_mid(:, 2), '-.b');
fill(x2, inBetween_fdr, 'r', 'FaceAlpha',0.05, 'LineStyle','none');
fill(x2, inBetween_power, 'r', 'FaceAlpha',0.05, 'LineStyle','none');
fill(x2, inBetween_fdr_2, 'b', 'FaceAlpha',0.05, 'LineStyle','none');
fill(x2, inBetween_power_2, 'b', 'FaceAlpha',0.05, 'LineStyle','none');
hold off;

axis(ax);
set(gca,'XTick',xlab,'FontSize',20);
set(gca,'YTick',xlab,'FontSize',20);
line = refline(0,option.q);
set(line, 'LineStyle', ':', 'Color', 'black');
xlabel('$\log (\nu)$' ,'FontSize',20);
legend('FDR for SK','Power for SK','FDR for SK+','Power for SK+','FontSize',15);
set(gca,'LooseInset',get(gca,'TightInset'))

saveas(gca,'plot/figure_4a.png');


%% calculate Wst

fdr_sk = zeros(option.tests, num_nu, 2);
power_sk = zeros(option.tests, num_nu, 2);

gamma_true = D * beta_true;
dataframe = rawvalue;
for i = 1: option.tests
    result_each_time = dataframe{i};
    for j = 1: num_nu
        stats = result_each_time{j};
        Z = stats.Z;
        t_Z = stats.t_Z;
        r = stats.r;
        t_r = stats.t_r;
        Z_tilde = t_Z .* (real(r) == real(t_r));
        %W = Z .* sign(Z - t_Z);
        %W = max(Z, t_Z) .* sign(Z - t_Z);
        W = Z .* sign(Z - Z_tilde);
        %W = max(Z, Z_tilde) .* sign(Z - Z_tilde);
        S = knockoffs.select(W, option.q, 'knockoff');
        [fdr_sk(i, j , 1), power_sk(i, j , 1)] ...
            = split_knockoffs.private.simu_eval(gamma_true, S);
        S = knockoffs.select(W, option.q, 'knockoff+');
        [fdr_sk(i, j , 2), power_sk(i, j , 2)] ...
            = split_knockoffs.private.simu_eval(gamma_true, S);
    end
end

%% make plots for FDR and Power


x = expo;
fdr = fdr_sk;
fdr = reshape(fdr, [option.tests, num_nu, 2]);
fdr_mid = reshape(mean(fdr, 1), [num_nu, 2]);

power = power_sk;
power = reshape(power, [option.tests, num_nu, 2]);
power_mid = reshape(mean(power, 1), [num_nu, 2]);

sd_fdr = reshape(std(fdr, 0, 1), [num_nu, 2]);
fdr_top = fdr_mid + sd_fdr * upper_bound;
fdr_bot = fdr_mid + sd_fdr * lower_bound;

sd_power = reshape(std(power, 0, 1), [num_nu, 2]);
power_top = power_mid + sd_power * upper_bound;
power_bot = power_mid + sd_power * lower_bound;

x2 = [x, fliplr(x)];

inBetween_fdr = [fdr_top(:, 1)', fliplr(fdr_bot(:, 1)')];
inBetween_power = [power_top(:, 1)', fliplr(power_bot(:, 1)')];
inBetween_fdr_2 = [fdr_top(:, 2)', fliplr(fdr_bot(:, 2)')];
inBetween_power_2 = [power_top(:, 2)', fliplr(power_bot(:, 2)')];

fig = figure();
hold on;
set(fig, 'DefaultTextInterpreter', 'latex');
grid on
plot(x, fdr_mid(:, 1), 'r');
plot(x, power_mid(:, 1), '-.r');
plot(x, fdr_mid(:, 2), 'b');
plot(x, power_mid(:, 2), '-.b');
fill(x2, inBetween_fdr, 'r', 'FaceAlpha',0.05, 'LineStyle','none');
fill(x2, inBetween_power, 'r', 'FaceAlpha',0.05, 'LineStyle','none');
fill(x2, inBetween_fdr_2, 'b', 'FaceAlpha',0.05, 'LineStyle','none');
fill(x2, inBetween_power_2, 'b', 'FaceAlpha',0.05, 'LineStyle','none');
hold off;

axis(ax);
set(gca,'XTick',xlab,'FontSize',20);
set(gca,'YTick',xlab,'FontSize',20);
line = refline(0,option.q);
set(line, 'LineStyle', ':', 'Color', 'black');
xlabel('$\log (\nu)$' ,'FontSize',20);
legend('FDR for SK','Power for SK','FDR for SK+','Power for SK+','FontSize',15);
set(gca,'LooseInset',get(gca,'TightInset'))

saveas(gca,'plot/figure_4b.png');


%% calculate Wbc

fdr_sk = zeros(option.tests, num_nu, 2);
power_sk = zeros(option.tests, num_nu, 2);

gamma_true = D * beta_true;
dataframe = rawvalue;
for i = 1: option.tests
    result_each_time = dataframe{i};
    for j = 1: num_nu
        stats = result_each_time{j};
        Z = stats.Z;
        t_Z = stats.t_Z;
        r = stats.r;
        t_r = stats.t_r;
        Z_tilde = t_Z .* (real(r) == real(t_r));
        %W = Z .* sign(Z - t_Z);
        W = max(Z, t_Z) .* sign(Z - t_Z);
        %W = Z .* sign(Z - Z_tilde);
        %W = max(Z, Z_tilde) .* sign(Z - Z_tilde);
        S = knockoffs.select(W, option.q, 'knockoff');
        [fdr_sk(i, j , 1), power_sk(i, j , 1)] ...
            = split_knockoffs.private.simu_eval(gamma_true, S);
        S = knockoffs.select(W, option.q, 'knockoff+');
        [fdr_sk(i, j , 2), power_sk(i, j , 2)] ...
            = split_knockoffs.private.simu_eval(gamma_true, S);
    end
end


%% make plots for FDR and Power



x = expo;
fdr = fdr_sk;
fdr = reshape(fdr, [option.tests, num_nu, 2]);
fdr_mid = reshape(mean(fdr, 1), [num_nu, 2]);

power = power_sk;
power = reshape(power, [option.tests, num_nu, 2]);
power_mid = reshape(mean(power, 1), [num_nu, 2]);

sd_fdr = reshape(std(fdr, 0, 1), [num_nu, 2]);
fdr_top = fdr_mid + sd_fdr * upper_bound;
fdr_bot = fdr_mid + sd_fdr * lower_bound;

sd_power = reshape(std(power, 0, 1), [num_nu, 2]);
power_top = power_mid + sd_power * upper_bound;
power_bot = power_mid + sd_power * lower_bound;

x2 = [x, fliplr(x)];

inBetween_fdr = [fdr_top(:, 1)', fliplr(fdr_bot(:, 1)')];
inBetween_power = [power_top(:, 1)', fliplr(power_bot(:, 1)')];
inBetween_fdr_2 = [fdr_top(:, 2)', fliplr(fdr_bot(:, 2)')];
inBetween_power_2 = [power_top(:, 2)', fliplr(power_bot(:, 2)')];

fig = figure();
hold on;
set(fig, 'DefaultTextInterpreter', 'latex');
grid on
plot(x, fdr_mid(:, 1), 'r');
plot(x, power_mid(:, 1), '-.r');
plot(x, fdr_mid(:, 2), 'b');
plot(x, power_mid(:, 2), '-.b');
fill(x2, inBetween_fdr, 'r', 'FaceAlpha',0.05, 'LineStyle','none');
fill(x2, inBetween_power, 'r', 'FaceAlpha',0.05, 'LineStyle','none');
fill(x2, inBetween_fdr_2, 'b', 'FaceAlpha',0.05, 'LineStyle','none');
fill(x2, inBetween_power_2, 'b', 'FaceAlpha',0.05, 'LineStyle','none');
hold off;

axis(ax);
set(gca,'XTick',xlab,'FontSize',20);
set(gca,'YTick',xlab,'FontSize',20);
line = refline(0,option.q);
set(line, 'LineStyle', ':', 'Color', 'black');
xlabel('$\log (\nu)$' ,'FontSize',20);
legend('FDR for SK','Power for SK','FDR for SK+','Power for SK+','FontSize',15);
set(gca,'LooseInset',get(gca,'TightInset'))

saveas(gca,'plot/figure_4c.png');


clearvars