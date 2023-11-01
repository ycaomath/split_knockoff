% This script reproduces Figure 6 of our paper.
%% choose parameters

k = 20; % the sparsity level
A = 0.5; % the signal noise ratio
n = 500;% the sample size
p = 100;% the number of features
c = 0.5; % the feature correlation

option = struct;
option.q = 0.2;
option.beta = 'cv_all';
option.tests = 200; % the number of simulation instances
option.normalize = true;
option.lambda = 10.^[0: -0.01: -6];

% set random seed
option.rng = 1;

%% generate 3 types of transformation

% generate D1, D2, and D3
D_G = zeros(p-1, p);

for i = 1:(p-1)
    D_G(i, i) = 1;
    D_G(i, i+1) = -1;
end

D_1 = eye(p);
D_2 = D_G;
D_3 = [eye(p); D_G];
D_s = {D_1, D_2, D_3};

test_frac = 10;
rawvalue = cell(3, test_frac - 2);


%% start the simulation experiments
for frac = 1: test_frac - 2
    option.frac = frac / test_frac;
    for i = 1: 3
        % choose the respective D for each example
        fprintf('Running simulation experiments for the %d-th split fraction and D_%d.\n', frac, i);
        D = D_s{i};
        simu_data = split_knockoffs.private.simu_unit_cv(n, p, D, A, c, k, option);
        rawvalue{i, frac} = simu_data.raw;
        save(sprintf('%s/result/temp',pwd));
    end
end


%% start calculating results on 3 types of W statistics

beta_true = zeros(p, 1);
for i = 1: k
    beta_true(i, 1) = A;
    if rem(i, 3) == 1
        beta_true(i, 1) = 0;
    end
end

fdr_sk = zeros(3, option.tests, 3, 2, test_frac - 2);
power_sk = zeros(3, option.tests, 3, 2, test_frac - 2);

for frac = 1: test_frac - 2
    for di = 1: 3
        D = D_s{di};
        gamma_true = D * beta_true;
        dataframe = rawvalue{di, frac};
        for i = 1: option.tests
            stats = dataframe{i};
            Z = stats.Z;
            t_Z = stats.t_Z;
            r = stats.r;
            t_r = stats.t_r;
            Z_tilde = t_Z .* (real(r) == real(t_r));

            % Ws
            W = Z .* sign(Z - t_Z);
            S = knockoffs.select(W, option.q, 'knockoff');
            [fdr_sk(di, i, 1, 1, frac), power_sk(di, i, 1, 1, frac)] ...
                = split_knockoffs.private.simu_eval(gamma_true, S);
            S = knockoffs.select(W, option.q, 'knockoff+');
            [fdr_sk(di, i, 1, 2, frac), power_sk(di, i, 1, 2, frac)] ...
                = split_knockoffs.private.simu_eval(gamma_true, S);

            % Wst
            W = Z .* sign(Z - Z_tilde);
            S = knockoffs.select(W, option.q, 'knockoff');
            [fdr_sk(di, i, 2, 1, frac), power_sk(di, i, 2, 1, frac)] ...
                = split_knockoffs.private.simu_eval(gamma_true, S);
            S = knockoffs.select(W, option.q, 'knockoff+');
            [fdr_sk(di, i, 2, 2, frac), power_sk(di, i, 2, 2, frac)] ...
                = split_knockoffs.private.simu_eval(gamma_true, S);

            % Wbc
            W = max(Z, t_Z) .* sign(Z - t_Z);
            S = knockoffs.select(W, option.q, 'knockoff');
            [fdr_sk(di, i, 3, 1, frac), power_sk(di, i, 3, 1, frac)] ...
                = split_knockoffs.private.simu_eval(gamma_true, S);
            S = knockoffs.select(W, option.q, 'knockoff+');
            [fdr_sk(di, i, 3, 2, frac), power_sk(di, i, 3, 2, frac)] ...
                = split_knockoffs.private.simu_eval(gamma_true, S);
        end
    end
end

mean_fdr_sk = reshape(mean(fdr_sk, 2), [3, 3, 2, test_frac - 2]);
mean_power_sk = reshape(mean(power_sk, 2), [3, 3, 2, test_frac - 2]);
sd_fdr_sk = reshape(std(fdr_sk, 0,  2), [3, 3, 2, test_frac - 2]);
sd_power_sk = reshape(std(power_sk, 0,  2), [3, 3, 2, test_frac - 2]);

save(sprintf('%s/result/Figure_6', pwd));
%% draw plots

t_value = tinv([0.1, 0.9], option.tests - 1);
lower_bound = t_value(1);
upper_bound = t_value(2);

x = [0.1: 0.1: 0.8];
mark = 'abcdefghi';
for i = 1: 3
    for di = 1: 3
        fdr = reshape(fdr_sk(di, :, i, :, :), [option.tests, 2, test_frac - 2]);
        power = reshape(power_sk(di, :, i, :, :), [option.tests, 2, test_frac - 2]);
        fdr_mid = reshape(mean(fdr, 1), [2, test_frac - 2]);
        power_mid = reshape(mean(power, 1), [2, test_frac - 2]);
        
        sd_fdr = reshape(std(fdr, 0, 1), [2, test_frac - 2]);
        fdr_top = fdr_mid + sd_fdr * upper_bound;
        fdr_bot = fdr_mid + sd_fdr * lower_bound;

        sd_power = reshape(std(power, 0, 1), [2, test_frac - 2]);
        power_top = power_mid + sd_power * upper_bound;
        power_bot = power_mid + sd_power * lower_bound;
        
        x2 = [x, fliplr(x)];
    
        inBetween_fdr = [fdr_top(1, :), fliplr(fdr_bot(1, :))];
        inBetween_power = [power_top(1, :), fliplr(power_bot(1, :))];
        inBetween_fdr_2 = [fdr_top(2, :), fliplr(fdr_bot(2, :))];
        inBetween_power_2 = [power_top(2, :), fliplr(power_bot(2, :))];

        fig = figure();
        xlab = [0.1:0.1:0.8];
        ylab = [0:0.2:1];
        ax = [0.1,0.8,0,1];
        hold on;
        set(fig, 'DefaultTextInterpreter', 'latex');
        grid on
        plot(x, fdr_mid(1, :), 'r');
        plot(x, power_mid(1, :), '-.r');
        plot(x, fdr_mid(2, :), 'b');
        plot(x, power_mid(2, :), '-.b');
        fill(x2, inBetween_fdr, 'r', 'FaceAlpha',0.05, 'LineStyle','none');
        fill(x2, inBetween_power, 'r', 'FaceAlpha',0.05, 'LineStyle','none');
        fill(x2, inBetween_fdr_2, 'b', 'FaceAlpha',0.05, 'LineStyle','none');
        fill(x2, inBetween_power_2, 'b', 'FaceAlpha',0.05, 'LineStyle','none');
        hold off;

        axis(ax);
        set(gca,'XTick',xlab,'FontSize',20);
        set(gca,'YTick',ylab,'FontSize',20);
        line = refline(0,option.q);
        set(line, 'LineStyle', ':', 'Color', 'black');
        xlabel('$\frac{\# \mathcal{D}_1}{\# \mathcal{D}}$' ,'FontSize',20);
        legend('FDR for SK','Power for SK','FDR for SK+','Power for SK+','FontSize',15);
        set(gca,'LooseInset',get(gca,'TightInset'))
        
        saveas(gca,sprintf('plot/figure_6%s', mark((i-1)*3+di)), 'png');
    end
end


%%
clearvars -except mean_fdr* mean_power* sd_fdr* sd_power*
