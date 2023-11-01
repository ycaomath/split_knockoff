% This script reproduces Figure 8 of our paper.
%% choose parameters

k = 20; % the sparsity level
n = 500;% the sample size
p = 100;% the number of features
c = 0.5; % the feature correlation

option = struct;
option.q = 0.2;
option.beta = 'cv_all';
option.frac = 2/5;
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

expo = [-0.5:0.1:0.5];
A_s = 10.^expo;
num_A = length(A_s);
rawvalue = cell(3, num_A);
fdr_knock = zeros(2, option.tests, 2, num_A);
power_knock = zeros(2, option.tests, 2, num_A);


%% start the simulation experiments
for test_A = 1: num_A
    A = A_s(test_A);
    for i = 1: 3
        % choose the respective D for each example
        fprintf('Running simulation experiments for the %d-th SNR and D_%d.\n', test_A, i);
        D = D_s{i};
        simu_data = split_knockoffs.private.simu_unit_cv(n, p, D, A, c, k, option);
        rawvalue{i, test_A} = simu_data.raw;
        if i <= 2
            fdr_knock(i,  :, :, test_A) = simu_data.fdr_k;
            power_knock(i,  :, :, test_A) = simu_data.power_k;
        end
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

fdr_sk = zeros(3, option.tests, 2, num_A);
power_sk = zeros(3, option.tests, 2, num_A);

for test_A = 1: num_A
    for di = 1: 3
        D = D_s{di};
        gamma_true = D * beta_true;
        dataframe = rawvalue{di, test_A};
        for i = 1: option.tests
            stats = dataframe{i};
            Z = stats.Z;
            t_Z = stats.t_Z;
            r = stats.r;
            t_r = stats.t_r;
            Z_tilde = t_Z .* (real(r) == real(t_r));

            % Wst
            W = Z .* sign(Z - Z_tilde);
            S = knockoffs.select(W, option.q, 'knockoff');
            [fdr_sk(di, i, 1, test_A), power_sk(di, i, 1, test_A)] ...
                = split_knockoffs.private.simu_eval(gamma_true, S);
            S = knockoffs.select(W, option.q, 'knockoff+');
            [fdr_sk(di, i, 2, test_A), power_sk(di, i, 2, test_A)] ...
                = split_knockoffs.private.simu_eval(gamma_true, S);

        end
    end
end

mean_fdr_sk = reshape(mean(fdr_sk, 2), [3, 2, num_A]);
mean_power_sk = reshape(mean(power_sk, 2), [3, 2, num_A]);
sd_fdr_sk = reshape(std(fdr_sk, 0,  2), [3, 2, num_A]);
sd_power_sk = reshape(std(power_sk, 0,  2), [3, 2, num_A]);

mean_fdr_knockoff = reshape(mean(fdr_knock, 2), [2, 2, num_A]);
mean_power_knockoff = reshape(mean(power_knock, 2), [2, 2, num_A]);
sd_fdr_knockoff = reshape(std(fdr_knock, 0,  2), [2, 2, num_A]);
sd_power_knockoff = reshape(std(power_knock, 0,  2), [2, 2, num_A]);

save(sprintf('%s/result/Figure_8', pwd));
%% draw plots

t_value = tinv([0.1, 0.9], option.tests - 1);
lower_bound = t_value(1);
upper_bound = t_value(2);

x = expo;
mark = 'abcdef';
for i = 1: 2
    for di = 1: 3
        fdr = reshape(fdr_sk(di, :, i, :), [option.tests, num_A]);
        power = reshape(power_sk(di, :, i, :), [option.tests, num_A]);
        fdr_mid = mean(fdr, 1);
        power_mid = mean(power, 1);
        
        sd_fdr = std(fdr, 0, 1);
        fdr_top = fdr_mid + sd_fdr * upper_bound;
        fdr_bot = fdr_mid + sd_fdr * lower_bound;

        sd_power = std(power, 0, 1);
        power_top = power_mid + sd_power * upper_bound;
        power_bot = power_mid + sd_power * lower_bound;
        
        
        
        x2 = [x, fliplr(x)];
    
        inBetween_fdr = [fdr_top, fliplr(fdr_bot)];
        inBetween_power = [power_top, fliplr(power_bot)];
        if di <= 2
            fdr_mid_knock = reshape(mean_fdr_knockoff(di, i, :), [1, num_A]);
            power_mid_knock = reshape(mean_power_knockoff(di, i, :), [1, num_A]);
            
            sd_fdr_knock = reshape(sd_fdr_knockoff(di, i, :), [1, num_A]);
            fdr_top_knock = fdr_mid_knock + sd_fdr_knock * upper_bound;
            fdr_bot_knock = fdr_mid_knock + sd_fdr_knock * lower_bound;
    
            sd_power_knock = reshape(sd_power_knockoff(di, i, :), [1, num_A]);
            power_top_knock = power_mid_knock + sd_power_knock * upper_bound;
            power_bot_knock = power_mid_knock + sd_power_knock * lower_bound;
            inBetween_fdr_2 = [fdr_top_knock, fliplr(fdr_bot_knock)];
            inBetween_power_2 = [power_top_knock, fliplr(power_bot_knock)];
        end

        fig = figure();
        xlab = x;
        ylab = [0:0.2:1];
        ax = [min(x),max(x),0,1];
        hold on;
        set(fig, 'DefaultTextInterpreter', 'latex');
        grid on
        plot(x, fdr_mid, '-.r');
        plot(x, power_mid, 'r');

        if di <= 2
            plot(x, fdr_mid_knock, '-.b');
            plot(x, power_mid_knock, 'b');
        end
        
        fill(x2, inBetween_fdr, 'r', 'FaceAlpha',0.05, 'LineStyle','none');
        fill(x2, inBetween_power, 'r', 'FaceAlpha',0.05, 'LineStyle','none');
        
        if di <= 2
            fill(x2, inBetween_fdr_2, 'b', 'FaceAlpha',0.05, 'LineStyle','none');
            fill(x2, inBetween_power_2, 'b', 'FaceAlpha',0.05, 'LineStyle','none');
        end
        hold off;

        axis(ax);
        set(gca,'XTick',xlab,'FontSize',20);
        set(gca,'YTick',ylab,'FontSize',20);
        line = refline(0,option.q);
        set(line, 'LineStyle', ':', 'Color', 'black');
        xlabel('$\log(A)$' ,'FontSize',20);
        if di <= 2 && i == 1
            legend('FDR for SK','Power for SK','FDR for Knockoff','Power for Knockoff','FontSize',15);
        elseif di <= 2 && i == 2
            legend('FDR for SK+','Power for SK+','FDR for Knockoff+','Power for Knockoff+','FontSize',15);
        elseif di == 3 && i == 1
            legend('FDR for SK','Power for SK','FontSize',15);
        elseif di == 3 && i == 2
            legend('FDR for SK+','Power for SK+','FontSize',15);
        end
        set(gca,'LooseInset',get(gca,'TightInset'))
        
        saveas(gca,sprintf('plot/figure_8%s', mark((i-1)*3+di)), 'png');
    end
end


%%
clearvars -except mean_fdr* mean_power* sd_fdr* sd_power*
