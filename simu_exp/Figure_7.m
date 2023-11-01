% This script reproduces Figure 7 of our paper.
%% choose parameters

k = 20; % the sparsity level
A = 1; % the signal noise ratio
n = 500;% the sample size
p = 100;% the number of features
c = 0.5; % the feature correlation

option = struct;
option.q = 0.2;
option.beta = 'cv_all';
option.frac = 2/5;
option.tests = 20; % the number of simulation instances
option.normalize = true;
option.lambda = 10.^[0: -0.01: -6];
randomsplits = 100; % the number of random splits

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

rawvalue = cell(3, randomsplits);


%% start the simulation experiments
for split = 1: randomsplits
    option.rng = split;
    for i = 1: 3
        % choose the respective D for each example
        fprintf('Running simulation experiments for the %d-th random split and D_%d.\n', split, i);
        D = D_s{i};
        simu_data = split_knockoffs.private.simu_unit_cv(n, p, D, A, c, k, option);
        rawvalue{i, split} = simu_data.raw;
        save(sprintf('%s/result/temp',pwd));
    end
end

save(sprintf('%s/result/Figure_7', pwd));

%% draw plots for the W = Ws

mark = 'abc';
tests = option.tests;

beta_true = zeros(p, 1);
for i = 1: k
    beta_true(i, 1) = A;
    if rem(i, 3) == 1
        beta_true(i, 1) = 0;
    end
end


for di = 1: 3
    D = D_s{di};
    gamma_true = D * beta_true;
    supp = find(gamma_true ~= 0);    
    possible_set = [];
    for split = 1: randomsplits
        dataframe = rawvalue{di, split};
        for test  = 1: tests
            W = dataframe{test}.Ws;
            result = knockoffs.select(W, option.q, 'knockoff');
            possible_set = union(possible_set, result);
        end
    end

    choices  = length(possible_set);
    freq = zeros(choices, 1);
    for split = 1: randomsplits
        dataframe = rawvalue{di, split};
        for test  = 1: tests
            W = dataframe{test}.Ws;
            result = knockoffs.select(W, option.q, 'knockoff');
            for j = 1: choices
                if ismember(possible_set(j),  result) == true
                    freq(j) = freq(j) + 1/randomsplits/tests;
                end
            end
        end
    end

    [~, I] = sort(freq, 'descend');

    feature_freq = cell(choices, 2);
    for j = 1: choices
        feature_freq{j, 1} = ismember(I(j), supp);
        feature_freq{j, 2} = freq(I(j));
    end


    display = 50;
    nonnull = zeros(50, 1);
    null = zeros(50, 1);
    for j = 1: display
        if feature_freq{j, 1} == 1
            nonnull(j) = feature_freq{j, 2};
        else
            null(j) = feature_freq{j, 2};
        end
    end

    fig = figure();
    hold on
    set(fig, 'DefaultTextInterpreter', 'latex');
    bar(null, 'blue');
    bar(nonnull, 'red');
    ylim([0, 1])
    xticks([ ])
    legend('nulls', 'nonnulls','FontSize', 18)
    hold off
    set(gca,'LooseInset',get(gca,'TightInset'))

    saveas(gca,sprintf('plot/figure_7%s', mark(di)), 'png');
end

%% draw plots for the W = Wst

mark = 'def';


for di = 1: 3
    D = D_s{di};
    gamma_true = D * beta_true;
    supp = find(gamma_true ~= 0);
    possible_set = [];
    for split = 1: randomsplits
        dataframe = rawvalue{di, split};
        for test  = 1: tests
            W = dataframe{test}.Wst;
            result = knockoffs.select(W, option.q, 'knockoff');
            possible_set = union(possible_set, result);
        end
    end

    choices  = length(possible_set);
    freq = zeros(choices, 1);
    for split = 1: randomsplits
        dataframe = rawvalue{di, split};
        for test  = 1: tests
            W = dataframe{test}.Wst;
            result = knockoffs.select(W, option.q, 'knockoff');
            for j = 1: choices
                if ismember(possible_set(j),  result) == true
                    freq(j) = freq(j) + 1/randomsplits/tests;
                end
            end
        end
    end

    [~, I] = sort(freq, 'descend');

    feature_freq = cell(choices, 2);
    for j = 1: choices
        feature_freq{j, 1} = ismember(I(j), supp);
        feature_freq{j, 2} = freq(I(j));
    end



    display = 50;
    nonnull = zeros(50, 1);
    null = zeros(50, 1);
    for j = 1: display
        if feature_freq{j, 1} == 1
            nonnull(j) = feature_freq{j, 2};
        else
            null(j) = feature_freq{j, 2};
        end
    end

    fig = figure();
    hold on
    set(fig, 'DefaultTextInterpreter', 'latex');
    bar(null, 'blue');
    bar(nonnull, 'red');
    ylim([0, 1])
    xticks([ ])
    legend('nulls', 'nonnulls','FontSize', 18)
    hold off
    set(gca,'LooseInset',get(gca,'TightInset'))

    saveas(gca,sprintf('plot/figure_7%s', mark(di)), 'png');
end

%% draw plots for the W = Wbc

mark = 'ghi';


for di = 1: 3
    D = D_s{di};
    gamma_true = D * beta_true;
    supp = find(gamma_true ~= 0);
    possible_set = [];
    for split = 1: randomsplits
        dataframe = rawvalue{di, split};
        for test  = 1: tests
            W = dataframe{test}.Wbc;
            result = knockoffs.select(W, option.q, 'knockoff');
            possible_set = union(possible_set, result);
        end
    end

    choices  = length(possible_set);
    freq = zeros(choices, 1);
    for split = 1: randomsplits
        dataframe = rawvalue{di, split};
        for test  = 1: tests
            W = dataframe{test}.Wbc;
            result = knockoffs.select(W, option.q, 'knockoff');
            for j = 1: choices
                if ismember(possible_set(j),  result) == true
                    freq(j) = freq(j) + 1/randomsplits/tests;
                end
            end
        end
    end

    [~, I] = sort(freq, 'descend');

    feature_freq = cell(choices, 2);
    for j = 1: choices
        feature_freq{j, 1} = ismember(I(j), supp);
        feature_freq{j, 2} = freq(I(j));
    end



    display = 50;
    nonnull = zeros(50, 1);
    null = zeros(50, 1);
    for j = 1: display
        if feature_freq{j, 1} == 1
            nonnull(j) = feature_freq{j, 2};
        else
            null(j) = feature_freq{j, 2};
        end
    end

    fig = figure();
    hold on
    set(fig, 'DefaultTextInterpreter', 'latex');
    bar(null, 'blue');
    bar(nonnull, 'red');
    ylim([0, 1])
    xticks([ ])
    legend('nulls', 'nonnulls','FontSize', 18)
    hold off
    set(gca,'LooseInset',get(gca,'TightInset'))

    saveas(gca,sprintf('plot/figure_7%s', mark(di)), 'png');
end



%%
clearvars
