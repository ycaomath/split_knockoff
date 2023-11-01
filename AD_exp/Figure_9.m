% This script reproduces Figure 9 of our paper. The frequencies of selected
% connections by Split Knockoffs under cross validation optimal nu can be found
% in 'connection_freq' in a descending order. The whole calculation
% may takes hours to finish. 

%% choose parameters

% add path
root = pwd;
addpath(sprintf('%s/data/AALfeat', root));


data1 = 15;
data2 = 30;

option.q = 0.2;
option.beta = 'cv_all';
option.frac = 150/752;
option.normalize = true;
option.lambda = 10.^[0: -0.01: -6];
option.nu = 10.^[0: 0.4: 2];

num_nu = length(option.nu);
result_edge = cell(num_nu, 1);

%% generate the response vector y
% load 15 T data
y_AD = load(sprintf('ADAS_%d_AD.mat', data1)); 
y_AD = y_AD.ADAS;
y_MCI = load(sprintf('ADAS_%d_MCI.mat', data1));
y_MCI = y_MCI.ADAS;
y_NC = load(sprintf('ADAS_%d_NC.mat', data1));
y_NC = y_NC.ADAS;
y1 = [y_AD;y_MCI;y_NC];

% load 30 T data
y_AD = load(sprintf('ADAS_%d_AD.mat', data2)); 
y_AD = y_AD.ADAS;
y_MCI = load(sprintf('ADAS_%d_MCI.mat', data2));
y_MCI = y_MCI.ADAS;
y_NC = load(sprintf('ADAS_%d_NC.mat', data2));
y_NC = y_NC.ADAS;
y2 = [y_AD;y_MCI;y_NC];

y = [y1; y2];
index = find( y < 0);
y(index) = []; % rule out samples with invalid scores
y = split_knockoffs.private.normc(y);


%% generate the design matrix X
% load 15 T data
X_AD = load(sprintf('AAL_%d_AD_feature_TIV.mat', data1));
X_AD = X_AD.feature_TIV;
X_MCI = load(sprintf('AAL_%d_MCI_feature_TIV.mat', data1));
X_MCI = X_MCI.feature_TIV;
X_NC = load(sprintf('AAL_%d_NC_feature_TIV.mat', data1));
X_NC = X_NC.feature_TIV;
X_1 = [X_AD; X_MCI;X_NC];

% load 30 T data
X_AD = load(sprintf('AAL_%d_AD_feature_TIV.mat', data2));
X_AD = X_AD.feature_TIV;
X_MCI = load(sprintf('AAL_%d_MCI_feature_TIV.mat', data2));
X_MCI = X_MCI.feature_TIV;
X_NC = load(sprintf('AAL_%d_NC_feature_TIV.mat', data2));
X_NC = X_NC.feature_TIV;
X_2 = [X_AD; X_MCI;X_NC];
X = [X_1; X_2];
X = double(X);
[n,~] = size(X);
X(index, :) = []; % rule out samples with invalid scores
id = [1:1:90]; % select the Cerebrum
X = X(:,id);
X = split_knockoffs.private.normc(X);

%% set the transformation D
p = length(id);
connect = load('aalConect.mat');
Edges = connect.B.Edges;
Edges = Edges.Variables;
Edges = Edges(:,:);
Edges = Edges(ismember(Edges(:,1),id) & ismember(Edges(:,2),id),:);

m = size(Edges,1);
D = zeros(m,p);
for i = 1:m
    D(i,Edges(i,1)) = 1;
    D(i,Edges(i,2)) = -1;
end
D = D(:,id); % select the Cerebrum

%% conduct Split Knockoff

randoms  = 100;
results = cell(randoms, 1);
stat_nu = cell(randoms, 1);
for i = 1: randoms
    option.rng = i;
    [results{i}, stat_nu{i}] = split_knockoffs.filter(X, D, y, option);
    fprintf('%.0f%%...',i/randoms*100);
    if mod(i, 20) == 0 && i ~= randoms
        fprintf('\n');
    end
    if i == randoms
        fprintf('done.\n');
    end
end



%% calculate selection frequencies

possible_set = [];
for i = 1: randoms
    possible_set = union(possible_set, results{i}.sk);
end

choices  = length(possible_set);
freq = zeros(choices, 1);
for i = 1: randoms
    for j = 1: choices
        if ismember(possible_set(j),  results{i}.sk) == true
            freq(j) = freq(j) + 1/randoms;
        end
    end
end
[~, I] = sort(freq, 'descend');
connection_freq = cell(choices, 3);
for j = 1: choices
    vertices  = Edges(possible_set(I(j)), :);
    connection_freq{j, 1} = connect.aalLabel.Var2(vertices(1));
    connection_freq{j, 2} = connect.aalLabel.Var2(vertices(2));
    connection_freq{j, 3} = freq(I(j));
end
%% make a plot
freq_sorted = sort(freq, 'descend');
plot_freq = freq_sorted(1:10);
plot_freq = sort(plot_freq);
barh(plot_freq, 0.3)

% labels added according to 'connection_freq'
yticklabels({ 'HIP(L) & INS(L)',  'HIP(R) & PUT(R)','HIP(R) & PreCG(R)', 'MTG(L) & ROL(L)', 'HIP(R) & TPOmid(R)','HIP(L) & PreCG(L)', 'HIP(R) & INS(R)', 'HIP(L) & PUT(L)', 'HIP(R) & LING(R)', 'HIP(L) & LING(L)'});
xlabel('Frequency')

set(gca,'LooseInset',get(gca,'TightInset'))
saveas(gca, 'plot/Figure_9', 'png');

% save results
save(sprintf('%s/result/Figure_9', pwd));

clearvars -except connection_freq