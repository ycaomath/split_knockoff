% This script reproduces Table 3 of our paper. The selected regions of
% Split Knockoff can be found in 'connections', with the choice of cross
% validation optimal nu being given in 'optimal_nu'.

%% choose parameters

% add path
root = pwd;
addpath(sprintf('%s/data/AALfeat', root));


data1 = 15;
data2 = 30;

option.q = 0.2;
option.beta = 'cv_beta';
option.frac = 150/752;
option.normalize = true;
option.lambda = 10.^[0: -0.01: -6];
option.nu = 10.^[0: 0.1: 2];

% set random seed
option.rng = 1;

num_nu = length(option.nu);
connections = cell(num_nu, 1);

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

%% set transformation D

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
[results, stat_nu] = split_knockoffs.filter(X, D, y, option);
optimal_nu = zeros(num_nu, 1);

%% generate the table for selected connections
for i = 1: num_nu
    S = results{i}.sk;
    result = Edges(S, :);
    result(:, 3) = [];
    connection = connect.aalLabel.Var2(result);
    connections{i} = connection;
    optimal_nu(i) = stat_nu{i}.nu;
end

% save results
save(sprintf('%s/result/Table_3', pwd));

clearvars -except connections optimal_nu