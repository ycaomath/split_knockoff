%% add the path of this package
addpath(pwd)

%% check dependency
if ~exist('glmnet_matlab','dir')
    disp('Package glmnet_matlab is required!')
    addpath('../glmnet_matlab/')
end