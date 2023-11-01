function[fdr, power] = simu_eval(gamma_true, result)
% this function calculates the FDR and Power in simulation experiments.
%
% input arguments
% gamma_true: the true value of gamma
% result: the estimated support set of gamma
%
% output arguments
% fdr: false discovery rate of the estimated support set
% power: power of the estimated support set

fdr = sum(gamma_true(result) == 0) / max(length(result), 1);

power = sum(gamma_true(result) ~= 0) / sum(gamma_true ~= 0);

end