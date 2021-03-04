function[POWER, FDR] = Power_FDR_for_simulation(gamma_true, results)
    % calculte Power and FDR for simulations where gamma_true is known
    % preprocessing, transfer results to m*1 matrix
    m = length(gamma_true);
    result = zeros(m, 1);
    for j = 1: length(results)
        result(results(1, j), 1) = 1;
    end
    
    % calculate FDR
    FDR = 0;
    total = sum(result);

    for j = 1:m 
        if result(j, 1) == 1
            if gamma_true(j, 1) == 0
                FDR = FDR + 1;
            end
        end
    end

    if total == 0
        FDR = 0;
    else
        FDR = FDR / total;
    end
    
    % calculate Power
    POWER = 0;
    total_p = sum(sum(gamma_true ~= 0));

    for j = 1:m
        if result(j, 1) == 1
            if gamma_true(j, 1) ~= 0
                POWER = POWER + 1;
            end
        end
    end

    POWER = POWER / total_p;
end