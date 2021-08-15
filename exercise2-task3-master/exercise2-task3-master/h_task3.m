function h = h_task3(h_coeff,lamda,n)

    % the my_dirac calculates the dirac function, but instead of Inf has 1 at number=0.
    % Inputs:
    % 1. coefficients
    % 2. labda: shows decay
    % 3. time sample
    
    % Output:
    % 1. value of h_function
    
    h = 0;
    for ind_h=1:length(h_coeff)
        h = h + h_coeff(ind_h)*exp(-lamda*(ind_h-1))*my_dirac(n-ind_h);
    end
end