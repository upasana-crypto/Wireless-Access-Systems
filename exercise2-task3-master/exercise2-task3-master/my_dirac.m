function d = my_dirac(number)

    % the my_dirac calculates the dirac function, but instead of Inf has 1 at number=0.
    % Inputs:
    % 1. number:
    % Output:
    % 1. value of my_dirac
    
    if number == 0
        d = 1;
    else
        d = 0;
    end

end