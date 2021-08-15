function p = power_calc(complex_matrix)

    % the power_calculation calculates the power of the complex vector  
    % Inputs:
    % 1. complex vector that is a row vector
    % Output:
    % 1. power defined as : sum()/length
    
    assert(~isempty(complex_matrix),'the given matrix is empty');
    if size(complex_matrix,1)==1
        p = sum(abs(complex_matrix).^2)/length(complex_matrix);
    else
        temp_p = zeros(1,size(complex_matrix,1));
        for i=1:size(complex_matrix,1)
            for j = 1:size(complex_matrix,2)
                temp_p(i) = temp_p(i)+abs(complex_matrix(i,j))^2;
            end
        end
        p = temp_p/size(complex_matrix,2);
    end
end