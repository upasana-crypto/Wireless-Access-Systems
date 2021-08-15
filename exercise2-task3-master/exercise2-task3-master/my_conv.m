function c = my_conv(vec1,vec2)

    % my_conv calculates the 1-d convolution of the vectors vec1 and vec2.
    % The end result has the same length as vec1 --> vec1 and vec2 are NOT
    % interchangeable. Attention: The function cuts the END of the
    % convolution to make the length of the result equal to the length of
    % vec1
    % Inputs:
    % 1. vec1: row vector that represents the signal
    % 2. vec2: row vector that represents the channel
    
    % Output:
    % 1. conv 
    
    assert(~isempty(vec1),'the given vec1 is empty');
    assert(~isempty(vec2),'the given vec2 is empty');
    assert( size(vec1,1)==1,'the given vec1 is NOT a row vector');
    assert( size(vec2,1)==1,'the given vec2 is NOT a row vector');
    
    result = zeros(1,length(vec1)+length(vec2)-1);
    for ind_result=1:length(result)
        for ind_vec2=1:length(vec2)
            if ind_result-ind_vec2+1 > 0 && ind_result-ind_vec2+1 <= length(vec1)
                result(ind_result) = result(ind_result) + vec2(ind_vec2)*vec1(ind_result-ind_vec2+1);
            end
        end
    end
   
    % cutting the edge
    c = result(1:end - length(vec2)+1); 

end