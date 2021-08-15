function [ber,errors] = calc_ber_err(x_Rx,x_Tx)

    % calculate_ber compares 2 vectors: x_Rx and x_Tx and computes the
    % errros made. It also computes the ber by diving errors with the
    % number of elements of the vectors

    % Inputs:
    % 1. x_Rx: row vector that represents the received bits
    % 2. x_Tx: row vector that represents the sent bits
    
    % Output:
    % 1. ber
    % 2. errors: number of errors found
    
    assert(~isempty(x_Rx),'the given x_Rx is empty');
    assert(~isempty(x_Tx),'the given x_Tx is empty');
    assert( size(x_Rx,1)==1,'the given x_Rx is NOT a row vector');
    assert( size(x_Tx,1)==1,'the given x_Tx is NOT a row vector');
    assert( size(x_Rx,2) == size(x_Tx,2),'the given vectors do NOT have the same length ');
    
    result = 0;
    for i=1:length(x_Rx)
        if x_Rx(i) ~= x_Tx(i)
            result = result+1;
        end
    end
    
    errors = result;
    ber = result / size(x_Rx,2);

end