function x_hat=demapping(s, constellationType)

    % the demapping function maps the symbol sequence s with a given
    % constellation type to a bit sequence x_hat, which contains elements
    % from the {0,1}
    % Inputs:
    % 1. s: symbol sequence 
    % 2. constellation type: BPSK, QPSK, 16-QAM, 64-QAM
    % Output:
    % 1. x: bit sequence
    
    
    %% constellation type input control
    % M is the number of points in the signal constellation. 
    
    % It is assumed that the receiver knows the average symbol energy 
    E=1;
    
    if strcmp(constellationType, 'BPSK')
        %% BPSK
        
        M=2;
        k=log2(M);        

        % hardcoded scaling factor
        scaling_factor_BPSK = sqrt(E);
        
        % scaled demapper_BPSK
        Scaled_Demapper_BPSK{1} = scaling_factor_BPSK * (-1);
        Scaled_Demapper_BPSK{2} = scaling_factor_BPSK * 1;  
        
        % initialize x_hat vector
        x_hat = zeros(1, length(s));
        
        % Apply ML decoder
        for temp=1:k:length(s)
            dist_1 = abs((s(temp)-Scaled_Demapper_BPSK{1}));
            dist_2 = abs((s(temp)-Scaled_Demapper_BPSK{2}));
            if dist_1 <= dist_2
                x_hat(temp)=0;
            else
                x_hat(temp)=1;
            end            
        end        
        
        
    elseif strcmp(constellationType, 'QPSK')
        %% QPSK
        M=4;
        k=log2(M);
        
        % hardcoded scaling factor
        scaling_factor_QPSK = sqrt(E/2);
        
        % Definition of Demapper
        Scaled_Demapper_QPSK{1} = scaling_factor_QPSK * [-1;-1];
        Scaled_Demapper_QPSK{2} = scaling_factor_QPSK * [1;-1];
        Scaled_Demapper_QPSK{3} = scaling_factor_QPSK * [-1;1];
        Scaled_Demapper_QPSK{4} = scaling_factor_QPSK * [1;1];
        
        % initialize x_hat
        x_hat=zeros(1, 2 * (size (s,2)));
        
        % initialize dist vector, to calculate the different distances
        dist=zeros(1,M);
        
        % Apply ML decoder
        bit_index=0;
        for temp=1:size(s,2)
            for temp_j=1:M
                dist(temp_j) = sqrt( (s(1,temp) - Scaled_Demapper_QPSK{temp_j}(1) )^2 + (s(2,temp) - Scaled_Demapper_QPSK{temp_j}(2) )^2  );
            end
            [~,min_index] = min(dist);
            demapped_bits_from_current_symbol=de2bi(min_index-1,k, 'left-msb');
            for i=1:k
                bit_index=bit_index+1;
                x_hat(bit_index)=demapped_bits_from_current_symbol(i);
            end                
        end   

        
    elseif strcmp(constellationType, '16-QAM')
        %% 16-QAM
        M=16;
        k=log2(M);
        
        % hardcoded scaling factor
        scaling_factor_16_QAM = sqrt(E/10);
        
        % Definition of Demapper
        Scaled_Demapper_16_QAM{1} = scaling_factor_16_QAM * [-3;-3];
        Scaled_Demapper_16_QAM{2} = scaling_factor_16_QAM * [-1;-3];
        Scaled_Demapper_16_QAM{3} = scaling_factor_16_QAM * [-3;-1];
        Scaled_Demapper_16_QAM{4} = scaling_factor_16_QAM * [-1;-1];
        Scaled_Demapper_16_QAM{5} = scaling_factor_16_QAM * [3;-3];
        Scaled_Demapper_16_QAM{6} = scaling_factor_16_QAM * [1;-3];
        Scaled_Demapper_16_QAM{7} = scaling_factor_16_QAM * [3;-1];
        Scaled_Demapper_16_QAM{8} = scaling_factor_16_QAM * [1;-1];
        Scaled_Demapper_16_QAM{9} = scaling_factor_16_QAM * [-3;3];
        Scaled_Demapper_16_QAM{10} = scaling_factor_16_QAM * [-1;3];
        Scaled_Demapper_16_QAM{11} = scaling_factor_16_QAM * [-3;1];
        Scaled_Demapper_16_QAM{12} = scaling_factor_16_QAM * [-1;1];
        Scaled_Demapper_16_QAM{13} = scaling_factor_16_QAM * [3;3];
        Scaled_Demapper_16_QAM{14} = scaling_factor_16_QAM * [1;3];
        Scaled_Demapper_16_QAM{15} = scaling_factor_16_QAM * [3;1];
        Scaled_Demapper_16_QAM{16} = scaling_factor_16_QAM * [1;1]; 
        
        % initialize x_hat
        x_hat=zeros(1, 2 * (size (s,2)));
        
        % initialize dist vector, to calculate the different distances
        dist=zeros(1,M);
        
        % Apply ML decoder
        bit_index=0;
        for temp=1:size(s,2)
            for temp_j=1:M
                dist(temp_j) = sqrt( (s(1,temp) - Scaled_Demapper_16_QAM{temp_j}(1) )^2 + (s(2,temp) - Scaled_Demapper_16_QAM{temp_j}(2) )^2  );
            end
            [~,min_index] = min(dist);
            demapped_bits_from_current_symbol=de2bi(min_index-1,k, 'left-msb');
            for i=1:k
                bit_index=bit_index+1;
                x_hat(bit_index)=demapped_bits_from_current_symbol(i);
            end                
        end        
        
        
        
    elseif strcmp(constellationType, '64-QAM')
        %% 64-QAM
        
        M=64;
        
        % k contains the number of mapped bits per symbol
        k = log2(M);
        
        % hardcode scaling factor
        scaling_factor_64_QAM = sqrt(E/42);
        
        % Definition of Demapper
        Scaled_Demapper_64_QAM{1} = scaling_factor_64_QAM * [-7;-7];
        Scaled_Demapper_64_QAM{2} = scaling_factor_64_QAM * [-5;-7];
        Scaled_Demapper_64_QAM{3} = scaling_factor_64_QAM * [-7;-5];
        Scaled_Demapper_64_QAM{4} = scaling_factor_64_QAM * [-5;-5];
        Scaled_Demapper_64_QAM{5} = scaling_factor_64_QAM * [-1;-7];
        Scaled_Demapper_64_QAM{6} = scaling_factor_64_QAM * [-3;-7];
        Scaled_Demapper_64_QAM{7} = scaling_factor_64_QAM * [-1;-5];
        Scaled_Demapper_64_QAM{8} = scaling_factor_64_QAM * [-3;-5];
        Scaled_Demapper_64_QAM{9} = scaling_factor_64_QAM * [-7;-1];
        Scaled_Demapper_64_QAM{10} = scaling_factor_64_QAM * [-5;-1];
        Scaled_Demapper_64_QAM{11} = scaling_factor_64_QAM * [-7;-3];
        Scaled_Demapper_64_QAM{12} = scaling_factor_64_QAM * [-5;-3];
        Scaled_Demapper_64_QAM{13} = scaling_factor_64_QAM * [-1;-1];
        Scaled_Demapper_64_QAM{14} = scaling_factor_64_QAM * [-3;-1];
        Scaled_Demapper_64_QAM{15} = scaling_factor_64_QAM * [-1;-3];
        Scaled_Demapper_64_QAM{16} = scaling_factor_64_QAM * [-3;-3];        
        Scaled_Demapper_64_QAM{17} = scaling_factor_64_QAM * [7;-7];
        Scaled_Demapper_64_QAM{18} = scaling_factor_64_QAM * [5;-7];
        Scaled_Demapper_64_QAM{19} = scaling_factor_64_QAM * [7;-5];
        Scaled_Demapper_64_QAM{20} = scaling_factor_64_QAM * [5;-5];        
        Scaled_Demapper_64_QAM{21} = scaling_factor_64_QAM * [1;-7];
        Scaled_Demapper_64_QAM{22} = scaling_factor_64_QAM * [3;-7];
        Scaled_Demapper_64_QAM{23} = scaling_factor_64_QAM * [1;-5];
        Scaled_Demapper_64_QAM{24} = scaling_factor_64_QAM * [3;-5];        
        Scaled_Demapper_64_QAM{25} = scaling_factor_64_QAM * [7;-1];
        Scaled_Demapper_64_QAM{26} = scaling_factor_64_QAM * [5;-1];
        Scaled_Demapper_64_QAM{27} = scaling_factor_64_QAM * [7;-3];
        Scaled_Demapper_64_QAM{28} = scaling_factor_64_QAM * [5;-3];        
        Scaled_Demapper_64_QAM{29} = scaling_factor_64_QAM * [1;-1];
        Scaled_Demapper_64_QAM{30} = scaling_factor_64_QAM * [3;-1];
        Scaled_Demapper_64_QAM{31} = scaling_factor_64_QAM * [1;-3];
        Scaled_Demapper_64_QAM{32} = scaling_factor_64_QAM * [3;-3];        
        Scaled_Demapper_64_QAM{33} = scaling_factor_64_QAM * [-7;7];
        Scaled_Demapper_64_QAM{34} = scaling_factor_64_QAM * [-5;7];
        Scaled_Demapper_64_QAM{35} = scaling_factor_64_QAM * [-7;5];
        Scaled_Demapper_64_QAM{36} = scaling_factor_64_QAM * [-5;5];        
        Scaled_Demapper_64_QAM{37} = scaling_factor_64_QAM * [-1;7];
        Scaled_Demapper_64_QAM{38} = scaling_factor_64_QAM * [-3;7];
        Scaled_Demapper_64_QAM{39} = scaling_factor_64_QAM * [-1;5];
        Scaled_Demapper_64_QAM{40} = scaling_factor_64_QAM * [-3;5];        
        Scaled_Demapper_64_QAM{41} = scaling_factor_64_QAM * [-7;1];
        Scaled_Demapper_64_QAM{42} = scaling_factor_64_QAM * [-5;1];
        Scaled_Demapper_64_QAM{43} = scaling_factor_64_QAM * [-7;3];
        Scaled_Demapper_64_QAM{44} = scaling_factor_64_QAM * [-5;3];        
        Scaled_Demapper_64_QAM{45} = scaling_factor_64_QAM * [-1;1];
        Scaled_Demapper_64_QAM{46} = scaling_factor_64_QAM * [-3;1];
        Scaled_Demapper_64_QAM{47} = scaling_factor_64_QAM * [-1;3];        
        Scaled_Demapper_64_QAM{48} = scaling_factor_64_QAM * [-3;3];        
        Scaled_Demapper_64_QAM{49} = scaling_factor_64_QAM * [7;7];
        Scaled_Demapper_64_QAM{50} = scaling_factor_64_QAM * [5;7];
        Scaled_Demapper_64_QAM{51} = scaling_factor_64_QAM * [7;5];
        Scaled_Demapper_64_QAM{52} = scaling_factor_64_QAM * [5;5];        
        Scaled_Demapper_64_QAM{53} = scaling_factor_64_QAM * [1;7];
        Scaled_Demapper_64_QAM{54} = scaling_factor_64_QAM * [3;7];
        Scaled_Demapper_64_QAM{55} = scaling_factor_64_QAM * [1;5];
        Scaled_Demapper_64_QAM{56} = scaling_factor_64_QAM * [3;5];        
        Scaled_Demapper_64_QAM{57} = scaling_factor_64_QAM * [7;1];
        Scaled_Demapper_64_QAM{58} = scaling_factor_64_QAM * [5;1];
        Scaled_Demapper_64_QAM{59} = scaling_factor_64_QAM * [7;3];
        Scaled_Demapper_64_QAM{60} = scaling_factor_64_QAM * [5;3];        
        Scaled_Demapper_64_QAM{61} = scaling_factor_64_QAM * [1;1];
        Scaled_Demapper_64_QAM{62} = scaling_factor_64_QAM * [1;3];
        Scaled_Demapper_64_QAM{63} = scaling_factor_64_QAM * [3;1];
        Scaled_Demapper_64_QAM{64} = scaling_factor_64_QAM * [3;3];
        
        % initialize x_hat
        x_hat=zeros(1, 2 * (size (s,2)));
        
        % initialize dist vector, to calculate the different distances
        dist=zeros(1,M);
        
        % Apply ML decoder
        bit_index=0;
        for temp=1:size(s,2)
            for temp_j=1:M
                dist(temp_j) = sqrt( (s(1,temp) - Scaled_Demapper_64_QAM{temp_j}(1) )^2 + (s(2,temp) - Scaled_Demapper_64_QAM{temp_j}(2) )^2  );
            end
            [~,min_index] = min(dist);
            demapped_bits_from_current_symbol=de2bi(min_index-1,k, 'left-msb');
            for i=1:k
                bit_index=bit_index+1;
                x_hat(bit_index)=demapped_bits_from_current_symbol(i);
            end                
        end
    else
        error('The given signal constellation is not supported. The supported constellations are: BPSK, QPSK, 16-QAM, 64-QAM');        
    end
    
end