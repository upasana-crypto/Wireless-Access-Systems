function s=mapping(x, constellationType)

    % the mapping function maps the bit sequence x, which contains elements
    % from the {0,1} to symbols using the desired constellation type.
    % Inputs:
    % 1. x: bit sequence
    % 2. constellation type: BPSK, QPSK, 16-QAM, 64-QAM
    % Output:
    % 1. s: symbol sequence
    
    
    % As mentioned in the exercise, average symbol energy is 1
    E=1;
    
    %% constellation type input control
    % M is the number of points in the signal constellation.     

    if strcmp(constellationType, 'BPSK')
               
        M=2;
        
        % k contains the number of mapped bits per symbol
        k = log2(M);
        
        % hardcode scaling factor
        scaling_factor_BPSK = sqrt(E);
        
        %% The following part of the code deals with the creation of the constellation map
        % This part is not essential for the exercise itself
        % but can help a great deal for debugging
        
        % hardcode the constellation 
%         Const_BPSK = scaling_factor_BPSK * [ -1, 1; 0, 0];
%         
%         figure(1);
%         plot(Const_BPSK(1,:),Const_BPSK(2,:), 'ko','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceColor','k', 'MarkerSize',5)
%         title('Scaled constellation for BPSK and Gray Code');
%         pl_BPSK = gca;
%         pl_BPSK.XAxisLocation = 'origin';
%         pl_BPSK.YAxisLocation = 'origin';
%         pl_BPSK.XLim = [-1.5 1.5];
%         pl_BPSK.YLim = [-1.5 1.5];
%         box off;
%         xticks([-1 1]);
%         yticks([-1 1]);
%         hold on;
        
        % hardcode gray code on the constellation
%         text(Const_BPSK(1,1)+.02,Const_BPSK(2,1)+.02,'0', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_BPSK(1,2)+.02,Const_BPSK(2,2)+.02,'1', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         hold off;
        
        
        %% definition of Mapper as a lookup table
        Mapper_BPSK{1} = -1;
        Mapper_BPSK{2} = 1;
        
        %% start of mapping process
        if ( mod(length(x),k) ~=0 )
            error("The given bitstream cannot be mapped entirely to symbols. Please change the length of input bitstream ");
        end
        
        % initialize s
        s_unscaled= zeros(1,length(x)/k);
        
        for temp=1:k:length(x)            
            % the last +1 is because matrices in matlab start from 1
            decimal_index= x(temp)+1;            
            s_unscaled(temp) = Mapper_BPSK{decimal_index};            
        end
        
        s= scaling_factor_BPSK * s_unscaled;

        
    elseif strcmp(constellationType, 'QPSK')
        M=4;
        
        % k contains the number of mapped bits per symbol
        k = log2(M);
        
        % hardcode scaling factor
        scaling_factor_QPSK = sqrt(E/2);
        
        %% The following part of the code deals with the creation of the constellation map
        % This part is not essential for the exercise itself
        % but can help a great deal for debugging
        
        % hardcode the constellation 
%         Const_QPSK = scaling_factor_QPSK * [ -1, 1, 1, -1; -1, -1, 1, 1];
%         
%         figure(2);
%         plot(Const_QPSK(1,:),Const_QPSK(2,:), 'ko','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceColor','k', 'MarkerSize',5)
%         title('Scaled constellation for QPSK and Gray Code');
%         pl_QPSK = gca;
%         pl_QPSK.XAxisLocation = 'origin';
%         pl_QPSK.YAxisLocation = 'origin';
%         pl_QPSK.XLim = [-1.5 1.5];
%         pl_QPSK.YLim = [-1.5 1.5];
%         box off;
%         xticks([-1 1]);
%         yticks([-1 1]);
%         hold on;
        
        % hardcode gray code on the constellation
%         text(Const_QPSK(1,1)+.02,Const_QPSK(2,1)+.02,'00', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_QPSK(1,2)+.02,Const_QPSK(2,2)+.02,'01', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_QPSK(1,3)+.02,Const_QPSK(2,3)+.02,'11', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_QPSK(1,4)+.02,Const_QPSK(2,4)+.02,'10', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         hold off;
        
        %% definition of Mapper
        Mapper_QPSK{1} = [-1;-1];
        Mapper_QPSK{2} = [1;-1];
        Mapper_QPSK{3} = [-1;1];
        Mapper_QPSK{4} = [1;1];
        
        %% start of mapping process
        if ( mod(length(x),k) ~=0 )
            error("The given bitstream cannot be mapped entirely to symbols. Please change the length of input bitstream ");
        end
        
        % initialize s
        s_unscaled= zeros(2,length(x)/k);
        
        for temp=1:k:length(x)            
            % the last +1 is because matrices in matlab start from 1
            decimal_index= x(temp)*2+x(temp+1)+1;            
            s_unscaled(:, fix(temp/k) +1) = Mapper_QPSK{decimal_index};            
        end
        
        s= scaling_factor_QPSK * s_unscaled;      
        
        
    elseif strcmp(constellationType, '16-QAM')
        M=16;
        
        % k contains the number of mapped bits per symbol
        k = log2(M);
        
        % hardcode scaling factor
        scaling_factor_16_QAM = sqrt(E/10);
        
        %% The following part of the code deals with the creation of the constellation map
        % This part is not essential for the exercise itself
        % but can help a great deal for debugging
        
        % hardcode the constellation 
%         Const_16_QAM = scaling_factor_16_QAM * [ -3, -1, -3, -1, -3, -1, -3, -1, 1, 3, 1, 3, 1, 3, 1, 3; -3, -3, -1, -1, 1, 1, 3, 3, -3, -3, -1, -1, 1, 1, 3, 3];
%         
%         figure(3);
%         plot(Const_16_QAM(1,:),Const_16_QAM(2,:), 'ko','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceColor','k', 'MarkerSize',5)
%         title('Scaled constellation for 16-QAM and Gray Code');
%         pl_16_QAM = gca;
%         pl_16_QAM.XAxisLocation = 'origin';
%         pl_16_QAM.YAxisLocation = 'origin';
%         pl_16_QAM.XLim = [-1.5 1.5];
%         pl_16_QAM.YLim = [-1.5 1.5];
%         box off;
%         xticks([-1 1]);
%         yticks([-1 1]);
%         hold on;
        
        % hardcode gray code on the constellation
%         text(Const_16_QAM(1,1)+.02,Const_16_QAM(2,1)+.02,'0000', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_16_QAM(1,2)+.02,Const_16_QAM(2,2)+.02,'0001', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_16_QAM(1,3)+.02,Const_16_QAM(2,3)+.02,'0010', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_16_QAM(1,4)+.02,Const_16_QAM(2,4)+.02,'0011', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         
%         text(Const_16_QAM(1,5)+.02,Const_16_QAM(2,5)+.02,'1010', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_16_QAM(1,6)+.02,Const_16_QAM(2,6)+.02,'1011', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_16_QAM(1,7)+.02,Const_16_QAM(2,7)+.02,'1000', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_16_QAM(1,8)+.02,Const_16_QAM(2,8)+.02,'1001', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         
%         text(Const_16_QAM(1,9)+.02,Const_16_QAM(2,9)+.02,'0101', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_16_QAM(1,10)+.02,Const_16_QAM(2,10)+.02,'0100', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_16_QAM(1,11)+.02,Const_16_QAM(2,11)+.02,'0111', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_16_QAM(1,12)+.02,Const_16_QAM(2,12)+.02,'0110', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         
%         text(Const_16_QAM(1,13)+.02,Const_16_QAM(2,13)+.02,'1111', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_16_QAM(1,14)+.02,Const_16_QAM(2,14)+.02,'1110', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_16_QAM(1,15)+.02,Const_16_QAM(2,15)+.02,'1101', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_16_QAM(1,16)+.02,Const_16_QAM(2,16)+.02,'1100', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         hold off;
        
        %% definition of Mapper
        Mapper_16_QAM{1} = [-3;-3];
        Mapper_16_QAM{2} = [-1;-3];
        Mapper_16_QAM{3} = [-3;-1];
        Mapper_16_QAM{4} = [-1;-1];
        Mapper_16_QAM{5} = [3;-3];
        Mapper_16_QAM{6} = [1;-3];
        Mapper_16_QAM{7} = [3;-1];
        Mapper_16_QAM{8} = [1;-1];
        Mapper_16_QAM{9} = [-3;3];
        Mapper_16_QAM{10} = [-1;3];
        Mapper_16_QAM{11} = [-3;1];
        Mapper_16_QAM{12} = [-1;1];
        Mapper_16_QAM{13} = [3;3];
        Mapper_16_QAM{14} = [1;3];
        Mapper_16_QAM{15} = [3;1];
        Mapper_16_QAM{16} = [1;1];
        
        %% start of mapping process
        if ( mod(length(x),k) ~=0 )
            error("The given bitstream cannot be mapped entirely to symbols. Please change the length of input bitstream ");
        end
        
        % initialize s
        s_unscaled= zeros(2,length(x)/k);
        
        for temp=1:k:length(x)
            
            % the last +1 is because matrices in matlab start from 1
            decimal_index= x(temp)*8+x(temp+1)*4+x(temp+2)*2+x(temp+3)+1;            
            s_unscaled(:, fix(temp/k) +1) = Mapper_16_QAM{decimal_index};            
        end
        
        s= scaling_factor_16_QAM * s_unscaled;
        
        
    elseif strcmp(constellationType, '64-QAM')
        M=64;
        
        % k contains the number of mapped bits per symbol
        k = log2(M);
        
        % hardcode scaling factor
        scaling_factor_64_QAM = sqrt(E/42);
        
        %% The following part of the code deals with the creation of the constellation map
        % This part is not essential for the exercise itself
        % but can help a great deal for debugging
%         Const_64_QAM = scaling_factor_64_QAM * [ -7, -5, -3, -1, -1, -3, -5, -7, -7, -5, -3, -1, -1, -3, -5, -7, -7, -5, -3, -1, -1, -3, -5, -7, -7, -5, -3, -1, -1, -3, -5, -7, 1, 3, 5, 7, 7, 5, 3, 1, 1, 3, 5, 7, 7, 5, 3, 1, 1, 3, 5, 7, 7, 5, 3, 1, 1, 3, 5, 7, 7, 5, 3, 1;
%                                                  -7, -7, -7, -7, -5, -5, -5, -5, -3, -3, -3, -3, -1, -1, -1, -1, 1, 1, 1, 1, 3, 3, 3, 3, 5, 5, 5, 5, 7, 7, 7, 7, -7, -7, -7, -7, -5, -5, -5, -5, -3, -3, -3, -3, -1, -1, -1, -1, 1, 1, 1, 1, 3, 3, 3, 3, 5, 5, 5, 5, 7, 7, 7, 7];
%         
%         figure(4);
%         plot(Const_64_QAM(1,:),Const_64_QAM(2,:), 'ko','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceColor','k', 'MarkerSize',5)
%         title('Scaled constellation for 64-QAM and Gray Code');
%         pl_64_QAM = gca;
%         pl_64_QAM.XAxisLocation = 'origin';
%         pl_64_QAM.YAxisLocation = 'origin';
%         pl_64_QAM.XLim = [-1.5 1.5];
%         pl_64_QAM.YLim = [-1.5 1.5];
%         box off;
%         xticks([-1 1]);
%         yticks([-1 1]);
%         hold on;
%         
%         text(Const_64_QAM(1,1)+.02,Const_64_QAM(2,1)+.02,'000000', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,2)+.02,Const_64_QAM(2,2)+.02,'000001', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,3)+.02,Const_64_QAM(2,3)+.02,'000101', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,4)+.02,Const_64_QAM(2,4)+.02,'000100', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');        
%         text(Const_64_QAM(1,5)+.02,Const_64_QAM(2,5)+.02,'000110', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,6)+.02,Const_64_QAM(2,6)+.02,'000111', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,7)+.02,Const_64_QAM(2,7)+.02,'000011', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,8)+.02,Const_64_QAM(2,8)+.02,'000010', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');        
%         text(Const_64_QAM(1,9)+.02,Const_64_QAM(2,9)+.02,'001010', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,10)+.02,Const_64_QAM(2,10)+.02,'001011', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,11)+.02,Const_64_QAM(2,11)+.02,'001111', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,12)+.02,Const_64_QAM(2,12)+.02,'001110', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');        
%         text(Const_64_QAM(1,13)+.02,Const_64_QAM(2,13)+.02,'001100', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,14)+.02,Const_64_QAM(2,14)+.02,'001101', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,15)+.02,Const_64_QAM(2,15)+.02,'001001', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,16)+.02,Const_64_QAM(2,16)+.02,'001000', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         
%         text(Const_64_QAM(1,17)+.02,Const_64_QAM(2,17)+.02,'101000', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,18)+.02,Const_64_QAM(2,18)+.02,'101001', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,19)+.02,Const_64_QAM(2,19)+.02,'101101', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,20)+.02,Const_64_QAM(2,20)+.02,'101100', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');        
%         text(Const_64_QAM(1,21)+.02,Const_64_QAM(2,21)+.02,'101110', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,22)+.02,Const_64_QAM(2,22)+.02,'101111', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,23)+.02,Const_64_QAM(2,23)+.02,'101011', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,24)+.02,Const_64_QAM(2,24)+.02,'101010', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');        
%         text(Const_64_QAM(1,25)+.02,Const_64_QAM(2,25)+.02,'100010', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,26)+.02,Const_64_QAM(2,26)+.02,'100011', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,27)+.02,Const_64_QAM(2,27)+.02,'100111', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,28)+.02,Const_64_QAM(2,28)+.02,'100110', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');        
%         text(Const_64_QAM(1,29)+.02,Const_64_QAM(2,29)+.02,'100100', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,30)+.02,Const_64_QAM(2,30)+.02,'100101', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,31)+.02,Const_64_QAM(2,31)+.02,'100001', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,32)+.02,Const_64_QAM(2,32)+.02,'100000', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         
%         text(Const_64_QAM(1,33)+.02,Const_64_QAM(2,33)+.02,'010100', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,34)+.02,Const_64_QAM(2,34)+.02,'010101', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,35)+.02,Const_64_QAM(2,35)+.02,'010001', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,36)+.02,Const_64_QAM(2,36)+.02,'010000', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');        
%         text(Const_64_QAM(1,37)+.02,Const_64_QAM(2,37)+.02,'010010', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,38)+.02,Const_64_QAM(2,38)+.02,'010011', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,39)+.02,Const_64_QAM(2,39)+.02,'010111', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,40)+.02,Const_64_QAM(2,40)+.02,'010110', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');        
%         text(Const_64_QAM(1,41)+.02,Const_64_QAM(2,41)+.02,'011110', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,42)+.02,Const_64_QAM(2,42)+.02,'011111', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,43)+.02,Const_64_QAM(2,43)+.02,'011011', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,44)+.02,Const_64_QAM(2,44)+.02,'011010', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');        
%         text(Const_64_QAM(1,45)+.02,Const_64_QAM(2,45)+.02,'011000', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,46)+.02,Const_64_QAM(2,46)+.02,'011001', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,47)+.02,Const_64_QAM(2,47)+.02,'011101', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,48)+.02,Const_64_QAM(2,48)+.02,'011100', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         
%         text(Const_64_QAM(1,49)+.02,Const_64_QAM(2,49)+.02,'111100', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,50)+.02,Const_64_QAM(2,50)+.02,'111101', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,51)+.02,Const_64_QAM(2,51)+.02,'111001', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,52)+.02,Const_64_QAM(2,52)+.02,'111000', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');        
%         text(Const_64_QAM(1,53)+.02,Const_64_QAM(2,53)+.02,'111010', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,54)+.02,Const_64_QAM(2,54)+.02,'111011', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,55)+.02,Const_64_QAM(2,55)+.02,'111111', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,56)+.02,Const_64_QAM(2,56)+.02,'111110', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');        
%         text(Const_64_QAM(1,57)+.02,Const_64_QAM(2,57)+.02,'110110', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,58)+.02,Const_64_QAM(2,58)+.02,'110111', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,59)+.02,Const_64_QAM(2,59)+.02,'110011', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,60)+.02,Const_64_QAM(2,60)+.02,'110010', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');        
%         text(Const_64_QAM(1,61)+.02,Const_64_QAM(2,61)+.02,'110000', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,62)+.02,Const_64_QAM(2,62)+.02,'110001', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,63)+.02,Const_64_QAM(2,63)+.02,'110101', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%         text(Const_64_QAM(1,64)+.02,Const_64_QAM(2,64)+.02,'110100', 'horiz','center','vert','bottom', 'Fontsize', 12, 'Color', 'b');
%               
        %% definition of Mapper
        Mapper_64_QAM{1} = [-7;-7];
        Mapper_64_QAM{2} = [-5;-7];
        Mapper_64_QAM{3} = [-7;-5];
        Mapper_64_QAM{4} = [-5;-5];
        Mapper_64_QAM{5} = [-1;-7];
        Mapper_64_QAM{6} = [-3;-7];
        Mapper_64_QAM{7} = [-1;-5];
        Mapper_64_QAM{8} = [-3;-5];
        Mapper_64_QAM{9} = [-7;-1];
        Mapper_64_QAM{10} = [-5;-1];
        Mapper_64_QAM{11} = [-7;-3];
        Mapper_64_QAM{12} = [-5;-3];
        Mapper_64_QAM{13} = [-1;-1];
        Mapper_64_QAM{14} = [-3;-1];
        Mapper_64_QAM{15} = [-1;-3];
        Mapper_64_QAM{16} = [-3;-3];        
        Mapper_64_QAM{17} = [7;-7];
        Mapper_64_QAM{18} = [5;-7];
        Mapper_64_QAM{19} = [7;-5];
        Mapper_64_QAM{20} = [5;-5];        
        Mapper_64_QAM{21} = [1;-7];
        Mapper_64_QAM{22} = [3;-7];
        Mapper_64_QAM{23} = [1;-5];
        Mapper_64_QAM{24} = [3;-5];        
        Mapper_64_QAM{25} = [7;-1];
        Mapper_64_QAM{26} = [5;-1];
        Mapper_64_QAM{27} = [7;-3];
        Mapper_64_QAM{28} = [5;-3];        
        Mapper_64_QAM{29} = [1;-1];
        Mapper_64_QAM{30} = [3;-1];
        Mapper_64_QAM{31} = [1;-3];
        Mapper_64_QAM{32} = [3;-3];        
        Mapper_64_QAM{33} = [-7;7];
        Mapper_64_QAM{34} = [-5;7];
        Mapper_64_QAM{35} = [-7;5];
        Mapper_64_QAM{36} = [-5;5];        
        Mapper_64_QAM{37} = [-1;7];
        Mapper_64_QAM{38} = [-3;7];
        Mapper_64_QAM{39} = [-1;5];
        Mapper_64_QAM{40} = [-3;5];        
        Mapper_64_QAM{41} = [-7;1];
        Mapper_64_QAM{42} = [-5;1];
        Mapper_64_QAM{43} = [-7;3];
        Mapper_64_QAM{44} = [-5;3];        
        Mapper_64_QAM{45} = [-1;1];
        Mapper_64_QAM{46} = [-3;1];
        Mapper_64_QAM{47} = [-1;3];        
        Mapper_64_QAM{48} = [-3;3];        
        Mapper_64_QAM{49} = [7;7];
        Mapper_64_QAM{50} = [5;7];
        Mapper_64_QAM{51} = [7;5];
        Mapper_64_QAM{52} = [5;5];        
        Mapper_64_QAM{53} = [1;7];
        Mapper_64_QAM{54} = [3;7];
        Mapper_64_QAM{55} = [1;5];
        Mapper_64_QAM{56} = [3;5];        
        Mapper_64_QAM{57} = [7;1];
        Mapper_64_QAM{58} = [5;1];
        Mapper_64_QAM{59} = [7;3];
        Mapper_64_QAM{60} = [5;3];        
        Mapper_64_QAM{61} = [1;1];
        Mapper_64_QAM{62} = [1;3];
        Mapper_64_QAM{63} = [3;1];
        Mapper_64_QAM{64} = [3;3];
        
         %% start of mapping process
        if ( mod(length(x),k) ~=0 )
            error('The given bitstream cannot be mapped entirely to symbols. Please change the length of input bitstream');
        end
        
        % initialize s
        s_unscaled= zeros(2,length(x)/k);
        
        for temp=1:k:length(x)            
            % the last +1 is because matrices in matlab start from 1
            decimal_index= x(temp)*32+x(temp+1)*16+x(temp+2)*8+x(temp+3)*4+x(temp+4)*2+x(temp+5)+1;            
            s_unscaled(:, fix(temp/k) +1) = Mapper_64_QAM{decimal_index};            
        end
        
        s= scaling_factor_64_QAM * s_unscaled;
    else
        error('The given signal constellation is not supported. The supported constellations are: BPSK, QPSK, 16-QAM, 64-QAM');
    end  
 
end