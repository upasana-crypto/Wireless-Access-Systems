clear all;
close all;
clc;

tic;
% N: the amount of symbols in each block
% N is fixed from the exercise 
N = 256;

% R is the rate
% R is fixed from the exercise
R_nominator = 1;
R_denominator = 2;
R = R_nominator/R_denominator;

% an array with all the memory depths of the encoder
% the branch impulse responses are obtained from Table1
m_array = [2,4,6];

% the number of branch impulse responses will be R_denominator
% for example:
g_octal_matr = zeros(length(m_array),R_denominator);
for index_m_array=1:length(m_array)
    if index_m_array == 1
        g_octal_matr(index_m_array,1) = 5;
        g_octal_matr(index_m_array,2) = 7;
    elseif index_m_array == 2
        g_octal_matr(index_m_array,1) = 27;
        g_octal_matr(index_m_array,2) = 31;
    elseif index_m_array == 3
        g_octal_matr(index_m_array,1) = 117;
        g_octal_matr(index_m_array,2) = 155;
    end
end

g_cell = cell(1,length(m_array));

for index_g_cell=1:length(g_cell)
    g_per_m = zeros(R_denominator,m_array(index_g_cell)+1);
    for index=1:R_denominator
        g_per_m(index,:) = oct2poly(g_octal_matr(index_g_cell,index));        
    end
    g_cell{index_g_cell} = g_per_m;
end

% mod_scheme contains which modulation schemes will be simulated
%mod_scheme = {'BPSK'};
%mod_scheme = {'QPSK'};
%mod_scheme = {'16-QAM'};
%mod_scheme = {'64-QAM'};
%mod_scheme = {'BPSK','QPSK'};
mod_scheme = {'BPSK', '64-QAM'};
%mod_scheme = {'BPSK', 'QPSK', '16-QAM', '64-QAM'};

% BER_coded_all_schemes will contain a [N_symbols_sim X Eb_No] matrix for EACH
% modulation scheme when coding is used
BER_coded_all_schemes = cell(length(m_array),length(mod_scheme));

% BER_uncoded_all_schemes will contain a [N_symbols_sim X Eb_No] matrix for EACH 
% modulation scheme when NO-coding is used
BER_uncoded_all_schemes = cell(1,length(mod_scheme));

% range of Eb_No in dB that will be used in the simulation
Eb_No_dB = -20:1.5:20;

% range of Eb_No in decimal values that will be used in simulation
Eb_No = 10.^(Eb_No_dB/10);

% iterate for N_symbols_sim to achieve more accurate statistical
% results for each Eb/No value
N_symbols_sim = 10000;

% iterate across all values of m
for index_m_array = 1: length(m_array)

    % extract the matrix containing the Branch Impulse Responses from cell
    temp_g_matrix = g_cell{index_m_array};
    
    % iterate across all modulation schemes
    for index = 1:length(mod_scheme)
        if strcmp(mod_scheme(index),'BPSK')
            % for BPSK --> a=1
            a = 1;        
        elseif strcmp(mod_scheme{index},'QPSK')        
            % for QPSK --> a=2
            a = 2;        
        elseif strcmp(cellstr(mod_scheme(index)),'16-QAM')        
            % for 16-QAM --> a=4
            a = 4;        
        elseif strcmp(mod_scheme(index),'64-QAM')        
            % for 64-QAM --> a=6
            a = 6;        
        else
            error('modulation scheme not supported');
        end

        % calculate Nb
        Nb = N*a*R -m_array(index_m_array);

        % Initialize the BER matrix
        % this matrix will contain all the BER's across all symbols and Eb/No
        % values
        BER_per_mod_scheme_coded = zeros(N_symbols_sim,length(Eb_No_dB));

        % iterate across different SNR values
        for index_Eb_No=1:length(Eb_No_dB)

            % iterate for all N_symbols_sim
            for index_symbol=1:N_symbols_sim

                % generate bit sequence that is the MESSAGE
                x_info = randi([0 1],1,Nb);

                % append m trailing bits at the end
                x = [x_info zeros(1,m_array(index_m_array)) ];

                %% encoding with the convolution encoder

                % 1. Initializing
                internal_state = zeros(1,m_array(index_m_array));
                x_encoded_bits = zeros(1,length(x)*R_denominator);

                % 2. for each input bit
                for i=1:length(x)

                    % load new bit input
                    new_bit = x(i);

                    % current state of encoder
                    current_state = [new_bit internal_state];

                    %% Encoder implementation

                    % calculate the output per branch and output it
                    for branch=1:R_denominator
                        operands = current_state.*temp_g_matrix(branch,:);
                        temp = operands(1);
                        for index2=2:length(operands)
                            temp = mod(temp+operands(index2),2);            
                        end
                        x_encoded_bits(R_denominator*i-1*mod(branch,R_denominator)) = temp;            
                    end      

                    %% update internal state
                    for k=(m_array(index_m_array)-1):-1:1
                        internal_state(k+1) = internal_state(k);
                    end
                    internal_state(1) = new_bit;

                end

                %% Interleaver implementation

                % B stands for the dimension of the interleaver
                N_block = 16;
                M_block = 16;

                % test how many times the interleaver must be used and if padding is
                % needed
                number_of_interleaving = floor(length(x_encoded_bits)/(N_block*M_block));
                remainder_of_interleaving = mod(length(x_encoded_bits),N_block*M_block);

                if remainder_of_interleaving == 0
                    x_interleaved_encoded = zeros(1,length(x_encoded_bits));
                else
                    error("The block interleaver couldn't work for the length of the input bitstream");
                end

                % pass the x_encoded_bits through the interleaver block by block
                for index3 = 1: number_of_interleaving

                    % writing into the interleaver
                    interleaver = reshape(x_encoded_bits(1+(index3-1)*N_block*M_block:index3*N_block*M_block),[N_block,M_block])';

                    % output of the interleaver
                    x_interleaved_encoded(1+(index3-1)*N_block*M_block:index3*N_block*M_block) = reshape(interleaver,[1,N_block*M_block]);        
                end
                %% Mapper for coded bits 

                % coded bits: call mapping function
                symbol_sequence_mapped_coded = mapping(x_interleaved_encoded,mod_scheme(index));

                %% Adding noise for CODED bits

                % 1. calculating average symbol energy for the coded stream
                sum=0;
                if strcmp(mod_scheme(index),'BPSK')
                    for index5 =1:N
                        sum = sum + symbol_sequence_mapped_coded(index5)^2;   
                    end
                else
                    for index5 =1:N
                        sum = sum + symbol_sequence_mapped_coded(1,index5)^2 + symbol_sequence_mapped_coded(2,index5)^2;
                    end
                end
                E_symb_avg = sum/N;

                % 2. calculate Eb
                Eb = E_symb_avg/a;

                % 3. calculate noise variance
                sigma_square = Eb/Eb_No(index_Eb_No);

                % 4. generate noise
                noise = sqrt(sigma_square/2)*randn(size(symbol_sequence_mapped_coded));

                % 5. adding noise
                symbol_sequence_noisy_coded = symbol_sequence_mapped_coded + noise;

                %% Demapper for CODED bits

                % coded bits: calling demapping function
                x_demapped_coded= demapping(symbol_sequence_noisy_coded,mod_scheme(index));

                %% De-interleaver implementation

                % pass the interleaved bits through the de-interleaver block by block
                for index4 = 1: number_of_interleaving

                    % writting into the de-interleaver
                    de_inteleaver = reshape(x_demapped_coded((1+(index4-1)*N_block*M_block:index4*N_block*M_block)),[N_block,M_block]);

                    % reading from the de-interleaver
                    x_de_interleaved(1+(index4-1)*N_block*M_block:index4*N_block*M_block) = reshape(de_inteleaver',[1,N_block*M_block]);

                end

                %% Decoding using a Viterbi decoder (for the CODED bits)
                % using the code from the exercise

                % Generate the trellis from the polynomials
                tr = poly2trellis( [m_array(index_m_array)+1],[g_octal_matr(index_m_array,1) g_octal_matr(index_m_array,2)]);

                % set the traceback length. If the code is 1/2, a typical value for it
                % is about 5 times m
                tracebacklength = 5*m_array(index_m_array);

                % call the vitdec function.
                x_decoded = vitdec(x_de_interleaved,tr,tracebacklength,'term','hard');

                % Attention: the output of vitdec is the decoded block of bits
                % INCLUDING the trailing bits

                % drop the m trailing bits
                x_estimated_info = x_decoded(1:length(x_info));

                % count bit error rate in the current block
                errors_in_current_block_symbol_coded = 0;
                for index6=1:length(x_info)
                    if x_estimated_info(index6) ~= x_info(index6)
                        errors_in_current_block_symbol_coded = errors_in_current_block_symbol_coded + 1;
                    end
                end
                BER_per_mod_scheme_coded(index_symbol,index_Eb_No) = errors_in_current_block_symbol_coded/length(x_info);          

            end
        end
        BER_coded_all_schemes{index_m_array,index} = BER_per_mod_scheme_coded;
    end

    %% UNCODED BITS section --> this part should be executed ONLY once
    if index_m_array == 1
        
        % iterate across all modulation schemes
        for index = 1:length(mod_scheme)
            if strcmp(mod_scheme(index),'BPSK')
                % for BPSK --> a=1
                a = 1;        
            elseif strcmp(mod_scheme{index},'QPSK')        
                % for QPSK --> a=2
                a = 2;        
            elseif strcmp(cellstr(mod_scheme(index)),'16-QAM')        
                % for 16-QAM --> a=4
                a = 4;        
            elseif strcmp(mod_scheme(index),'64-QAM')        
                % for 64-QAM --> a=6
                a = 6;        
            else
                error('modulation scheme not supported');
            end
        
            Nb_uncoded = N*a;          

            BER_per_mod_scheme_uncoded = zeros(N_symbols_sim,length(Eb_No_dB));
            % iterate across different SNR values
            for index_Eb_No=1:length(Eb_No_dB)

                % iterate for all N_symbols_sim
                for index_symbol=1:N_symbols_sim
                    
                    x_info_uncoded = randi([0 1],1,Nb_uncoded);
            
                    %% Mapper for uncoded bits
            
                    % uncoded_bits: call mapping function
                    symbol_sequence_mapped_uncoded = mapping(x_info_uncoded,mod_scheme(index));
            
                    %%  Adding noise for UNCODED bits
 
                    % 1. calculating average symbol energy for the coded stream
                    sum_uncoded = 0;
                    if strcmp(mod_scheme(index),'BPSK')
                        for index5_unc =1:N
                            sum_uncoded = sum_uncoded + symbol_sequence_mapped_uncoded(index5_unc)^2;   
                        end
                    else
                        for index5_unc =1:N
                            sum_uncoded = sum_uncoded + symbol_sequence_mapped_uncoded(1,index5_unc)^2 + symbol_sequence_mapped_uncoded(2,index5_unc)^2;
                        end
                    end
                    
                    E_symb_avg_unc = sum_uncoded/N;
 
                    % 2. calculate Eb
                    Eb_unc = E_symb_avg_unc/a;
 
                    % 3. calculate noise variance
                    sigma_square_unc= Eb_unc/Eb_No(index_Eb_No);
 
                    % 4. generate noise
                    noise_unc = sqrt(sigma_square_unc/2)*randn(size(symbol_sequence_mapped_uncoded));
 
                    % 5. adding noise
                    symbol_sequence_noisy_uncoded = symbol_sequence_mapped_uncoded + noise_unc;
            
                    %% Demapping for uncoded bits
            
                    % uncoded bits: calling demapping function
                    x_demapped_uncoded = demapping(symbol_sequence_noisy_uncoded,mod_scheme(index));
                         
                    %% BER calculation for the Uncoded bits 
                    errors_in_current_block_symbol_uncoded = 0;
                    for index6_uncoded=1:length(x_info_uncoded)
                        if x_demapped_uncoded(index6_uncoded) ~= x_info_uncoded(index6_uncoded)
                            errors_in_current_block_symbol_uncoded = errors_in_current_block_symbol_uncoded +1;
                        end
                    end
                    BER_per_mod_scheme_uncoded(index_symbol,index_Eb_No) = errors_in_current_block_symbol_uncoded /length(x_info_uncoded);
                end
            end
            BER_uncoded_all_schemes{index} = BER_per_mod_scheme_uncoded;
        end      
    end     
end

%% Data Post-Processing for CODED bits --> find the average of BER across all symbols for each  mod scheme
Avg_BER_coded = cell(length(m_array),length(mod_scheme));

for index_m_array=1:length(m_array)
    for index=1:length(mod_scheme)
         Avg_BER_coded{index_m_array,index} = mean(BER_coded_all_schemes{index_m_array,index},1);        
    end
end

%% Data Post-Processing for UNCODED bits 
Avg_BER_uncoded = cell(1,length(mod_scheme));
for index=1:length(mod_scheme)
    Avg_BER_uncoded{index} = mean(BER_uncoded_all_schemes{index},1);
end

toc;
%% Plots

figure()
semilogy(Eb_No_dB, Avg_BER_coded{1,1}, 'b--*', 'LineWidth', 3, 'MarkerSize',10)
hold on;
semilogy(Eb_No_dB, Avg_BER_coded{2,1}, 'r--+', 'LineWidth', 3, 'MarkerSize',10);
hold on;
semilogy(Eb_No_dB, Avg_BER_coded{3,1}, 'g--o', 'LineWidth', 3, 'MarkerSize',10);
hold on;
title('Comparison across encoders with different memory depth (m) for BPSK');
xlabel('Eb/No in dB');
ylabel('BER in log scale');
xlim([min(Eb_No_dB) max(Eb_No_dB)]);
ylim([5*10^(-log10(N_symbols_sim)) 1]);
ax = gca;
set(gca,'FontSize',18);
ax.GridAlpha = 1;
ax.MinorGridAlpha = 0.55;
grid on;
leg = legend({'m = 2','m = 4','m = 6'},'Location','southwest');
set(leg,'FontSize',18);
hold off;

figure()
semilogy(Eb_No_dB, Avg_BER_coded{1,1}, 'b--*', 'LineWidth', 3, 'MarkerSize',10)
hold on;
semilogy(Eb_No_dB, Avg_BER_coded{2,1}, 'r--+', 'LineWidth', 3, 'MarkerSize',10);
hold on;
semilogy(Eb_No_dB, Avg_BER_coded{3,1}, 'g--o', 'LineWidth', 3, 'MarkerSize',10);
hold on;
semilogy(Eb_No_dB, Avg_BER_uncoded{1}, 'k--d', 'LineWidth', 3, 'MarkerSize',10);
title('Comparison across encoders with different memory depth (m) and uncoded for BPSK');
xlabel('Eb/No in dB');
ylabel('BER in log scale');
xlim([min(Eb_No_dB) max(Eb_No_dB)]);
ylim([5*10^(-log10(N_symbols_sim)) 1]);
ax = gca;
set(gca,'FontSize',18);
ax.GridAlpha = 1;
ax.MinorGridAlpha = 0.55;
grid on;
leg = legend({'m = 2','m = 4','m = 6','uncoded'},'Location','southwest');
set(leg,'FontSize',18);
hold off;

figure()
semilogy(Eb_No_dB, Avg_BER_coded{1,2}, 'b--*', 'LineWidth', 3, 'MarkerSize',10)
hold on;
semilogy(Eb_No_dB, Avg_BER_coded{2,2}, 'r--+', 'LineWidth', 3, 'MarkerSize',10);
hold on;
semilogy(Eb_No_dB, Avg_BER_coded{3,2}, 'g--o', 'LineWidth', 3, 'MarkerSize',10);
hold on;
title('Comparison across encoders with different memory depth (m) for 64-QAM');
xlabel('Eb/No in dB');
ylabel('BER in log scale');
xlim([min(Eb_No_dB) max(Eb_No_dB)]);
ylim([5*10^(-log10(N_symbols_sim)) 1]);
ax = gca;
set(gca,'FontSize',18);
ax.GridAlpha = 1;
ax.MinorGridAlpha = 0.55;
grid on;
leg = legend({'m = 2','m = 4','m = 6'},'Location','southwest');
set(leg,'FontSize',18);
hold off;

figure()
semilogy(Eb_No_dB, Avg_BER_coded{1,2}, 'b--*', 'LineWidth', 3, 'MarkerSize',10)
hold on;
semilogy(Eb_No_dB, Avg_BER_coded{2,2}, 'r--+', 'LineWidth', 3, 'MarkerSize',10);
hold on;
semilogy(Eb_No_dB, Avg_BER_coded{3,2}, 'g--o', 'LineWidth', 3, 'MarkerSize',10);
hold on;
semilogy(Eb_No_dB, Avg_BER_uncoded{2}, 'k--d', 'LineWidth', 3, 'MarkerSize',10);
title('Comparison across encoders with different memory depth (m) and uncoded for 64-QAM');
xlabel('Eb/No in dB');
ylabel('BER in log scale');
xlim([min(Eb_No_dB) max(Eb_No_dB)]);
ylim([5*10^(-log10(N_symbols_sim)) 1]);
ax = gca;
set(gca,'FontSize',18);
ax.GridAlpha = 1;
ax.MinorGridAlpha = 0.55;
grid on;
leg = legend({'m = 2','m = 4','m = 6','uncoded'},'Location','southwest');
hold off;