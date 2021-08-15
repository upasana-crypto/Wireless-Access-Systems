clear all;
close all;
clc;

% start the timer
tic;

%% start of figure 6 implementation

% # of simulations
N_sim = 100;

% #of symbols in OFDM
N = 256;

% fft / ifft points
Nfft = N;

N_data = 8*N;
N_train = 2*N;

% cyclic prefix
Ncp = 10;

% range of Eb/No
Eb_N0_dB = -20:1.5:30;
Eb_N0_dec = 10.^(Eb_N0_dB/10);

%% multi-path channel
L = 10;
lamda = 1/5;

% draw random coefficients from Gaussian distribution
h_coeff = sqrt(1/2)*randn(1,L)+1i*sqrt(1/2)*randn(1,L);

% construct h
h = zeros(1,length(h_coeff));
for ind_h = 1:length(h_coeff)
    h(ind_h)= h(ind_h)+h_coeff(ind_h)*exp(-lamda*(ind_h-1));    
end

% calculate power of h
P_h = power_calc(h);

% from exercise --> power(h) =1
P_h_wanted = 1;

% construct h_scaled
h_sc = zeros(size(h));
for ind_h_sc=1: length(h)
    h_sc(ind_h_sc) = sqrt(P_h_wanted/P_h) * h(ind_h_sc);
end

% verify that power(h_sc) = 1
P_h_sc = power_calc(h_sc);

% calculate H for debugging purposes
h_usc_fft = fft(h_sc,Nfft);

% calculate power of h_fft unscaled
P_h_usc_fft = power_calc(h_usc_fft);

% scale h_fft
h_sc_fft = sqrt(P_h_sc/P_h_usc_fft)*h_usc_fft;

% verify that scaling succeded
P_h_sc_fft = power_calc(h_sc_fft);

%% interleaver parameters
N_block = 16;
M_block = 16;

%% choice of the encoder
% R is the rate
R_nominator = 1;
R_denominator = 2;
R = R_nominator/R_denominator;

% m is the memory depth of the encoder
m = 4;

% the number of branch impulse responses will be R_denominator
g_octal_matr = zeros(1,R_denominator);
g_octal_matr(1) = 27;
g_octal_matr(2) = 31;

g_matrix = zeros(R_denominator,m+1);
for index=1:R_denominator
    g_matrix(index,:) = oct2poly(g_octal_matr(index));
end

mod_schemes = {'BPSK','QPSK','16-QAM','64-QAM'};
%mod_schemes = {'16-QAM'};
%mod_schemes = {'BPSK'};
%mod_schemes = {'BPSK','QPSK'};

BER_coded = zeros(length(mod_schemes), length(Eb_N0_dB));
BER_uncoded = zeros(length(mod_schemes), length(Eb_N0_dB));
BER_coded_NI = zeros(length(mod_schemes), length(Eb_N0_dB));

for ind_m_s = 1: length(mod_schemes)
    
    if strcmp(mod_schemes(ind_m_s),'BPSK')
        % for BPSK --> a=1
        a = 1;        
    elseif strcmp(mod_schemes{ind_m_s},'QPSK')        
        % for QPSK --> a=2
        a = 2;        
    elseif strcmp(cellstr(mod_schemes(ind_m_s)),'16-QAM')        
        % for 16-QAM --> a=4
        a = 4;        
    elseif strcmp(mod_schemes(ind_m_s),'64-QAM')        
        % for 64-QAM --> a=6
        a = 6;        
    else
        error('modulation scheme not supported');
    end
    
    % calculate Nb for coded scenario
    Nb_coded_data = N_data*a*R -m;
    
    % calculate Nb for the uncoded scenario
    Nb_uncoded_data = N_data*a;
    
    % initialization
    ber_coded_interl = zeros(1,N_sim);
    ber_uncoded = zeros(1,N_sim);
    ber_coded_NI = zeros(1,N_sim);
    
    for ind_Eb_N0_dB=1:length(Eb_N0_dB)
        
        for ind_N_sim=1:N_sim
            
            %% coded & interleaver version
            
            % generate bit sequence that is the MESSAGE
            x_in_cod = randi([0 1],1,Nb_coded_data);
            
            % append m trailing bits at the end
            x_in_cod_TrailBits = [x_in_cod zeros(1,m) ];
            
            % encoding with the convolution encoder
            % 1. Initializing
            internal_state = zeros(1,m);
            x_encoded_bits = zeros(1,length(x_in_cod_TrailBits)*R_denominator);
            
            % 2. for each input bit
            for i=1:length(x_in_cod_TrailBits)
                
                % load new bit input
                new_bit = x_in_cod_TrailBits(i);
                
                % current state of encoder
                current_state = [new_bit internal_state];
                
                % Encoder implementation
                
                % calculate the output per branch and output it
                for branch=1:R_denominator
                    operands = current_state.*g_matrix(branch,:);
                    temp = operands(1);
                    for index2=2:length(operands)
                        temp = mod(temp+operands(index2),2);
                    end
                    x_encoded_bits(R_denominator*i-1*mod(branch,R_denominator)) = temp;
                end
                
                % update internal state
                for k=(m-1):-1:1
                    internal_state(k+1) = internal_state(k);
                end
                internal_state(1) = new_bit;
                
            end
            
            % Interleaver implementation
            
            % test how many times the interleaver must be used and if padding is
            % needed
            number_of_interleaving = floor(length(x_encoded_bits)/(N_block*M_block));
            remainder_of_interleaving = mod(length(x_encoded_bits),N_block*M_block);
            
            if remainder_of_interleaving == 0
                x_cod_intrl = zeros(1,length(x_encoded_bits));
            else
                error("The block interleaver couldn't work for the length of the input bitstream");
            end
            
            % pass the x_encoded_bits through the interleaver block by block
            for ind_intrl = 1: number_of_interleaving
                
                % writing into the interleaver
                interleaver = reshape(x_encoded_bits(1+(ind_intrl-1)*N_block*M_block:ind_intrl*N_block*M_block),[N_block,M_block])';
                
                % output of the interleaver
                x_cod_intrl(1+(ind_intrl-1)*N_block*M_block:ind_intrl*N_block*M_block) = reshape(interleaver,[1,N_block*M_block]);
            end
            
            % mapper
            s_data_mapped_vec = mapping(x_cod_intrl,mod_schemes(ind_m_s));
            
            % prepare to write symbols as complex numbers
            if strcmp(mod_schemes{ind_m_s},'BPSK')
                % in BPSK case nothing has to be done
                s_data_mapped_row = s_data_mapped_vec;
            else
                % multiply 2nd row with i and add it to 1st row
                s_data_mapped_row = [1,1i]*s_data_mapped_vec;
            end      
            
            % s_data_mapped_per_symbol
            s_data_mapped_per_symbol = transpose(reshape(s_data_mapped_row,N,length(s_data_mapped_row)/N));
            
            % adding the training symbols
            s_train = ones(1,N_train);
            
            % creating the 10-OFDM-symbol block
            s_10_OFDM_Bl_row = [s_train,s_data_mapped_row];
            
            %%
            
            % prepare s_10_OFDM_Bl_row for ifft
            s_10_OFDM_Bl_vector = transpose(reshape(s_10_OFDM_Bl_row,N,length(s_10_OFDM_Bl_row)/N));
            
            % calculate power per row before ifft
            P_s_10_OFDM_Bl_vector = power_calc(s_10_OFDM_Bl_vector);
                        
            % calculate ifft of s_10_OFDM_Bl_vector -- 1 row = 1 symbol
            S_10_OFDM_Bl_unsc_vec = zeros(size(s_10_OFDM_Bl_vector));
            for ind_s10OFDM_Bl_rec_unsc_vec = 1: size(s_10_OFDM_Bl_vector,1)
                S_10_OFDM_Bl_unsc_vec(ind_s10OFDM_Bl_rec_unsc_vec,:) = ifft(s_10_OFDM_Bl_vector(ind_s10OFDM_Bl_rec_unsc_vec,:),Nfft);                
            end
            
            % calculate power per row AFTER ifft
            P_S_10_OFDM_Bl_unsc_vec = power_calc(S_10_OFDM_Bl_unsc_vec);
            
            % scaling so that power is conserved
            S_10_OFDM_Bl_sc_vec = zeros(size(S_10_OFDM_Bl_unsc_vec));
            for ind_S_10_OFDM_Bl_rec_sc_vec = 1: size(S_10_OFDM_Bl_unsc_vec,1)
                S_10_OFDM_Bl_sc_vec(ind_S_10_OFDM_Bl_rec_sc_vec,:) = sqrt((P_s_10_OFDM_Bl_vector(ind_S_10_OFDM_Bl_rec_sc_vec)/P_S_10_OFDM_Bl_unsc_vec(ind_S_10_OFDM_Bl_rec_sc_vec))) * S_10_OFDM_Bl_unsc_vec(ind_S_10_OFDM_Bl_rec_sc_vec,:);
            end
            
            % verify that power AFTER ifft is conserved
            P_S_10_OFDM_Bl_sc_vec = power_calc(S_10_OFDM_Bl_sc_vec);
                        
            % verify that the power is preserved
            for ind_verify = 1:length(P_S_10_OFDM_Bl_unsc_vec)
                assert(abs(P_S_10_OFDM_Bl_sc_vec(ind_verify) - P_s_10_OFDM_Bl_vector(ind_verify)) < 10^(-6),'power conservation for first ifft transform is wrong');
            end
            
            % adding C.P per symbol (row)
            S_10_OFDM_Bl_sc_CP_vec = zeros((N_train+N_data)/N,N+Ncp);
            for ind_S_10_OFDM_Bl_sc_CP_vec=1:size(S_10_OFDM_Bl_sc_CP_vec,1)
                S_10_OFDM_Bl_sc_CP_vec(ind_S_10_OFDM_Bl_sc_CP_vec,:) = [S_10_OFDM_Bl_sc_vec(ind_S_10_OFDM_Bl_sc_CP_vec,end-Ncp+1:end) S_10_OFDM_Bl_sc_vec(ind_S_10_OFDM_Bl_sc_CP_vec,:)];
            end
                        
            % filtering the signal with the multi-path channel impulse
            % response per symbol
            S_10_OFDM_Bl_sc_CP_ch_vec = zeros(size(S_10_OFDM_Bl_sc_CP_vec));
            for ind_S_10_OFDM_Bl_sc_CP_ch_vec = 1: size(S_10_OFDM_Bl_sc_CP_vec,1)
                S_10_OFDM_Bl_sc_CP_ch_vec(ind_S_10_OFDM_Bl_sc_CP_ch_vec,:) = my_conv(S_10_OFDM_Bl_sc_CP_vec(ind_S_10_OFDM_Bl_sc_CP_ch_vec,:),h_sc);
            end
            
            S_10_OFDM_Bl_sc_CP_ch_noisy_vec = zeros(size(S_10_OFDM_Bl_sc_CP_ch_vec));


            % ORIGINAL IDEA FOR Eb/No calculation
            % initialize Eb vector per symbol
            for ind_noise = 1: size(S_10_OFDM_Bl_sc_CP_ch_vec,1)
                
                % calculate Eb per symbol
                Eb_S_10_OFDM_Bl_sc_CP_ch_per_row = power_calc(S_10_OFDM_Bl_sc_CP_ch_vec(ind_noise,:))/a;        

                % calculate noise variance
                sigma_sq_per_row = Eb_S_10_OFDM_Bl_sc_CP_ch_per_row / Eb_N0_dec(ind_Eb_N0_dB);
            
                % create complex noise
                if strcmp(mod_schemes(ind_m_s),'BPSK')                    
                    
                    n = sqrt(sigma_sq_per_row/2)*randn(1,size(S_10_OFDM_Bl_sc_CP_ch_vec,2));
                    %n = zeros(1,size(S_10_OFDM_Bl_sc_CP_ch_vec,2));
                else
                    n = sqrt(sigma_sq_per_row/2)*randn(1,size(S_10_OFDM_Bl_sc_CP_ch_vec,2)) + 1i*sqrt(sigma_sq_per_row/2)*randn(1,size(S_10_OFDM_Bl_sc_CP_ch_vec,2));
                end
                % adding noise
                S_10_OFDM_Bl_sc_CP_ch_noisy_vec(ind_noise,:) = S_10_OFDM_Bl_sc_CP_ch_vec(ind_noise,:) + n;
                
            end
 
%             %% ALTERNATIVE IDEA FOR Eb/No calculation
%             % initialize Eb vector per symbol
%             for ind_noise = 1: size(S_10_OFDM_Bl_sc_CP_ch_vec,1)
% 
%                 % calculate Eb per symbol
%                 Eb_S_10_OFDM_Bl_sc_CP_ch_per_row = power_calc(s_10_OFDM_Bl_vector(ind_noise,:))/a;
% 
%                 % calculate noise variance
%                 sigma_sq_per_row = Eb_S_10_OFDM_Bl_sc_CP_ch_per_row / Eb_N0_dec(ind_Eb_N0_dB);
% 
%                 %create complex noise
%                 if strcmp(mod_schemes(ind_m_s),'BPSK')
%                     n = sqrt(sigma_sq_per_row/2)*randn(1,size(S_10_OFDM_Bl_sc_CP_ch_vec,2));
%                 else
%                     n = sqrt(sigma_sq_per_row/2)*randn(1,size(S_10_OFDM_Bl_sc_CP_ch_vec,2)) + 1i*sqrt(sigma_sq_per_row/2)*randn(1,size(S_10_OFDM_Bl_sc_CP_ch_vec,2));
%                 end
%                 % adding noise
%                 S_10_OFDM_Bl_sc_CP_ch_noisy_vec(ind_noise,:) = S_10_OFDM_Bl_sc_CP_ch_vec(ind_noise,:) + n;
% 
%             end
             
             
             
             
             % remove CP
             S_10_OFDM_Bl_rec_no_CP_vec = zeros(size(S_10_OFDM_Bl_sc_vec));
             for ind_S_10_OFDM_Bl_rec_no_CP_vec=1:size(S_10_OFDM_Bl_rec_no_CP_vec,1)
                 S_10_OFDM_Bl_rec_no_CP_vec(ind_S_10_OFDM_Bl_rec_no_CP_vec,:) = S_10_OFDM_Bl_sc_CP_ch_noisy_vec(ind_S_10_OFDM_Bl_rec_no_CP_vec,Ncp+1:end);
             end
            
            % calculate power before applying fft
            P_S_10_OFDM_Bl_rec_no_CP_vector = power_calc(S_10_OFDM_Bl_rec_no_CP_vec);
            
            % apply fft --> 1 row = 1 symbol
            s_10_OFDM_Bl_rec_unsc_vec = zeros(size(S_10_OFDM_Bl_rec_no_CP_vec));            
            for ind_s10OFDM_Bl_rec_unsc_vec = 1: size(s_10_OFDM_Bl_rec_unsc_vec,1)
                s_10_OFDM_Bl_rec_unsc_vec(ind_s10OFDM_Bl_rec_unsc_vec,:) = fft(S_10_OFDM_Bl_rec_no_CP_vec(ind_s10OFDM_Bl_rec_unsc_vec,:),Nfft);                
            end     
            
            % calculate power of unscaled after applying fft
            P_s_10_OFDM_Bl_rec_unsc_vec = power_calc(s_10_OFDM_Bl_rec_unsc_vec);
            
            % scaling result of fft
            s_10_OFDM_Bl_rec_sc_vec = zeros(size(S_10_OFDM_Bl_unsc_vec));
            for ind_s_10_OFDM_Bl_rec_sc_vec = 1: size(S_10_OFDM_Bl_unsc_vec,1)
                s_10_OFDM_Bl_rec_sc_vec = sqrt((P_S_10_OFDM_Bl_rec_no_CP_vector(ind_s_10_OFDM_Bl_rec_sc_vec)/P_s_10_OFDM_Bl_rec_unsc_vec(ind_s_10_OFDM_Bl_rec_sc_vec)))*s_10_OFDM_Bl_rec_unsc_vec;
            end
                
            % calculate power of scaled after applying fft
            P_s_10_OFDM_Bl_rec_sc_vec = power_calc(s_10_OFDM_Bl_rec_sc_vec);
            
            % verify that power is preserved
            for ind_P_s_10_OFDM_Bl_rec_sc_vec = 1: size(s_10_OFDM_Bl_rec_sc_vec,1)
                assert(abs(P_S_10_OFDM_Bl_rec_no_CP_vector(ind_P_s_10_OFDM_Bl_rec_sc_vec) - P_s_10_OFDM_Bl_rec_sc_vec(ind_P_s_10_OFDM_Bl_rec_sc_vec)) < 10^(-6),'power conservation for first fft transform is wrong');
            end
            
            % use the training symbols to get the H estimates

            output_train_symb = s_10_OFDM_Bl_rec_sc_vec(1:(N_train/N),:);
            H_est_train = mean(output_train_symb,1);
            
            % equalization
            % 1. initialization
            s_data_rec_eq_vec = s_10_OFDM_Bl_rec_sc_vec(N_train/N +1:end,:);
            
            % 2. perform equalization
            for ind_s_data_rec_eq_vec = 1: size(s_data_rec_eq_vec,1)
                for ind_j = 1:size(s_data_rec_eq_vec,2)
                    s_data_rec_eq_vec(ind_s_data_rec_eq_vec,ind_j) = conj(H_est_train(ind_j))/abs(H_est_train(ind_j))^2 * s_10_OFDM_Bl_rec_sc_vec(ind_s_data_rec_eq_vec+N_train/N,ind_j);
                end
            end
            
            % convert vector to 1 row
            s_data_rec_eq_row = reshape(s_data_rec_eq_vec.',1,[]);
            
            % prepare to demap symbols -- convert the complex vector to real vector
            if strcmp(mod_schemes{ind_m_s},'BPSK')
                % in BPSK case nothing has to be done
                s_data_eq_vec = s_data_rec_eq_row;
            else
                % convert to real vector with 2 rows
                s_data_eq_vec = [real(s_data_rec_eq_row); imag(s_data_rec_eq_row)];
            end
            
            % demap
            x_data_demapped_coded = demapping(s_data_eq_vec,mod_schemes(ind_m_s));
            
            % apply the de-interleaver
            
            for ind_deint = 1: number_of_interleaving

                % writting into the de-interleaver
                de_inteleaver = reshape(x_data_demapped_coded((1+(ind_deint-1)*N_block*M_block:ind_deint*N_block*M_block)),[N_block,M_block]);

                % reading from the de-interleaver
                x_data_deinterleaved(1+(ind_deint-1)*N_block*M_block:ind_deint*N_block*M_block) = reshape(de_inteleaver',[1,N_block*M_block]);

            end
            
            % apply the viterbi decoder:
            
            % 1. Generate the trellis from the polynomials
            tr = poly2trellis( m+1,[g_octal_matr(1) g_octal_matr(2)]);

            % set the traceback length. If the code is 1/2, a typical value for it
            % is about 5 times m
            tracebacklength = 5*m;

            % call the vitdec function.
            x_data_decoded = vitdec(x_data_deinterleaved,tr,tracebacklength,'term','hard');

            % Attention: the output of vitdec is the decoded block of bits
            % INCLUDING the trailing bits

            % drop the m trailing bits
            x_data_estimated_info = x_data_decoded(1:length(x_in_cod));
            
            % count bit error rate in the current block
            [ber_coded_interl(ind_N_sim),errors]  = calc_ber_err(x_data_estimated_info,x_in_cod);           
            
            
      
            %% coded (no interleaver) version
            
            x_coded_NI = x_encoded_bits;
            
            % mapper
            s_data_mapped_NI_vec = mapping(x_coded_NI,mod_schemes(ind_m_s));
            
            % prepare to write symbols as complex numbers
            if strcmp(mod_schemes{ind_m_s},'BPSK')
                % in BPSK case nothing has to be done
                s_data_mapped_NI_row = s_data_mapped_NI_vec;
            else
                % multiply 2nd row with i and add it to 1st row
                s_data_mapped_NI_row = [1,1i]*s_data_mapped_NI_vec;
            end      
                     
            % adding the training symbols
            s_train_NI = ones(1,N_train);
            
            % creating the 10-OFDM-symbol block
            s_10_OFDM_Bl_NI_row = [s_train_NI,s_data_mapped_NI_row];
            
            %%
            
            % prepare s_10_OFDM_Bl_row for ifft
            s_10_OFDM_Bl_NI_vector = transpose(reshape(s_10_OFDM_Bl_NI_row,N,length(s_10_OFDM_Bl_NI_row)/N));
            
            % calculate power per row before ifft
            P_s_10_OFDM_Bl_NI_vector = power_calc(s_10_OFDM_Bl_NI_vector);
                        
            % calculate ifft of s_10_OFDM_Bl_vector -- 1 row = 1 symbol
            S_10_OFDM_Bl_unsc_NI_vec = zeros(size(s_10_OFDM_Bl_NI_vector));
            for ind_s10OFDM_Bl_rec_unsc_NI_vec = 1: size(s_10_OFDM_Bl_NI_vector,1)
                S_10_OFDM_Bl_unsc_NI_vec(ind_s10OFDM_Bl_rec_unsc_NI_vec,:) = ifft(s_10_OFDM_Bl_NI_vector(ind_s10OFDM_Bl_rec_unsc_NI_vec,:),Nfft);                
            end
            
            % calculate power per row AFTER ifft
            P_S_10_OFDM_Bl_unsc_NI_vec = power_calc(S_10_OFDM_Bl_unsc_NI_vec);
            
            % scaling so that power is conserved
            S_10_OFDM_Bl_sc_NI_vec = zeros(size(S_10_OFDM_Bl_unsc_NI_vec));
            for ind_S_10_OFDM_Bl_rec_sc_NI_vec = 1: size(S_10_OFDM_Bl_unsc_NI_vec,1)
                S_10_OFDM_Bl_sc_NI_vec(ind_S_10_OFDM_Bl_rec_sc_NI_vec,:) = sqrt((P_s_10_OFDM_Bl_NI_vector(ind_S_10_OFDM_Bl_rec_sc_NI_vec)/P_S_10_OFDM_Bl_unsc_NI_vec(ind_S_10_OFDM_Bl_rec_sc_NI_vec))) * S_10_OFDM_Bl_unsc_NI_vec(ind_S_10_OFDM_Bl_rec_sc_NI_vec,:);
            end
            
            % verify that power AFTER ifft is conserved
            P_S_10_OFDM_Bl_sc_NI_vec = power_calc(S_10_OFDM_Bl_sc_NI_vec);
                        
            % verify that the power is preserved
            for ind_verify_NI = 1:length(P_s_10_OFDM_Bl_NI_vector)
                assert(abs(P_S_10_OFDM_Bl_sc_NI_vec(ind_verify_NI) - P_s_10_OFDM_Bl_NI_vector(ind_verify_NI)) < 10^(-6),'power conservation for first ifft transform is wrong');
            end
            
            % adding C.P per symbol (row)
            S_10_OFDM_Bl_sc_CP_NI_vec = zeros((N_train+N_data)/N,N+Ncp);
            for ind_S_10_OFDM_Bl_sc_CP_NI_vec=1:size(S_10_OFDM_Bl_sc_CP_NI_vec,1)
                S_10_OFDM_Bl_sc_CP_NI_vec(ind_S_10_OFDM_Bl_sc_CP_NI_vec,:) = [S_10_OFDM_Bl_sc_NI_vec(ind_S_10_OFDM_Bl_sc_CP_NI_vec,end-Ncp+1:end) S_10_OFDM_Bl_sc_NI_vec(ind_S_10_OFDM_Bl_sc_CP_NI_vec,:)];
            end
                        
            % filtering the signal with the multi-path channel impulse
            % response per symbol
            S_10_OFDM_Bl_sc_CP_ch_NI_vec = zeros(size(S_10_OFDM_Bl_sc_CP_NI_vec));
            for ind_S_10_OFDM_Bl_sc_CP_ch_NI_vec = 1: size(S_10_OFDM_Bl_sc_CP_NI_vec,1)
                S_10_OFDM_Bl_sc_CP_ch_NI_vec(ind_S_10_OFDM_Bl_sc_CP_ch_NI_vec,:) = my_conv(S_10_OFDM_Bl_sc_CP_NI_vec(ind_S_10_OFDM_Bl_sc_CP_ch_NI_vec,:),h_sc);
            end
            
            S_10_OFDM_Bl_sc_CP_ch_noisy_NI_vec = zeros(size(S_10_OFDM_Bl_sc_CP_ch_NI_vec));

            % ORIGINAL IDEA FOR Eb/No calculation
            % initialize Eb vector per symbol
            for ind_noise_NI = 1: size(S_10_OFDM_Bl_sc_CP_ch_NI_vec,1)
                
                % calculate Eb per symbol
                Eb_S_10_OFDM_Bl_sc_CP_ch_per_row_NI = power_calc(S_10_OFDM_Bl_sc_CP_ch_NI_vec(ind_noise_NI,:))/a;        

                % calculate noise variance
                sigma_sq_per_row_NI = Eb_S_10_OFDM_Bl_sc_CP_ch_per_row_NI / Eb_N0_dec(ind_Eb_N0_dB);
            
                % create complex noise
                if strcmp(mod_schemes(ind_m_s),'BPSK')                    
                    
                    n_NI = sqrt(sigma_sq_per_row_NI/2)*randn(1,size(S_10_OFDM_Bl_sc_CP_ch_NI_vec,2));
                    %n_NI = zeros(1,size(S_10_OFDM_Bl_sc_CP_ch_NI_vec,2));
                else
                    n_NI = sqrt(sigma_sq_per_row_NI/2)*randn(1,size(S_10_OFDM_Bl_sc_CP_ch_NI_vec,2)) + 1i*sqrt(sigma_sq_per_row_NI/2)*randn(1,size(S_10_OFDM_Bl_sc_CP_ch_NI_vec,2));
                    %n_NI = zeros(1,size(S_10_OFDM_Bl_sc_CP_ch_NI_vec,2));
                end
                % adding noise
                S_10_OFDM_Bl_sc_CP_ch_noisy_NI_vec(ind_noise_NI,:) = S_10_OFDM_Bl_sc_CP_ch_NI_vec(ind_noise_NI,:) + n_NI;
                
            end
 
%             %% ALTERNATIVE IDEA FOR Eb/No calculation
%             % initialize Eb vector per symbol
%             for ind_noise = 1: size(S_10_OFDM_Bl_sc_CP_ch_vec,1)
% 
%                 % calculate Eb per symbol
%                 Eb_S_10_OFDM_Bl_sc_CP_ch_per_row = power_calc(s_10_OFDM_Bl_vector(ind_noise,:))/a;
% 
%                 % calculate noise variance
%                 sigma_sq_per_row = Eb_S_10_OFDM_Bl_sc_CP_ch_per_row / Eb_N0_dec(ind_Eb_N0_dB);
% 
%                 %create complex noise
%                 if strcmp(mod_schemes(ind_m_s),'BPSK')
%                     n = sqrt(sigma_sq_per_row/2)*randn(1,size(S_10_OFDM_Bl_sc_CP_ch_vec,2));
%                 else
%                     n = sqrt(sigma_sq_per_row/2)*randn(1,size(S_10_OFDM_Bl_sc_CP_ch_vec,2)) + 1i*sqrt(sigma_sq_per_row/2)*randn(1,size(S_10_OFDM_Bl_sc_CP_ch_vec,2));
%                 end
%                 % adding noise
%                 S_10_OFDM_Bl_sc_CP_ch_noisy_vec(ind_noise,:) = S_10_OFDM_Bl_sc_CP_ch_vec(ind_noise,:) + n;
% 
%             end
             
             
             % remove CP
             S_10_OFDM_Bl_rec_no_CP_NI_vec = zeros(size(S_10_OFDM_Bl_sc_NI_vec));
             for ind_S_10_OFDM_Bl_rec_no_CP_NI_vec=1:size(S_10_OFDM_Bl_rec_no_CP_NI_vec,1)
                 S_10_OFDM_Bl_rec_no_CP_NI_vec(ind_S_10_OFDM_Bl_rec_no_CP_NI_vec,:) = S_10_OFDM_Bl_sc_CP_ch_noisy_NI_vec(ind_S_10_OFDM_Bl_rec_no_CP_NI_vec,Ncp+1:end);
             end
            
            % calculate power before applying fft
            P_S_10_OFDM_Bl_rec_no_CP_NI_vector = power_calc(S_10_OFDM_Bl_rec_no_CP_NI_vec);
            
            % apply fft --> 1 row = 1 symbol
            s_10_OFDM_Bl_rec_unsc_NI_vec = zeros(size(S_10_OFDM_Bl_rec_no_CP_NI_vec));            
            for ind_s10OFDM_Bl_rec_unsc_NI_vec = 1: size(s_10_OFDM_Bl_rec_unsc_NI_vec,1)
                s_10_OFDM_Bl_rec_unsc_NI_vec(ind_s10OFDM_Bl_rec_unsc_NI_vec,:) = fft(S_10_OFDM_Bl_rec_no_CP_NI_vec(ind_s10OFDM_Bl_rec_unsc_NI_vec,:),Nfft);                
            end     
            
            % calculate power of unscaled after applying fft
            P_s_10_OFDM_Bl_rec_unsc_NI_vec = power_calc(s_10_OFDM_Bl_rec_unsc_NI_vec);
            
            % scaling result of fft
            s_10_OFDM_Bl_rec_sc_NI_vec = zeros(size(s_10_OFDM_Bl_rec_unsc_NI_vec));
            for ind_s_10_OFDM_Bl_rec_sc_NI_vec = 1: size(s_10_OFDM_Bl_rec_unsc_NI_vec,1)
                s_10_OFDM_Bl_rec_sc_NI_vec = sqrt((P_S_10_OFDM_Bl_rec_no_CP_NI_vector(ind_s_10_OFDM_Bl_rec_sc_NI_vec)/P_s_10_OFDM_Bl_rec_unsc_NI_vec(ind_s_10_OFDM_Bl_rec_sc_NI_vec)))*s_10_OFDM_Bl_rec_unsc_NI_vec;
            end
                
            % calculate power of scaled after applying fft
            P_s_10_OFDM_Bl_rec_sc_NI_vec = power_calc(s_10_OFDM_Bl_rec_sc_NI_vec);
            
            % verify that power is preserved
            for ind_P_s_10_OFDM_Bl_rec_sc_NI_vec = 1: size(s_10_OFDM_Bl_rec_sc_NI_vec,1)
                assert(abs(P_S_10_OFDM_Bl_rec_no_CP_NI_vector(ind_P_s_10_OFDM_Bl_rec_sc_NI_vec) - P_s_10_OFDM_Bl_rec_sc_NI_vec(ind_P_s_10_OFDM_Bl_rec_sc_NI_vec)) < 10^(-6),'power conservation for first fft transform is wrong');
            end
            
            % use the training symbols to get the H estimates

            output_train_symb_NI = s_10_OFDM_Bl_rec_sc_NI_vec(1:(N_train/N),:);
            H_est_train_NI = mean(output_train_symb_NI,1);
            
            % equalization
            % 1. initialization
            s_data_rec_eq_NI_vec = s_10_OFDM_Bl_rec_sc_NI_vec(N_train/N +1:end,:);
            
            % 2. perform equalization
            for ind_s_data_rec_eq_NI_vec = 1: size(s_data_rec_eq_NI_vec,1)
                for ind_j_NI = 1:size(s_data_rec_eq_NI_vec,2)
                    s_data_rec_eq_NI_vec(ind_s_data_rec_eq_NI_vec,ind_j_NI) = conj(H_est_train_NI(ind_j_NI))/abs(H_est_train_NI(ind_j_NI))^2 * s_10_OFDM_Bl_rec_sc_NI_vec(ind_s_data_rec_eq_NI_vec+N_train/N,ind_j_NI);
                end
            end
            
            % convert vector to 1 row
            s_data_rec_eq_NI_row = reshape(s_data_rec_eq_NI_vec.',1,[]);
            
            % prepare to demap symbols -- convert the complex vector to real vector
            if strcmp(mod_schemes{ind_m_s},'BPSK')
                % in BPSK case nothing has to be done
                s_data_eq_NI_vec = s_data_rec_eq_NI_row;
            else
                % convert to real vector with 2 rows
                s_data_eq_NI_vec = [real(s_data_rec_eq_NI_row); imag(s_data_rec_eq_NI_row)];
            end
            
            % demap
            x_data_demapped_coded_NI = demapping(s_data_eq_NI_vec,mod_schemes(ind_m_s));
            
            % apply the viterbi decoder:
            
            % 1. Generate the trellis from the polynomials
            tr = poly2trellis( m+1,[g_octal_matr(1) g_octal_matr(2)]);
            
            % set the traceback length. If the code is 1/2, a typical value for it
            % is about 5 times m
            tracebacklength = 5*m;
            
            % call the vitdec function.
            x_data_decoded_NI = vitdec(x_data_demapped_coded_NI,tr,tracebacklength,'term','hard');
            
            % Attention: the output of vitdec is the decoded block of bits
            % INCLUDING the trailing bits
            
            % drop the m trailing bits
            x_data_estimated_info_NI = x_data_decoded_NI(1:length(x_in_cod));
            
            % count bit error rate in the current block
            [ber_coded_NI(ind_N_sim),errors]  = calc_ber_err(x_data_estimated_info_NI,x_in_cod);          
            
            
            
            
            %% uncoded version
            x_in_uncod = randi([0 1],1,Nb_uncoded_data);
            
            % mapper
            s_data_mapped_uncoded_vec = mapping(x_in_uncod,mod_schemes(ind_m_s));
            
            % prepare to write symbols as complex numbers
            if strcmp(mod_schemes{ind_m_s},'BPSK')
                % in BPSK case nothing has to be done
                s_data_mapped_uncoded_row = s_data_mapped_uncoded_vec;
            else
                % multiply 2nd row with i and add it to 1st row
                s_data_mapped_uncoded_row = [1,1i]*s_data_mapped_uncoded_vec;
            end      
            
            % adding the training symbols
            s_train_uncoded = ones(1,N_train);
            
            % creating the 10-OFDM-symbol block
            s_10_OFDM_Bl_uncoded_row = [s_train_uncoded ,s_data_mapped_uncoded_row];
            
            % prepare s_10_OFDM_Bl_uncoded_row for ifft
            s_10_OFDM_Bl_uncoded_vector = transpose(reshape(s_10_OFDM_Bl_uncoded_row,N,length(s_10_OFDM_Bl_uncoded_row)/N));
            
            % calculate power per row before ifft
            P_s_10_OFDM_Bl_uncoded_vector = power_calc(s_10_OFDM_Bl_uncoded_vector);
                        
            % calculate ifft of s_10_OFDM_Bl_uncoded_vector -- 1 row = 1 symbol
            S_10_OFDM_Bl_unsc_uncoded_vec = zeros(size(s_10_OFDM_Bl_uncoded_vector));
            for ind_s10OFDM_Bl_rec_unsc_uncoded_vec = 1: size(s_10_OFDM_Bl_uncoded_vector,1)
                S_10_OFDM_Bl_unsc_uncoded_vec(ind_s10OFDM_Bl_rec_unsc_uncoded_vec ,:) = ifft(s_10_OFDM_Bl_uncoded_vector(ind_s10OFDM_Bl_rec_unsc_uncoded_vec ,:),Nfft);                
            end
            
            % calculate power per row AFTER ifft
            P_S_10_OFDM_Bl_unsc_uncoded_vec = power_calc(S_10_OFDM_Bl_unsc_uncoded_vec);
            
            % scaling so that power is conserved
            S_10_OFDM_Bl_sc_uncoded_vec = zeros(size(S_10_OFDM_Bl_unsc_uncoded_vec));
            for ind_S_10_OFDM_Bl_rec_sc_uncoded_vec = 1: size(S_10_OFDM_Bl_unsc_uncoded_vec,1)
                S_10_OFDM_Bl_sc_uncoded_vec(ind_S_10_OFDM_Bl_rec_sc_uncoded_vec,:) = sqrt((P_s_10_OFDM_Bl_uncoded_vector(ind_S_10_OFDM_Bl_rec_sc_uncoded_vec)/P_S_10_OFDM_Bl_unsc_uncoded_vec(ind_S_10_OFDM_Bl_rec_sc_uncoded_vec))) * S_10_OFDM_Bl_unsc_uncoded_vec(ind_S_10_OFDM_Bl_rec_sc_uncoded_vec,:);
            end
            
            % verify that power AFTER ifft is conserved
            P_S_10_OFDM_Bl_sc_uncoded_vec = power_calc(S_10_OFDM_Bl_sc_uncoded_vec);
                        
            % verify that the power is preserved
            for ind_verify = 1:length(P_S_10_OFDM_Bl_unsc_uncoded_vec)
                assert(abs(P_S_10_OFDM_Bl_sc_uncoded_vec(ind_verify) - P_s_10_OFDM_Bl_uncoded_vector(ind_verify)) < 10^(-6),'power conservation for first ifft transform is wrong');
            end
            
            % adding C.P per symbol (row)
            S_10_OFDM_Bl_sc_CP_uncoded_vec = zeros((N_train+N_data)/N,N+Ncp);
            for ind_S_10_OFDM_Bl_sc_CP_uncoded_vec=1:size(S_10_OFDM_Bl_sc_CP_uncoded_vec,1)
                S_10_OFDM_Bl_sc_CP_uncoded_vec(ind_S_10_OFDM_Bl_sc_CP_uncoded_vec,:) = [S_10_OFDM_Bl_sc_uncoded_vec(ind_S_10_OFDM_Bl_sc_CP_uncoded_vec,end-Ncp+1:end) S_10_OFDM_Bl_sc_uncoded_vec(ind_S_10_OFDM_Bl_sc_CP_uncoded_vec,:)];
            end
                        
            % filtering the signal with the multi-path channel impulse
            % response per symbol
            S_10_OFDM_Bl_sc_CP_ch_uncoded_vec = zeros(size(S_10_OFDM_Bl_sc_CP_uncoded_vec));
            for ind_S_10_OFDM_Bl_sc_CP_ch_uncoded_vec = 1: size(S_10_OFDM_Bl_sc_CP_uncoded_vec,1)
                S_10_OFDM_Bl_sc_CP_ch_uncoded_vec(ind_S_10_OFDM_Bl_sc_CP_ch_uncoded_vec,:) = my_conv(S_10_OFDM_Bl_sc_CP_uncoded_vec(ind_S_10_OFDM_Bl_sc_CP_ch_uncoded_vec,:),h_sc);
            end
            
            S_10_OFDM_Bl_sc_CP_ch_noisy_uncoded_vec = zeros(size(S_10_OFDM_Bl_sc_CP_ch_uncoded_vec));


            % ORIGINAL IDEA FOR Eb/No calculation
            % initialize Eb vector per symbol
            for ind_noise_uncoded = 1: size(S_10_OFDM_Bl_sc_CP_ch_uncoded_vec,1)
                
                % calculate Eb per symbol
                Eb_S_10_OFDM_Bl_sc_CP_ch_per_row_unncoded = power_calc(S_10_OFDM_Bl_sc_CP_ch_uncoded_vec(ind_noise_uncoded,:))/a;        

                % calculate noise variance
                sigma_sq_per_row_uncoded = Eb_S_10_OFDM_Bl_sc_CP_ch_per_row_unncoded / Eb_N0_dec(ind_Eb_N0_dB);
            
                % create complex noise
                if strcmp(mod_schemes(ind_m_s),'BPSK')
                                     
                    n = sqrt(sigma_sq_per_row_uncoded/2)*randn(1,size(S_10_OFDM_Bl_sc_CP_ch_uncoded_vec,2));
                    % n = zeros(1,size(S_10_OFDM_Bl_sc_CP_ch_uncoded_vec,2));
                else
                    n = sqrt(sigma_sq_per_row_uncoded/2)*randn(1,size(S_10_OFDM_Bl_sc_CP_ch_uncoded_vec,2)) + 1i*sqrt(sigma_sq_per_row_uncoded/2)*randn(1,size(S_10_OFDM_Bl_sc_CP_ch_uncoded_vec,2));
                    %n = zeros(1,size(S_10_OFDM_Bl_sc_CP_ch_uncoded_vec,2));
                end
                % adding noise
                S_10_OFDM_Bl_sc_CP_ch_noisy_uncoded_vec(ind_noise_uncoded,:) = S_10_OFDM_Bl_sc_CP_ch_uncoded_vec(ind_noise_uncoded,:) + n;
                
            end
           
            % remove CP
            S_10_OFDM_Bl_rec_no_CP_uncoded_vec = zeros(size(S_10_OFDM_Bl_sc_uncoded_vec));
            for ind_S_10_OFDM_Bl_rec_no_CP_uncoded_vec=1:size(S_10_OFDM_Bl_rec_no_CP_uncoded_vec,1)
                S_10_OFDM_Bl_rec_no_CP_uncoded_vec(ind_S_10_OFDM_Bl_rec_no_CP_uncoded_vec,:) = S_10_OFDM_Bl_sc_CP_ch_noisy_uncoded_vec(ind_S_10_OFDM_Bl_rec_no_CP_uncoded_vec,Ncp+1:end);
            end
            
            % calculate power before applying fft
            P_S_10_OFDM_Bl_rec_no_CP_uncoded_vector = power_calc(S_10_OFDM_Bl_rec_no_CP_uncoded_vec);
            
            % apply fft --> 1 row = 1 symbol
            s_10_OFDM_Bl_rec_unsc_uncoded_vec = zeros(size(S_10_OFDM_Bl_rec_no_CP_uncoded_vec));            
            for ind_s10OFDM_Bl_rec_unsc_uncoded_vec = 1: size(s_10_OFDM_Bl_rec_unsc_uncoded_vec,1)
                s_10_OFDM_Bl_rec_unsc_uncoded_vec(ind_s10OFDM_Bl_rec_unsc_uncoded_vec,:) = fft(S_10_OFDM_Bl_rec_no_CP_uncoded_vec(ind_s10OFDM_Bl_rec_unsc_uncoded_vec,:),Nfft);                
            end     
            
            % calculate power of unscaled after applying fft
            P_s_10_OFDM_Bl_rec_unsc_uncoded_vec = power_calc(s_10_OFDM_Bl_rec_unsc_uncoded_vec);
            
            % scaling result of fft
            s_10_OFDM_Bl_rec_sc_uncoded_vec = zeros(size(s_10_OFDM_Bl_rec_unsc_uncoded_vec));
            for ind_s_10_OFDM_Bl_rec_sc_uncoded_vec = 1: size(s_10_OFDM_Bl_rec_unsc_uncoded_vec,1)
                s_10_OFDM_Bl_rec_sc_uncoded_vec = sqrt((P_S_10_OFDM_Bl_rec_no_CP_uncoded_vector(ind_s_10_OFDM_Bl_rec_sc_uncoded_vec)/P_s_10_OFDM_Bl_rec_unsc_uncoded_vec(ind_s_10_OFDM_Bl_rec_sc_uncoded_vec)))*s_10_OFDM_Bl_rec_unsc_uncoded_vec;
            end
                
            % calculate power of scaled after applying fft
            P_s_10_OFDM_Bl_rec_sc_uncoded_vec = power_calc(s_10_OFDM_Bl_rec_sc_uncoded_vec);
            
            % verify that power is preserved
            for ind_P_s_10_OFDM_Bl_rec_sc_uncoded_vec = 1: size(s_10_OFDM_Bl_rec_sc_uncoded_vec,1)
                assert(abs(P_S_10_OFDM_Bl_rec_no_CP_uncoded_vector(ind_P_s_10_OFDM_Bl_rec_sc_uncoded_vec) - P_s_10_OFDM_Bl_rec_sc_uncoded_vec(ind_P_s_10_OFDM_Bl_rec_sc_uncoded_vec)) < 10^(-6),'power conservation for first UNCODED fft transform is wrong ');
            end
            
            % use the training symbols to get the H estimates

            output_train_symb_uncoded = s_10_OFDM_Bl_rec_sc_uncoded_vec(1:(N_train/N),:);
            H_est_train_uncoded = mean(output_train_symb_uncoded,1);
            
            % equalization
            % 1. initialization
            s_data_rec_eq_uncoded_vec = s_10_OFDM_Bl_rec_sc_uncoded_vec(N_train/N +1:end,:);
            
            % 2. perform equalization
            for ind_s_data_rec_eq_uncoded_vec = 1: size(s_data_rec_eq_uncoded_vec,1)
                for ind_j_uncoded = 1:size(s_data_rec_eq_uncoded_vec,2)
                    s_data_rec_eq_uncoded_vec(ind_s_data_rec_eq_uncoded_vec,ind_j_uncoded) = conj(H_est_train_uncoded(ind_j_uncoded))/abs(H_est_train_uncoded(ind_j_uncoded))^2 * s_10_OFDM_Bl_rec_sc_uncoded_vec(ind_s_data_rec_eq_uncoded_vec+N_train/N,ind_j_uncoded);
                end
            end
            
            % convert vector to 1 row
            s_data_rec_eq_uncoded_row = reshape(s_data_rec_eq_uncoded_vec.',1,[]);
            
            % prepare to demap symbols -- convert the complex vector to real vector
            if strcmp(mod_schemes{ind_m_s},'BPSK')
                % in BPSK case nothing has to be done
                s_data_eq_uncoded_vec = s_data_rec_eq_uncoded_row;
            else
                % convert to real vector with 2 rows
                s_data_eq_uncoded_vec = [real(s_data_rec_eq_uncoded_row); imag(s_data_rec_eq_uncoded_row)];
            end
            
            % demap
            x_data_demapped_uncoded = demapping(s_data_eq_uncoded_vec,mod_schemes(ind_m_s));
            
            % count bit error rate in the current block
            [ber_uncoded(ind_N_sim),errors]  = calc_ber_err(x_data_demapped_uncoded,x_in_uncod);  
            
            
        end
        BER_coded(ind_m_s,ind_Eb_N0_dB) = mean(ber_coded_interl);
        BER_uncoded(ind_m_s,ind_Eb_N0_dB) = mean(ber_uncoded);
        BER_coded_NI(ind_m_s,ind_Eb_N0_dB) = mean(ber_coded_NI);        
    end  
end

toc;

%% Compare coded AND interl BER for different mod schemes

figure();
semilogy(Eb_N0_dB,BER_coded(1,:), 'b--*', 'LineWidth', 3, 'MarkerSize',10);
xlabel('Eb/N0 in dB');
ylabel('BER in log scale');
xlim([min(Eb_N0_dB) max(Eb_N0_dB)]);
ylim([5*10^(-log10(8*N_sim)) 1]);
ax = gca;
set(gca,'FontSize',18);
ax.GridAlpha = 1;
ax.MinorGridAlpha = 0.55;
grid on;
hold on;
semilogy(Eb_N0_dB, BER_coded(2,:), 'r--+', 'LineWidth', 3, 'MarkerSize',10);
hold on;
semilogy(Eb_N0_dB, BER_coded(3,:), 'k--x', 'LineWidth', 3, 'MarkerSize',10);
hold on;
semilogy(Eb_N0_dB, BER_coded(4,:), 'g--o', 'LineWidth', 3, 'MarkerSize',10);
leg = legend({'BPSK','QPSK','16-QAM','64-QAM'},'Location','southwest');
set(leg,'FontSize',18);
title('Comparison for BER coded and interleaved across different mod schemes');
hold off;

figure();
semilogy(Eb_N0_dB,BER_uncoded(1,:), 'b--*', 'LineWidth', 3, 'MarkerSize',10);
xlabel('Eb/N0 in dB');
ylabel('BER in log scale');
xlim([min(Eb_N0_dB) max(Eb_N0_dB)]);
ylim([5*10^(-log10(8*N_sim)) 1]);
ax = gca;
set(gca,'FontSize',18);
ax.GridAlpha = 1;
ax.MinorGridAlpha = 0.55;
grid on;
hold on;
semilogy(Eb_N0_dB, BER_uncoded(2,:), 'r--+', 'LineWidth', 3, 'MarkerSize',10);
hold on;
semilogy(Eb_N0_dB, BER_uncoded(3,:), 'k--x', 'LineWidth', 3, 'MarkerSize',10);
hold on;
semilogy(Eb_N0_dB, BER_uncoded(4,:), 'g--o', 'LineWidth', 3, 'MarkerSize',10);
leg = legend({'BPSK','QPSK','16-QAM','64-QAM'},'Location','southwest');
set(leg,'FontSize',18);
title('Comparison for BER uncoded across different mod schemes');
hold off;

figure();
semilogy(Eb_N0_dB,BER_coded_NI(1,:), 'b--*', 'LineWidth', 3, 'MarkerSize',10);
xlabel('Eb/N0 in dB');
ylabel('BER in log scale');
xlim([min(Eb_N0_dB) max(Eb_N0_dB)]);
ylim([5*10^(-log10(8*N_sim)) 1]);
ax = gca;
set(gca,'FontSize',18);
ax.GridAlpha = 1;
ax.MinorGridAlpha = 0.55;
grid on;
hold on;
semilogy(Eb_N0_dB, BER_coded_NI(2,:), 'r--+', 'LineWidth', 3, 'MarkerSize',10);
hold on;
semilogy(Eb_N0_dB, BER_coded_NI(3,:), 'k--x', 'LineWidth', 3, 'MarkerSize',10);
hold on;
semilogy(Eb_N0_dB, BER_coded_NI(4,:), 'g--o', 'LineWidth', 3, 'MarkerSize',10);
leg = legend({'BPSK','QPSK','16-QAM','64-QAM'},'Location','southwest');
set(leg,'FontSize',18);
title('Comparison for BER coded but NO interleaving across different mod schemes');
hold off;

figure();
semilogy(Eb_N0_dB,BER_coded(1,:), 'b--*', 'LineWidth', 3, 'MarkerSize',10);
xlabel('Eb/N0 in dB');
ylabel('BER in log scale');
xlim([min(Eb_N0_dB) max(Eb_N0_dB)]);
ylim([5*10^(-log10(8*N_sim)) 1]);
ax = gca;
set(gca,'FontSize',18);
ax.GridAlpha = 1;
ax.MinorGridAlpha = 0.55;
grid on;
hold on;
semilogy(Eb_N0_dB, BER_coded_NI(1,:), 'r--+', 'LineWidth', 3, 'MarkerSize',10);
hold on;
semilogy(Eb_N0_dB, BER_uncoded(1,:), 'k--x', 'LineWidth', 3, 'MarkerSize',10);
leg = legend({'coded & inter','coded & NO inter','uncoded'},'Location','southwest');
set(leg,'FontSize',18);
title('Comparison for coded&inter, coded&NO interv and uncoded for BPSK');
hold off;

figure();
semilogy(Eb_N0_dB,BER_coded(2,:), 'b--*', 'LineWidth', 3, 'MarkerSize',10);
xlabel('Eb/N0 in dB');
ylabel('BER in log scale');
xlim([min(Eb_N0_dB) max(Eb_N0_dB)]);
ylim([5*10^(-log10(8*N_sim)) 1]);
ax = gca;
set(gca,'FontSize',18);
ax.GridAlpha = 1;
ax.MinorGridAlpha = 0.55;
grid on;
hold on;
semilogy(Eb_N0_dB, BER_coded_NI(2,:), 'r--+', 'LineWidth', 3, 'MarkerSize',10);
hold on;
semilogy(Eb_N0_dB, BER_uncoded(2,:), 'k--x', 'LineWidth', 3, 'MarkerSize',10);
leg = legend({'coded & inter','coded & NO inter','uncoded'},'Location','southwest');
set(leg,'FontSize',18);
title('Comparison for coded&inter, coded&NO interv and uncoded for QPSK');
hold off;

figure();
semilogy(Eb_N0_dB,BER_coded(3,:), 'b--*', 'LineWidth', 3, 'MarkerSize',10);
xlabel('Eb/N0 in dB');
ylabel('BER in log scale');
xlim([min(Eb_N0_dB) max(Eb_N0_dB)]);
ylim([5*10^(-log10(8*N_sim)) 1]);
ax = gca;
set(gca,'FontSize',18);
ax.GridAlpha = 1;
ax.MinorGridAlpha = 0.55;
grid on;
hold on;
semilogy(Eb_N0_dB, BER_coded_NI(3,:), 'r--+', 'LineWidth', 3, 'MarkerSize',10);
hold on;
semilogy(Eb_N0_dB, BER_uncoded(3,:), 'k--x', 'LineWidth', 3, 'MarkerSize',10);
leg = legend({'coded & inter','coded & NO inter','uncoded'},'Location','southwest');
set(leg,'FontSize',18);
title('Comparison for coded&inter, coded&NO interv and uncoded for 16-QAM');
hold off;

figure();
semilogy(Eb_N0_dB,BER_coded(4,:), 'b--*', 'LineWidth', 3, 'MarkerSize',10);
xlabel('Eb/N0 in dB');
ylabel('BER in log scale');
xlim([min(Eb_N0_dB) max(Eb_N0_dB)]);
ylim([5*10^(-log10(8*N_sim)) 1]);
ax = gca;
set(gca,'FontSize',18);
ax.GridAlpha = 1;
ax.MinorGridAlpha = 0.55;
grid on;
hold on;
semilogy(Eb_N0_dB, BER_coded_NI(4,:), 'r--+', 'LineWidth', 3, 'MarkerSize',10);
hold on;
semilogy(Eb_N0_dB, BER_uncoded(4,:), 'k--x', 'LineWidth', 3, 'MarkerSize',10);
leg = legend({'coded & inter','coded & NO inter','uncoded'},'Location','southwest');
set(leg,'FontSize',18);
title('Comparison for coded&inter, coded&NO interv and uncoded for 64-QAM');
hold off;
