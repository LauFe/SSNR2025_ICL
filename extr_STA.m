function [out] = extr_STA(sig, pulses, fs, diff_mode)
%EXTR_STA computes the spike-triggered average of the input EMG signal
%(sig) based on the pulses of the reference motor units (MU) (pulses) 35 ms
%window centered at each action potential.
%INPUTS:
% - sig: EMG signals stored in a cell following the channel arrangement
%        {n_channel_rows, n_channel_cols}(1 x sig_length)
% - pulses: MU pulses arranged in {1 x nMU}(1 x n_pulses)
% - fs: sampling frequency [Hz]
% - diff_mode: flag to apply DEMUSE differential high pass filter
%OUTPUT:
% - out: spike-triggered average of each MU's pulses for each EMG channel:
%        {n_channel_rows, n_channel_cols, n_MU}(1 x MUAP_length)
flag_perch = 0;
if size(sig,1) > size(sig,2)
    sig = sig';
end
% Preallocate memory
out = cell(size(sig,1),length(pulses(1,:)));
win = round(0.01*fs); %round(0.0125*fs);
win2 = round(0.005*fs);
if diff_mode
    % Differential filter
    Nf = 50;
    Fpass = 450;
    Fstop = 500;
    d = designfilt('differentiatorfir','FilterOrder',Nf, ...
    'PassbandFrequency',Fpass,'StopbandFrequency',Fstop, ...
    'SampleRate',fs);
end
% Preallocate arrays
num_mu = size(pulses, 2);
num_r = size(sig, 1);

out = cell(num_r, num_mu);
pulses2 = pulses;
% For every MU
for mu = 1:size(pulses,2)
        for r = 1:size(sig,1)
           
            if flag_perch == 0
                if diff_mode 
                    aux = zeros(length(pulses{mu}), 2*win+Nf+1);
                else
                    aux = zeros(length(pulses{mu}), 2*win+1);
                end
                % For every MU's pulses
                for pulse = 1:length(pulses{mu})
                    if ~isempty(sig(r,:))
                        % Check that window falls within the signal size
                        l_idx = pulses{mu}(pulse)-win;
                        h_idx = pulses{mu}(pulse)+win;

                        l_idx2 = pulses{mu}(pulse)-win2;
                        h_idx2 = pulses{mu}(pulse)+win2;

                        if diff_mode
                            % Check that window falls within the signal size
                            l_idx = pulses{mu}(pulse)-win-Nf/2;
                            h_idx = pulses{mu}(pulse)+win+Nf/2;
                        end

                        %if l_idx > 1 && h_idx < length(sig{r,c})
                        if l_idx > 1 && h_idx < length(sig(r,:))
                            %aux(pulse,:) = sig{r,c}(l_idx: h_idx);
                            aux(pulse,:) = sig(r,l_idx: h_idx);
                            get_sig = sig(r,l_idx2: h_idx2);
                            win_tmp = l_idx: h_idx;
                            tmp_index = find(abs(get_sig) >= 0.15*max(abs(get_sig)));
                            new_pulse = win_tmp(tmp_index(1));
                            pulses2{mu}(pulse) = new_pulse;
                            
                        elseif l_idx < 1 % The window is too large on the left side
                            range_idx = l_idx : h_idx;
                            % Take from the first sample up to the upper limit and use the first value to fill the window vector
                          
                            aux(pulse,range_idx > 0) = sig(r,1: h_idx);
                            aux(pulse,range_idx < 0) = ones(1,sum(range_idx < 0)).*sig(r,1);
                           
                        %elseif h_idx > length(sig{r,c}) % The window is too large on the right side
                        elseif h_idx > length(sig(r,:))
                            range_idx = l_idx : h_idx;
                            % Take from the lower limit up to the last sample and use that value to fill the window vector
                            
                            aux(pulse,range_idx <= length(sig(r,:))) = sig(r,l_idx: end);
                            aux(pulse,range_idx > length(sig(r,:))) = ones(1,sum(range_idx > length(sig(r,:)))).*sig(r,end);
                        end
                        

                    end

                    if diff_mode
                        % Time differentiation (first derivative)
                        aux(pulse,:) = filter(d, aux(pulse,:));
                    end
                end
            else
                if diff_mode 
                    aux = zeros(length(pulses{mu}), 2*win+Nf+1);
                else
                    aux = zeros(length(pulses{r,mu}), 2*win+1);
                end
                % For every MU's pulses
                for pulse = 1:length(pulses{r,mu})

                    %if ~isempty(sig{r,c})
                    if ~isempty(sig(r,:))
                        % Check that window falls within the signal size
                        l_idx = pulses{r,mu}(pulse)-win;
                        h_idx = pulses{r,mu}(pulse)+win;

                        if diff_mode
                            % Check that window falls within the signal size
                            l_idx = pulses{r,mu}(pulse)-win-Nf/2;
                            h_idx = pulses{r,mu}(pulse)+win+Nf/2;
                        end

                        %if l_idx > 1 && h_idx < length(sig{r,c})
                        if l_idx > 1 && h_idx < length(sig(r,:))
                            %aux(pulse,:) = sig{r,c}(l_idx: h_idx);
                            aux(pulse,:) = sig(r,l_idx: h_idx);

                        elseif l_idx < 1 % The window is too large on the left side
                            range_idx = l_idx : h_idx;
                            % Take from the first sample up to the upper limit and use the first value to fill the window vector
                            %aux(pulse,range_idx > 0) = sig{r,c}(1: h_idx);
                            aux(pulse,range_idx > 0) = sig(r,1: h_idx);
                            %aux(pulse,range_idx < 0) = ones(1,sum(range_idx < 0)).*sig{r,c}(1);
                            aux(pulse,range_idx < 0) = ones(1,sum(range_idx < 0)).*sig(r,1);

                        %elseif h_idx > length(sig{r,c}) % The window is too large on the right side
                            elseif h_idx > length(sig(r,:))
                            range_idx = l_idx : h_idx;
                            % Take from the lower limit up to the last sample and use that value to fill the window vector
                            %aux(pulse,range_idx <= length(sig{r,c})) = sig{r,c}(l_idx: end);
                            %aux(pulse,range_idx > length(sig{r,c})) = ones(1,sum(range_idx > length(sig{r,c}))).*sig{r,c}(1,end);
                            aux(pulse,range_idx <= length(sig(r,:))) = sig(r,l_idx: end);
                            aux(pulse,range_idx > length(sig(r,:))) = ones(1,sum(range_idx > length(sig(r,:)))).*sig(r,end);
                        end
                    end

                    if diff_mode
                        % Time differentiation (first derivative)
                        aux(pulse,:) = filter(d,aux(pulse,:));
                    end
                end
            end
            if diff_mode
                % Discard  initial filter transient and delay
                aux(:,[1:Nf/2 end-Nf/2:end]) = [];
            end
                        
            % Compute the spike triggered average
            %out{r,c,mu} = mean(aux);
            if size(aux,1) == 1
                out{r,mu} = aux;
            else
                out{r,mu} = mean(aux);
            end
        end
    
end
end
