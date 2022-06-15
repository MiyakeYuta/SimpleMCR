%% load test signal para for randome waveform generation
% for generation of designated waveform move %% disignated waveform
% parameter
load 'tsg_para2'
sig_sz = tsg_para.signal_size;
pk_num = tsg_para.pk_num;
pk_amp_max = tsg_para.peak_amp_max;
pk_sigma_rng = tsg_para.peak_sigma_range;  % full width
csv_norm_cst = tsg_para.csv_norm_const;
num_sig = tsg_para.number_signals;
noise_amp = tsg_para.noise_amp;
noise_add = tsg_para.noise_add;
%% Shift data parameter
% if you want to make some specific waveform, just input values in the
% following parameters.
rng("shuffle");
center_pos = randi([1, sig_sz], 1, pk_num);
pk_amp = rand(1, pk_num) * pk_amp_max;
pk_sigma = randi([pk_sigma_rng(1), pk_sigma_rng(2)], 1, pk_num);
% rng("shuffle");

%% Designated waveform parameter
%load 'tsg_para_save'
sig_num = tsg_para_save.sig_num;
pk_num = tsg_para_save.pk_num;
pk_amp = tsg_para_save.pk_amp;
center_pos = tsg_para_save.center_pos;
pk_sigma = tsg_para_save.pk_sigma;
num_sig = tsg_para_save.num_sig;
sig_sz = tsg_para_save.sig_size;
noise_amp = tsg_para_save.noise_amp;
noise_add = tsg_para_save.noise_add;

%% make the test signal
close all
for l = 1:num_sig
    l
    tst_pk(1:sig_sz) = 0;

    for i = 1:pk_num
        for j = 1:sig_sz
            % Lorentz distribution, pk_sigma : half of full width needed, not to
            % use the probability distribution function (normalization is
            % different)
%             tst_pk(j) = tst_pk(j) + pk_amp(i)  / (1 + ((j - center_pos(i)) / (0.5 * pk_sigma(i)))^2);
            % Gaussian distribution
            tst_pk(j) = tst_pk(j) + pk_amp(i) * exp(-(j - center_pos(i))^2 / (2 * pk_sigma(i)^2));
        end
    end
  
    tst_sig_seq(:,l) = tst_pk';
    if noise_add == "on"
        noise = noise_amp * rand(sig_sz, 1);
        tst_sig_seq(:,l) = tst_sig_seq(:,l) + noise;
    end

end

for i = 1 : num_sig
   plot(tst_sig_seq(:,i));
   hold on;
end

signal = tst_sig_seq;
tsg_para_save.signal = signal;
filename = "signal data " + num2str(sig_num) + ".mat";
save(filename, "tsg_para_save");

%% 
pure_spec = [];
%% save the signal of non-shifted
% manually saved the files 3A, 3B, 3C,..., and load the data
signal = tsg_para_save.signal;
%% combine each signal
pure_spec = [pure_spec signal];
%% make mixture signal
close all
num_comp = 3;
num_spec = 40;
conc = rand(num_comp, num_spec);
mix_spec = pure_spec * conc;

figure
plot(mix_spec(:,1:end))

%% combine different spectrum sets
mix_spec = [mix_spec_A, mix_spec_B, mix_spec_C];
mix_conc = [conc_A, conc_B, conc_C];

plot(mix_spec(:, 1:end))
