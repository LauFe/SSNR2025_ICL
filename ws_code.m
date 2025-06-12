%% Laura Ferrante - Natural BionicS workshop for SSNR2025 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data: high density surface EMG (sEMG) and high density intramuscular (iEMG)
load('sEMG.mat');load('iEMG.mat');
% Frequency content of sEMG and iEMG. Estimate PSD using Welch's overlapped segment averaging estimatior.
figure
subplot(2,1,1)
[Pxx,F] = pwelch(sEMG.data,hamming(sEMG.fs),(sEMG.fs/2)/2,sEMG.fs*2+1,sEMG.fs,'power'); %sEMG
plot(F,Pxx)
xlim([0 500])
xlabel('Frequency (Hz)')
ylabel('PDS')
title('HD surface EMG data')
subplot(2,1,2)
[Pxx,F] = pwelch(iEMG.data(:,1:40),hamming(iEMG.fs),(iEMG.fs/2)/2,iEMG.fs*2+1,iEMG.fs,'power'); % iEMG
plot(F,Pxx)
xlim([0 2000])
xlabel('Frequency (Hz)')
ylabel('PDS')
title('HD intramuscular EMG data')

% filter HD-sEMG signals
sEMG.fnyq = sEMG.fs / 2;
sEMG.fc_high = 10;
sEMG.fc_low = 500;
fnotch = 50;
[b,a] = butter(4,[sEMG.fc_high sEMG.fc_low]/sEMG.fnyq,'bandpass');
sEMG.data_filt = filtfilt(b,a,sEMG.data); % data (samples x channels)
[b,a] = butter(4,[(fnotch-2)/sEMG.fnyq (fnotch+2)/sEMG.fnyq],'stop');
sEMG.data_filt = filtfilt(b,a,sEMG.data_filt); % data (samples x channels)

keyboard;
% Run EMGLAB motor unit decomposition
global EMGSIGNAL
EMGSIGNAL = struct('data',iEMG.data,'rate',iEMG.fs);
emglab;
% load annotation .eaf file

%% 
load('decomp.mat');
% MU 1 as example
n_MU = size(decomp.PinkyFlexion.MU,2);
nCh = size(iEMG.data,2);
% Raster plot
h0 = figure;
for ii = 1:n_MU
    hold on
    [~,h0] = plotRaster_single(h0,decomp.PinkyFlexion.MU(1,ii),iEMG.fs,'r','-',ii);
end
ylabel('MU#')
xlabel('Time [s]')

% Spike triggered average
iEMG.sta = cell(40,n_MU);
for ii = 1:40
    [STA_ch] = extr_STA(iEMG.data(:,ii)',decomp.PinkyFlexion.MU,iEMG.fs,0);
    iEMG.sta(ii,:) = STA_ch;
end
nPoints = size(cell2mat(iEMG.sta(:,1)),2);
X = cell2mat(iEMG.sta(:,1));
h = figure;
set(h,'color','w');
set(findall(h, '-property', 'fontsize'), 'fontsize', 18)
set(h,'units','points','position',[10,10,400,400]) 
            
[xx,yy] = meshgrid(1:nPoints,1:nCh);
pfig = pcolor(xx/iEMG.fs,yy,X);
xlabel('Time [s]');
ylabel('Ch#');
colormap('jet')
grid off
colorbar('eastoutside');
set(pfig, 'EdgeColor', 'none');

%% NMF - estimating the dimensionality latent space 
firing = zeros(n_MU,size(iEMG.data,1)); 
nmfdata = [];
figure
for jj=1:n_MU
    firing(jj,decomp.PinkyFlexion.MU{1,jj})= 1; 
    win = hann(round(0.4 * iEMG.fs));
    hold on
    plot([0:length(firing)-1]/iEMG.fs,fftfilt(win,firing'),'linewidth',1); 
    smoothed_spikes = fftfilt(win,firing');
end
nmf_spike_weight = max(smoothed_spikes',[],2);
nmf_spike_weight(nmf_spike_weight==0)=1; 
nmfdata = smoothed_spikes./repmat(nmf_spike_weight',size(smoothed_spikes,1),1);
nmfdata(nmfdata<0) = 0;
[wcoeff_1,~,latent_1,~,explained_1] = pca(nmfdata);