% Program utama untuk Dual Tree Discrete Wavelet Transform terhadap sinyal FPCG
% Nama  : Tuah Jihan Manullang
% prodi : S1 Teknik Telekomunikasi

% Environment 
warning off;
clear all;
close all;
clc;

%% Memilih folder untuk menyimpan
direk =  uigetdir('Choose a folder where you store the data');

if ~isequal(direk, 0)
    
    Nfiles = dir(fullfile(direk, '*.wav'));
    J = 4;
    for id = 4 
        fprintf('DT-CWT Level 4');
        for ix = 1:numel(Nfiles);
            
             % Import data ke Matlab
            fname = Nfiles(ix).name;
            dname = fullfile(direk, fname);
            [x, fs] = audioread(dname);
            
            fprintf('%d) File: %s', ix, fname);
            
            %% Input Signal
            x = x(1:512);
            
            % normalisasi data mentah agar berada pada -1 hingga +1 volt
            x = x ./ max(abs(x));
            
            % centering
            x = x - mean(x);
            
            % tambahkan noise acak N(0,1)
            % datan = wgn(length(data), 1, 0)';
            snrawgn = 5;
            datan = awgn(x, snrawgn, 'measured'); %Input Signal+Noise
            xnoise = x+datan;

%% Create the Analysis Filters for the two trees
% First-stage Filter coefficients - Tree R
Faf{1} = [0         0
   -0.0884   -0.0112
    0.0884    0.0112
    0.6959    0.0884
    0.6959    0.0884
    0.0884   -0.6959
   -0.0884    0.6959
    0.0112   -0.0884
    0.0112   -0.0884
         0         0];

% First-stage filter coefficients - Tree I     
Faf{2} = [ 0.0112  0
    0.0112         0
   -0.0884   -0.0884
    0.0884   -0.0884
    0.6959    0.6959
    0.6959   -0.6959
    0.0884    0.0884
   -0.0884    0.0884
         0    0.0112
         0   -0.0112];
 
 %Second and subsequent stages filter coefficients - Tree R
 af{1} = [ 0.0352         0
         0         0
   -0.0883   -0.1143
    0.2339         0
    0.7603    0.5875
    0.5875   -0.7603
         0    0.2339
   -0.1143    0.0883
         0         0
         0   -0.0352];

%Second and Subsequent stages filter coefficients - Tree I
af{2} = [0   -0.0352
         0         0
   -0.1143    0.0883
         0    0.2339
    0.5875   -0.7603
    0.7603    0.5875
    0.2339         0
   -0.0883   -0.1143
         0         0
    0.0352         0];

%% Hitung Global Threshold
% Menghitung Noise Level Estimation
param = 'h';
[nLevel, xs] = NoiselevelEstimation(xnoise);

%% Menghitung Threshold
NoiseVar = nLevel.^2;
lenData = length(xnoise);
Threshold = sqrt(2.*NoiseVar*log10(lenData));

% Proses Denoising secara global
wt = dddtree('cplxdt',xnoise,J,Faf,af);

% Ekstrak Koefisien Detail dari Tree R dan Tree I, pada setiap Level
outputindices = {[1 1]; [1 2]; [2 1]; [2 2]; [3 1]; [3 2]; [4 1]; [4 2]};
cD = dddtreecfs('e',wt,'ind',outputindices);
xapp = dddtreecfs('e',wt,'lowpass');

% proses denoising
Cnew = Denoise_Signal(cD, Threshold, J, param);

% proses menggabungkan data
n = 1;
for m = 1 : 2 : J * 2
    val = [];
    val(:, :, 1) = Cnew{1, m};
    val(:, :, 2) = Cnew{1, m + 1};
    wt.cfs{1, n} = val;
    n = n + 1;
end
appr = xapp.cfs{1, 5};
wt.cfs{1, n} = appr;

% Kenakan Transformasi Inverse wavelet Forward (inverse Wavelet biasa)
dden = idddtree(wt);

% Centering 
y = dden - mean(dden);

    %% Analisis Parameter
    % Menghitung error dengan MSE, SNR dan PSNR
    
    % Hitung MSE
    err1 = (norm(x(:)-y(:),2).^2)/numel(x);
    fprintf('>> The Mean-squared Error is %0.4f\n', err1);

    % Hitung SNR
    noiseampestimation = x-y;
    snr1 = 20*log10(rms(x)/rms(noiseampestimation));
    fprintf('>> The Signal Noise to ratio is %0.4f\n', snr1);
    
    % Hitung RMSE
    RMSE = sqrt(err1);
    fprintf('>> The RMSE is %0.4f\n', RMSE);

    %% Menampilkan hasil setiap langkah
    addpath('./plots');
    
    outfolder = 'Output Plots';
    if ~exist(outfolder, 'dir')
        mkdir(outfolder);
    end
    sname = fname(1:length(fname)-4);
    
    % plot sinyal asli
    foname = sprintf('%s_pcgAsli.jpg', sname);
    onam1 = fullfile(outfolder, foname);
    figure;
    ax2 = axes('Position',[0.14 0.17 0.78 0.74]);
    ax2.ActivePositionProperty = 'position';
    plot((1:length(x))/fs, x,'LineWidth', .5);
    xlabel('Waktu (sec)', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    ylabel('Amplitudo', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    title('Sinyal PCG - \it{Ground Truth}');
    set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 2, 'GridAlpha', 0.1);
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    export_fig (onam1, '-jpg', '-r200', '-a4', '-painters', '-transparent');
   
    % sinyal kena noise
    foname = sprintf('%s_pcgNoisy_Level_Dekomposisi%d_%s.jpg', sname, J, param);
    onam2 = fullfile(outfolder, foname);
    figure;
    ax2 = axes('Position',[0.14 0.17 0.78 0.74]);
    ax2.ActivePositionProperty = 'position';
    plot((1:length(xnoise))/fs, xnoise,'LineWidth', .5);
    xlabel('Waktu (sec)', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    ylabel('Amplitudo', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    title('Sinyal PCG - Tercampur Noise Acak');
    set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 2, 'GridAlpha', 0.1);
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    export_fig (onam2, '-jpg', '-r200', '-a4', '-painters', '-transparent');
    
    % denoised signal
    foname = sprintf('%s_pcgDenoised_Level_Dekomposisie%d_%s.jpg', sname, J,param);
    onam3 = fullfile(outfolder, foname);
    figure;
    ax2 = axes('Position',[0.14 0.17 0.78 0.74]);
    ax2.ActivePositionProperty = 'position';
    plot((1:length(y))/fs,y,'LineWidth', .5);
    xlabel('Waktu (sec)', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    ylabel('Amplitudo', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    title('Sinyal PCG - Hasil Denoising');
    set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 2, 'GridAlpha', 0.1);
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    export_fig (onam3, '-jpg', '-r200', '-a4', '-painters', '-transparent');
    
    % sinyal asli vs sinyal denoised
    foname = sprintf('%s_pcgAslivsDenoised_Level_Dekomposisi%d_%s.jpg', sname, J,param);
    onam4 = fullfile(outfolder, foname);
    figure;
    ax2 = axes('Position',[0.14 0.17 0.78 0.74]);
    ax2.ActivePositionProperty = 'position';
    plot((1:length(xnoise))/fs,xnoise,'LineWidth', .5); hold on;
    plot((1:length(y))/fs,y,'LineWidth', .5); hold off;    
    xlabel('Waktu (sec)', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    ylabel('Amplitudo', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    title('Sinyal PCG vs Hasil \it{Denoising}');
    set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 2, 'GridAlpha', 0.1);
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    export_fig (onam4, '-jpg', '-r200', '-a4', '-painters', '-transparent');
    legend('Noisy Signal','Denoised Signal')
        end
    end
end

