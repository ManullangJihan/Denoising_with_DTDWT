% Program utama untuk Dual Tree Discrete Wavelet Transform terhadap sinyal FPCG
% Nama  : Tuah Jihan
% prodi : S1 TT 

% Environment 
warning off;
clear all;
close all;
clc;


%% Memilih folder untuk menyimpan

[fname, pname] = uigetfile('*.wav', 'Pilih sebuah data PCG');
 
if ~isequal(fname, 0) || ~isequal(pname, 0)
  

%% Import data
pcgfile = fullfile(pname, fname);
[x, fs] = audioread(pcgfile);
fprintf('Processing: %s\n', fname);
    
%% Index data selection
% Gunakan data ntuk t detik saja
% t1 = 3;
% t2 = 6;
% N1 = round(t1 * fs);
% N2 = round(t2 * fs);
% x = x(N1 : N2-1);
x = x(1:512);
% Centering
x = x - mean(x);

% Normalisasi
x = x ./ max(abs(x));

%% Tambahkan Noise acak
snrawgn = 5;
datax = awgn(x,snrawgn,'measured');
xnoise = x + datax;

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
J = 4;
wt = dddtree('cplxdt',xnoise,J,Faf,af);

% Ekstrak Koefisien Detail dari Tree R dan Tree I, pada setiap Level
outputindices = {[1 1]; [1 2]; [2 1]; [2 2]; [3 1];[3 2]; [4 1]; [4 2]};
cD = dddtreecfs('e',wt,'ind',outputindices);
xapp = dddtreecfs('e',wt,'lowpass');

% Proses Denoising
Cnew = Denoise_Signal(cD, Threshold, J, param);

% Proses Menggabungkan Data
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

% Normalisasi
y = y ./ max(abs(y));

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
    
    % Denoised Signal
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
    
    % Noisy Signal vs Denoised Signal
    foname = sprintf('%s_pcgAslivsDenoised_Level_Dekomposisi%d_%s.jpg', sname, J,param);
    onam4 = fullfile(outfolder, foname);
    figure;
    ax2 = axes('Position',[0.14 0.17 0.78 0.74]);
    ax2.ActivePositionProperty = 'position';
    plot((1:length(xnoise))/fs,xnoise,'LineWidth', .5); hold on;
    plot((1:length(y))/fs,y,'r','LineWidth', 2); hold off;    
    xlabel('Waktu (sec)', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    ylabel('Amplitudo', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    title('Sinyal Noise vs Hasil \it{Denoising}');
    legend('Noisy Signal','Denoised Signal');
    set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 2, 'GridAlpha', 0.1);
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    export_fig (onam4, '-jpg', '-r200', '-a4', '-painters', '-transparent');
    end
%% END