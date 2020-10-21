%% Program Untuk Menghitung Noise Level Estimation
function [nLevel, xs] = NoiselevelEstimation(x)

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
 
 % Second and subsequent stages filter coefficients - Tree R
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

% Second and Subsequent stages filter coefficients - Tree I
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

%% Proses Dual-tree Discrete Wavelet Transform
J = 1;
wt = dddtree('cplxdt',x,J,Faf,af);
outputindices = {[1 1]; [1 2];};
% Ekstrak Koefisien Detail
cD1 = dddtreecfs('e',wt,'ind',outputindices);
% assignin('base', 'cD', cD1); % cara menampilkan isi dari variable yang
% ada di dalam fungsi ke workspace matlab
xr = cell2mat(cD1);
xs = median(abs(xr));
nLevel = xs / 0.6745;
