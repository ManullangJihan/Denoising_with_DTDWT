
function Df = Denoise_Signal(D, NoiselevelEstimation, J, param)
% fungsi ini bertujuan melakukan filtering/denoising menggunakan 3 metode
% yaitu hard thresholding dan soft thresholding

Df = D;
switch (param)
    case 'h' % hard thresholding
        % untuk yang real detail components
        for m = 1 : 2 : size(D, 2)
            x = D{1, m};
            if m >= 2 * J - 1
                xfr = x .* (abs(x) > NoiselevelEstimation(1));
            else
                xfr = zeros(size(x));
            end
            Df{1, m} = xfr;
        end
        
        % untuk yang imaginary detail components
        for m = 2 : 2 : size(D, 2)
            x = D{1, m};
            if m >= 2 * J - 1
                xfi = x .* (abs(x) > NoiselevelEstimation(2));
            else
                xfi = zeros(size(x));
            end
            Df{1, m} = xfi;
        end
        
    case 's' % soft thresholding
        % untuk yang real detail components
        for m = 1 : 2 : size(D, 2)
            x = D{1, m};
            if m >= 2 * J - 1
                temp = abs(x) - NoiselevelEstimation(1);
                temp = (temp + abs(temp))/2;
                xfr = sign(x) .* temp;
            else
                xfr = zeros(size(x));
            end
            Df{1, m} = xfr;
        end
        
        % untuk yang real detail components
        for m = 2 : 2 : size(D, 2)
            x = D{1, m};
            if m >= 2 * J - 1
                temp = abs(x) - NoiselevelEstimation(2);
                temp = (temp + abs(temp))/2;
                xfi = sign(x) .* temp;
            else
                xfi = zeros(size(x));
            end
            Df{1, m} = xfi;
        end
        
    otherwise
        error('Parameter Keliru\n');
end