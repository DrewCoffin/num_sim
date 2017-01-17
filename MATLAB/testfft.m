function[] = testfft(noiseamp, duration, sigamp, sigfreq, varargin)
Fs = 1000; %sampling frequency
T = 1/Fs;
L = 1000*duration;
t = (0:L-1)*T;
curve = sigamp*sin(2*pi*sigfreq*t);
%celldisp(varargin)
for i = [1:length(varargin)]
    if rem(i,2) == 0
        freq = varargin{i};
        amp = varargin{i-1};
        curve = curve + amp*sin(2*pi*freq*t);
    end
end
C = curve + noiseamp*randn(size(t));
Y = fft(C);
P2 = abs(Y/L); %two-sided spectrum
P1 = 2*P2(1:L/2+1); %single-sided spectrum
f = Fs*(0:(L/2))/L; %rescale frequency axis
figure
ax1 = subplot(2,1,1);
ax2 = subplot(2,1,2);
plot(ax1,t,C)
title(ax1,'Original')
plot(ax2,f,P1)
title(ax2,'Fourier Transform')
end