function [reconstructedSignal] = fFourierDecomposition(time,inputSignal,numberOfComponents)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
time = time - time(1);

L = length(time);
sampleT  = time(2)-time(1);
Fs = 1/sampleT;

fftSignal = fft(inputSignal);
fftMagnitude = abs(fftSignal/L);
fftMagnitude = fftMagnitude(1:floor(L/2+1));
fftMagnitude(2:end-1) = 2*fftMagnitude(2:end-1);
fftPhase = angle(fftSignal);

f = Fs*(0:(L/2))/L;

mainComponents = unique(fftMagnitude);
numOfMainComponents = numberOfComponents;
fourierSignal = zeros(size(inputSignal));

magniComponent = zeros(numOfMainComponents,1);
mainFreq = zeros(numOfMainComponents,1);
phaseComponent  = zeros(numOfMainComponents,1);

fourierSignalAnalytic = 0;

syms x;

for i = 1:numOfMainComponents
    magniComponent(i) = mainComponents(end-numOfMainComponents+i);
    indxFreq = find(fftMagnitude == magniComponent(i));
    mainFreq(i) = f(indxFreq);
    phaseComponent(i) = fftPhase(indxFreq);
    fourierSignal = magniComponent(i)*cos(2*pi*mainFreq(i)*time+phaseComponent(i))+fourierSignal;
    fourierSignalAnalytic = magniComponent(i)*cos(2*pi*mainFreq(i)*x+phaseComponent(i))+fourierSignalAnalytic;
end

reconstructedSignal.magnitude = magniComponent;
reconstructedSignal.freq = mainFreq;
reconstructedSignal.phase = phaseComponent;
reconstructedSignal.signal = fourierSignal;
reconstructedSignal.signalAnalytic = fourierSignalAnalytic;

figure, hold on, grid on
title('Original vs Fourier version')
plot(time,inputSignal)
plot(time,fourierSignal)
xlabel('Time [s]')
legend('InputSignal','Fourier approximation')

% figure
% plot(f,fftMagnitude) 
% title('Single-Sided Amplitude Spectrum of Input(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')

end

