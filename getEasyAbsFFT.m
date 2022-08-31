function AbsFFT = getEasyAbsFFT(data,windowSize,overlapSize,thisFreq,FreqBinNum,slopeCompensateOn)

if length(data) < 300
    error('data file is <300')
end

% disp(['windowSize = ',num2str(windowSize)])
% disp(['overlapSize = ',num2str(overlapSize)])

t = (1:floor((length(data)-windowSize + 1) / overlapSize))-1;
f = linspace(0,thisFreq/2,FreqBinNum+1);
[TimeBin,~] = meshgrid(t,f);
AbsFFT = TimeBin * 0;
%
for idx = 1:length(t)
    
    % StartingPoint = windowSize * (idx-1) + 1;
    StartingPoint = overlapSize * (idx-1) + 1;
    
    
    thisData = data(StartingPoint:StartingPoint+windowSize-1);%.* hamming(windowSize);
    
    if slopeCompensateOn
        %%Try Linear Fit Subtraction
        x = (1:length(thisData))';
        p = polyfit(x,thisData,1);
        fittedData = polyval(p,x);
        thisData = thisData - fittedData;
    else
        
        thisData = thisData - mean(thisData);
    end
    
    
    thisData = thisData .* hamming(windowSize);
    
    
    Fs = thisFreq;
    T = 1/Fs;
    L = FreqBinNum * 2;
    
    
    S = thisData;
    
    Y = fft(S,L);
    
    f = Fs*(0:(L/2))/L;
    P2 = abs(Y);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    AbsFFT(:,idx) = P1;
    
end