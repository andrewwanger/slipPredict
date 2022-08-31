%% Load stuff
load('slip_data_v2.mat')
%% Plot data to find start-end points
D = Differential_2_4;
figure, hold on, plot(1:length(D),D)
L = 20; %window
S = L/2; v = nan(1,length(D)-L);
for i = 1:length(D) - L
    v(i) = mean(D(i:i+L));
end
plot(1+S:length(v)+S,v)
% set figure params
ylabel('Cap Count'),xlabel('sample'),legend('North-South'),set(gcf,'Position',[201  369  1706  344])
figure, heatmap(flipud(AbsFFT(1:30,:)))
%% Segment visualization
%--------------------- Select Data ----------------------------------------
%index: {1}=slow[125,600], {2}=med[250,1200], {3}=fast[500,2400]
D_temp = slip_data_v2{3};
doPlot = 1; %1 = yes, 0 = no
doDispFrames = 0; %1 = yes, 0 = no
%--------------------- Plot Data ------------------------------------------
D = D_temp.data(3:end);C = D_temp.cycles;sensingMode=D_temp.sensingMode;windowSize=D_temp.windowSize;overlapSize=D_temp.overlapSize;overlapTime=D_temp.overlapTime;windowDuration=D_temp.windowDuration;thisFreq=D_temp.thisFreq;FreqBinNum=D_temp.FreqBinNum;slopeCompensateOn=D_temp.slopeCompensateOn;titleString=D_temp.titleString;freqLim=D_temp.freqLim;clim=D_temp.clim;s = D_temp.first_start;e = D_temp.first_end;d = D_temp.last_start;
% s_a=D_temp.AbsFFT_first_start;e_a=D_temp.AbsFFT_first_end;d_a=D_temp.AbsFFT_last_start;

if doPlot == 1
    figure
    AbsFFT = FFT_plot_yuri_thr(sensingMode,D,windowSize,overlapSize,overlapTime,windowDuration,thisFreq,FreqBinNum,slopeCompensateOn,titleString,freqLim,clim);
    clf
    hold on
    plot(1:length(D),D)
else
    figure
    AbsFFT = FFT_plot_yuri_thr(sensingMode,D,windowSize,overlapSize,overlapTime,windowDuration,thisFreq,FreqBinNum,slopeCompensateOn,titleString,freqLim,clim);
    close
end

L = 20; S = L/2;
v = nan(1,length(D)-L);
for i = 1:length(D) - L
    v(i) = mean(D(i:i+L));
end

if doPlot == 1
    plot(1+S:length(v)+S,v)
    X = nan(1,C);
end

slip_period = (d-s)/(C-1);
% slip_period_a = (d_a-s_a)/(C-1);

% clearvars ML_ds
% k = 0;
for i = 1:C
    if doPlot == 1
        plot([s,s]+round(slip_period*(i-1)),[-50,300],'g')
        plot([e,e]+round(slip_period*(i-1)),[-50,300],'r')
    end
    
    offset = (e-s)/8;
    x = round(s + slip_period*(i-1) + offset);
    offset = (e-s)/6;
    y = round(e + slip_period*(i-1) - offset);
    if doPlot == 1
        plot([x,x],[-50,300],'g')
        plot([y,y],[-50,300],'r')
    end
    a = x; b = y;
    
    %---------Make ML dataset----------
    fft_dat = getEasyAbsFFT(D(x:y),windowSize,1,thisFreq,FreqBinNum,slopeCompensateOn);
%     k = k + 1;
%     ML_ds{k}.fft_data = fft_dat;
%     ML_ds{k}.original_data = D(x:y);
%     ML_ds{k}.speed = D_temp.speed;
%     ML_ds{k}.material = D_temp.material;
    %----------------------------------
    
    temp = abs(lookup_fft-x);
    x = find(temp == min(temp));
    temp = abs(lookup_fft-y);
    y = find(temp == min(temp));
    if doDispFrames == 1
        disp(['x:',num2str(x),', y:',num2str(y),'.....[',num2str(a),',',num2str(b),']'])
    end
    
%     offset = 1;
%     x = round(s_a + slip_period_a*(i-1) + offset);
%     y = round(e_a + slip_period_a*(i-1) - offset);
%   
    if doPlot == 1
        X(i) = sum(sum(fft_dat(1:10,:))) / size(fft_dat,2);
        text(s+slip_period*(i-1),-50,num2str(round(X(i))))
        text(s+round(slip_period*(i-1)),250,num2str(i))
    end

%     AbsFFT(1:10,[s:e])
%     
%     round([2483:2966]/length(D)*100);
    
end
if doPlot == 1
    title(['speed: ',num2str(D_temp.speed),', average weight: ',num2str(mean(X))])
    ylabel('Cap Count')
    xlabel('sample')
    legend('North-South')
    set(gcf,'Position',[201  369  1706  344])
end

%% Making data files
%index: cherry   {1}=slow[125,600], {2}=med[250,1200], {3}=fast[500,2400],
%                {4}=[20,10], {5}=[40,100], {6}=[80,200], {7}=[160,400]
%                {8}=[20,50], {9}=[30,70], {10}=[50,120], {11}=[60,140]
%                {12}=[70,170], {13}=[90,240], {14}=[100,300]
%                {15}=[110,400], {16}=[120,520], {17}=[130,660]
%NEW index:
%       cherry   {1}=slow[10,700], {2}=med[15,700], {3}=fast[20,700]
%                {4}=slow[25,700], {5}=med[30,700], {6}=fast[35,700]
%                {7}=slow[40,700], {8}=med[45,700], {9}=fast[50,700]
%                {10}=slow[55,700], {11}=med[60,700], {12}=fast[65,700]
%                {13}=slow[70,700], {14}=med[75,700], {15}=fast[80,700]
%                {16}=slow[85,700], {17}=med[90,700], {18}=fast[95,700]
%                {19}=slow[100,700]
%       basswood {20}=slow[10,700], {21}=med[15,700], {22}=fast[20,700]
%                {23}=slow[25,700], {24}=med[30,700], {25}=fast[35,700]
%                {26}=slow[40,700], {27}=med[45,700], {28}=fast[50,700]
%                {29}=slow[55,700], {30}=med[60,700], {31}=fast[65,700]
%                {32}=slow[70,700], {33}=med[75,700], {34}=fast[80,700]
%                {35}=slow[85,700], {36}=med[90,700], {37}=fast[95,700]
%                {38}=slow[100,700]
%       acrylic  {39}=slow[10,700], {40}=med[15,700], {41}=fast[20,700]
%                {42}=slow[25,700], {43}=med[30,700], {44}=fast[35,700]
%                {45}=slow[40,700], {46}=med[45,700], {47}=fast[50,700]
%                {48}=slow[55,700], {49}=med[60,700], {50}=fast[65,700]
%                {51}=slow[70,700], {52}=med[75,700], {53}=fast[80,700]
%                {54}=slow[85,700], {55}=med[90,700], {56}=fast[95,700]
%                {57}=slow[100,700]

% data_to_read = [8,9,5,10,11,12,6,13]; data is within 300 frames

ind = 39;
slip_data_v2{ind}.data = Differential_2_4;
slip_data_v2{ind}.first_start = 3478;
slip_data_v2{ind}.first_end = 9582;
slip_data_v2{ind}.last_start = 60540;
slip_data_v2{ind}.cycles = 7;
slip_data_v2{ind}.material='acrylic';
slip_data_v2{ind}.speed = 10;

slip_data_v2{ind}.sensingMode=sensingMode;
slip_data_v2{ind}.windowSize=windowSize;
slip_data_v2{ind}.overlapSize=overlapSize;
slip_data_v2{ind}.overlapTime=overlapTime;
slip_data_v2{ind}.windowDuration=windowDuration;
slip_data_v2{ind}.thisFreq=thisFreq;
slip_data_v2{ind}.FreqBinNum=FreqBinNum;
slip_data_v2{ind}.slopeCompensateOn=slopeCompensateOn;
slip_data_v2{ind}.titleString=titleString;
slip_data_v2{ind}.freqLim=freqLim;
slip_data_v2{ind}.clim=clim;

%% Extra plots and stuff - raw amplitude
figure, hold on
x_dat = [];
y_dat = [];
for i = 39:45
    data = slip_data_v2{i}.data;
    speed = slip_data_v2{i}.speed;
    C = slip_data_v2{i}.cycles;
    s = slip_data_v2{i}.first_start;
    e = slip_data_v2{i}.first_end;
    d = slip_data_v2{i}.last_start;
    slip_period = (d-s)/(C-1);
    a = [];
    for j = 1:C
        offset = (e-s)/8;
        x = round(s + slip_period*(j-1) + offset);
        offset = (e-s)/6;
        y = round(e + slip_period*(j-1) - offset);
        a = [a;mean(data(x:y))];
%         a = [a;data(x:y)];
        scatter(repmat(speed+j/20,1,length(data(x:y))),data(x:y),'filled')
    end
    x_dat = [x_dat,repmat(speed,1,length(a))];
    y_dat = [y_dat,a'];
    
%     scatter(repmat(speed,1,length(a)),a,'filled')
end
ylim([0,200])

figure
x = x_dat; 
y1 = y_dat; 
scatter(x,y1,25,'b','*') 
P = polyfit(x,y1,1);
yfit = P(1)*x+P(2);
hold on;
plot(x,yfit,'r-.');

figure
makeReg(x_dat,y_dat,[0,max(y_dat)],slip_data_v2{i}.material) % both should be row vectors
% set(gcf,'Position',[2401.0000  486.6000  154.4000  147.2000])

%% Extra plots and stuff - fft band
select_band = 301;

% figure, hold on
x_dat = [];
y_dat = [];
for i = 1:19
    D = slip_data_v2{i}.data;
    speed = slip_data_v2{i}.speed;
    C = slip_data_v2{i}.cycles;
    s = slip_data_v2{i}.first_start;
    e = slip_data_v2{i}.first_end;
    d = slip_data_v2{i}.last_start;
    slip_period = (d-s)/(C-1);
    a = [];
    for j = 1:C
        offset = (e-s)/8;
        x = round(s + slip_period*(j-1) + offset);
        offset = (e-s)/6;
        y = round(e + slip_period*(j-1) - offset);
        
        fft_dat = getEasyAbsFFT(D(x:y),windowSize,1,thisFreq,FreqBinNum,slopeCompensateOn);
        
        a = [a,mean(fft_dat,2)];
%         a = [a,fft_dat];
    end
    x_dat = [x_dat,repmat(speed,1,size(a,2))];
    y_dat = [y_dat,a];
    
%     scatter(repmat(speed,1,length(a)),a,'filled')
end
%% Single Hz analysis
fit_hz = zeros(1,size(y_dat,1));
for i = 1:size(y_dat,1)
    [~,~,Rsq] = linfit(x_dat,y_dat(i,:));
    fit_hz(i) = Rsq;
end
disp(['min:',num2str(min(fit_hz)),', max:',num2str(max(fit_hz))])
%% Band Hz analysis
a = [];
b1 = [];
b2 = [];
k = 0;
slopes = [];
for j = 0:301
    fit_hz = zeros(1,size(y_dat,1)-j);
    for i = 1:size(y_dat,1)-j
        [slope,~,Rsq] = linfit(x_dat,sum(y_dat(i:i+j,:),1));
        fit_hz(i) = Rsq;
        slopes = [slopes slope];
    end
%     disp(['band size:',num2str(j+1),' [min,max]:[',num2str(min(fit_hz)),',',num2str(max(fit_hz)),']'])
    k = k + 1;
    a(k) = j+1;
    b1(k) = min(fit_hz);
    b2(k) = max(fit_hz);
end
figure,hold on,scatter(a,b2),scatter(a,b1,'s')
% set(gcf,'Position',[244.2000  357.0000  328.0000  196.0000])
set(gcf,'Position',[374.6000  292.2000  184.8000  153.6000])
title('fit values vs. frequency bands')
xlabel('frequency band (Hz)')
ylabel('R^2')
legend('maximum R^2','minimum R^2')
%% plot single Hz's
st = 281;
figure('units','normalized','outerposition',[0 0 1 1])
k = 0;
for i = st:st+19
    k = k + 1;
    subplot(4,5,k)
    makeReg(x_dat,y_dat(i,:),[0,max(y_dat(i,:))],slip_data_v2{1}.material) % both should be row vectors
    title([num2str(i),' Hz'])
end
%% plot bands of 10 Hz's
st = 1;
figure('units','normalized','outerposition',[0 0 1 1])
k = 0;
for i = 101:10:300
    k = k + 1;
    subplot(4,5,k)
    y_temp = sum(y_dat(i:i+10,:));
    
    makeReg(x_dat,y_temp,[0,max(y_temp)],slip_data_v2{1}.material) % both should be row vectors
    title([num2str(i),'-',num2str(i+10),' Hz'])
end