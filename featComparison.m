%%
%--------------------Read EEG Data for patient 01--------------------------
% Read first batch of data
parfor i = 1:20
    fileName = sprintf('chb01/chb01_%d.edf', i);
    [~,records(:,:,i-9)] = edfread(fileName);
end

% Each channel is represented by a column for filtering 
data = permute(records, [2, 1, 3]);

% Indices corresponding to seizure events 
seizureInd = [900*3+(749:759) 900*4+(367:373) 900*15+(433:443)...
    900*16+(254:266) 900*18+(430:452) 900*21+(82:105) 900*26+(465:491)];

% Create low pass filter 
Fs = 256;  % Sampling Frequency

N      = 50;      % Order
Fc     = 25;      % Cutoff Frequency
DpassU = 0.01;    % Upper Passband Ripple
DpassL = 0.01;    % Lower Passband Ripple
DstopU = 0.0001;  % Upper Stopband Attenuation
DstopL = 0.0001;  % Lower Stopband Attenuation

% Calculate the coefficients using the FIRCLS function.
b  = fircls(N, [0 Fc Fs/2]/(Fs/2), [1 0], [1+DpassU DstopU], [1-DpassL ...
            -DstopL]);
Hd = dfilt.dffir(b);

LPFiltered = filter(Hd, data);
records = permute(LPFiltered,[2 1 3]);
%%
rec3 = records(:,:,6);
rec4 = records(:,:,4);

matrix3 = ones(32, 900, 23);
matrix4 = ones(32, 900, 23);

wdSize = 1024;
level = 5;
wName = 'db4';

%% For third record
for j = 1:size(rec3, 2) / wdSize
    for k = 1:size(rec3, 1)
        signal = rec3(k,1 + (wdSize * (j - 1)):wdSize * j);
        wpt = wpdec(signal,level,wName);    % Wavelet packet tree
        for m = 1:2^level
            wpCoef = wpcoef(wpt, [level, m - 1]);
            matrix3(m, j, k) = log(sqrt(mean(wpCoef.*wpCoef)));
        end
    end
end
%%
for i = 1:23
    for j = 1:32
        f = figure;
        plot(matrix3(j, :, i))
        y = get(gca, 'ylim');
        
        %         % Add lines
        %         h1 = line([749 749], y);
        %         h2 = line([759 759], y);
        
        %         % Set properties of lines
        %         set([h1 h2],'Color','k','LineWidth',0.5)
        
        % Add a patch
        %         patch([749 759 759 749],[y(1) y(1) y(2) y(2)],'y')
        %         set(gca,'children',flipud(get(gca,'children')))
        
        xlabel('Window')
        ylabel('MSS of WP Coefficients')
        saveas(f, sprintf('patient01rec15feat/ch%dpacket%d.png', i, j));
        close(f)
    end
end
%%
% For fourth record

for j = 1:size(rec4, 2) / wdSize
    for k = 1:size(rec4, 1)
        signal = rec4(k,1 + (wdSize * (j - 1)):wdSize * j);
        wpt = wpdec(signal,level,wName);    % Wavelet packet tree
        for m = 1 :2^level
            wpCoef = wpcoef(wpt, [level, m - 1]);
            matrix4(m, j, k) = log(sqrt(mean(wpCoef.*wpCoef)));
        end
    end
end

%%
for i = 1:23
    for j = 1:32
        f = figure;
        plot(matrix4(j, :, i))
        
        xlabel('Window')
        ylabel('MSS of WP Coefficients')
        saveas(f, sprintf('patient01rec16feat/ch%dpacket%d.png', i, j));
        close(f)
    end
end

%% Comparison
for i = 1:23
    for j = 1:32
        f = figure;
        subplot(211)
        plot(matrix3(j, :, i))
        xlabel('Window')
        ylabel('MSS of WP Coefficients')
        
        subplot(212)
        plot(matrix4(j, :, i))
        xlabel('Window')
        ylabel('MSS of WP Coefficients')
        
        saveas(f, sprintf('patient01rec3and4comp/ch%dpacket%d.png', i, j));
        close(f)
    end
end

% %
% mat3 = permute(matrix3, [2 1 3]);
% mat4 = permute(matrix4, [2 1 3]);
% 
% var3 = var(mat3);
% var3 = permute(var3, [2, 3, 1]);
% var4 = var(mat4);
% var4 = permute(var4, [2, 3, 1]);
% 
% mean3 = mean(mat3);
% mean3 = permute(mean3, [2, 3, 1]);
% mean4 = mean(mat3);
% mean4 = permute(mean4, [2, 3, 1]);
% 
% coefVar3 =  var3./abs(mean3);
% coefVar4 =  var4./abs(mean4);
% 
% [M3, I3] = max(coefVar3);
% [M4, I4] = max(coefVar4);
% 
% [~, maxInd3] = sort(M3, 'descend');
% [~, maxInd4] = sort(M4, 'descend');
% 
% promChan = maxInd3(1:8);
% promWvPkt = I3(promChan);
% 
% Indices corresponding to seizure events 
% seizureInd = 433:443;
% seizureInd1 = 254:266;
% 
% featureVector3 = [];
% featureVector4 = [];
% for i = 1:length(promChan)
%     featureVector3 = [featureVector3; matrix3(promWvPkt(i),:,promChan(i))];
%     featureVector4 = [featureVector4; matrix4(promWvPkt(i),:,promChan(i))];
% end
% 
% Preparing the training dataset and target vectors
% targetVectorNonSeiz = zeros(1, 1800 - length([seizureInd (seizureInd1 + 900)]))';
% targetVectorSeiz = ones(1, length([seizureInd (seizureInd1 + 900)]))';
% 
% featureVectorSeiz = [featureVector3(:,seizureInd) featureVector4(:,seizureInd1)];
% featureVectorSeiz = featureVectorSeiz';
% 
% featureVectorNonSeiz = [featureVector3 featureVector4];
% featureVectorNonSeiz(:,[seizureInd (seizureInd1 + 900)]) = [];
% featureVectorNonSeiz = featureVectorNonSeiz';
% 
%% Create and train SVM model
% Randomly select 200 samples from non-seizure data
% randSampleInd = randi(size(featureVectorNonSeiz, 1), 1, 500);
% finalFeatureVector = [featureVectorSeiz;featureVectorNonSeiz(randSampleInd,:)];
% finalTargetVector = [targetVectorSeiz;targetVectorNonSeiz(randSampleInd,:)];
% 
% SVMmodel = fitcsvm(finalFeatureVector, finalTargetVector, 'Standardize',...
%     true,'KernelFunction','RBF','KernelScale','auto');

%% Prominent channels for the first patient
% Channels 2, 6, 22, 9
% Corresponding wavelet packets 13, 9, 6, 9
