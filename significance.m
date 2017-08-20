%%
%--------------------Read EEG Data for patient 03--------------------------
% Read first batch of data
parfor i = 1:10
    fileName = sprintf('chb03/chb03_%d.edf', i);
    [~,records(:,:,i)] = edfread(fileName);
end
%%
% Each channel is represented by a column for filtering 
% Only channel with seizure data are considered, therefore indexing used
% will br from 1-6 later
records = records(:,:,[1 2 3 4 5 6]);
data = permute(records, [2, 1, 3]);

% Indices corresponding to seizure events 
% Patient 01
% seizureInd = [900*2+(749:759) 900*3+(367:373) 900*14+(433:443)...
%     900*15+(254:266) 900*17+(430:452)];

% Indices corresponding to seizure events
% Patient 03
% seizureInd = [(91:104) 900*1+(183:199) 900*2+(108:125) 900*3+(541:554)...
%     900*4+(496:507) 900*5+(437:445)];

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
%% Consider samples of seizure and non-seizure events 
Fs = 256;

seiz1 = records(:,362*Fs+1:414*Fs,1);
seiz2 = records(:,731*Fs+1:796*Fs,2);
seiz3 = records(:,432*Fs+1:501*Fs,3);
seiz4 = records(:,2162*Fs+1:2214*Fs,4); 
seiz5 = records(:,1982*Fs+1:2029*Fs,5);
seiz6 = records(:,1725*Fs+1:1778*Fs,6);

seiz(:,:,1) =seiz1(:,1:12032);
seiz(:,:,2) =seiz2(:,1:12032);
seiz(:,:,3) =seiz3(:,1:12032);
seiz(:,:,4) =seiz4(:,1:12032);
seiz(:,:,5) =seiz5(:,1:12032);
seiz(:,:,6) =seiz6(:,1:12032);

nonSeiz(:,:,1) = records(:,1:47*Fs,1);
nonSeiz(:,:,2) = records(:,1500*Fs+1:1547*Fs,1);
nonSeiz(:,:,3) = records(:,100*Fs+1:147*Fs,2);
nonSeiz(:,:,4) = records(:,2000*Fs+1:2047*Fs,2);
nonSeiz(:,:,5) = records(:,100*Fs+1:147*Fs,3);
nonSeiz(:,:,6) = records(:,2100*Fs+1:2147*Fs,3);
nonSeiz(:,:,7) = records(:,1:47*Fs,4);
nonSeiz(:,:,8) = records(:,2500*Fs+1:2547*Fs,4);
nonSeiz(:,:,9) = records(:,700*Fs+1:747*Fs,5);
nonSeiz(:,:,10) = records(:,2500*Fs+1:2547*Fs,5);
nonSeiz(:,:,11) = records(:,1000*Fs+1:1047*Fs,6);
nonSeiz(:,:,12) = records(:,2000*Fs+1:2047*Fs,6);
%% For non seizure data - Feature computation
wdSize = 1024;
level = 5;
wName = 'db4';

for r = 1:size(nonSeiz, 3)
    for j = 1:size(nonSeiz, 2) / wdSize
        for k = 1:size(nonSeiz, 1)
            signal = nonSeiz(k,1 + (wdSize * (j - 1)):wdSize * j, r);
            wpt = wpdec(signal,level,wName);
            for m = 1:2^level
                wpCoef = wpcoef(wpt, [level, m - 1]);
                nonSeizFeatures(j, k, m, r) = log(sqrt(mean(wpCoef.*wpCoef)));
            end
        end
    end
end

%% For seizure data - Feature computation
for r = 1:size(seiz, 3)
    for j = 1:size(seiz, 2) / wdSize
        for k = 1:size(seiz, 1)
            signal = seiz(k,1 + (wdSize * (j - 1)):wdSize * j, r);
            wpt = wpdec(signal,level,wName);
            for m = 1:2^level
                wpCoef = wpcoef(wpt, [level, m - 1]);
                seizFeatures(j, k, m, r) = log(sqrt(mean(wpCoef.*wpCoef)));
            end
        end
    end
end

%% Anova testing for each feature
% A matrix of p values is generated for each feature and the feature with
% maximum number of p values less than 0.001 are considered as significant
% features
for i = 1:1:size(seizFeatures, 2)
    for j = 1:1:size(seizFeatures, 3)
        p = [];
        for k = 1:1:size(nonSeizFeatures, 4)
            for h = 1:1:size(seizFeatures, 4)
                p(k, h) = anova1([nonSeizFeatures(:,i, j, k)...
                    seizFeatures(:,i, j, h)], [], 'off');
                close all 
            end
        end
        signifMatrix = p < 0.001;
        sigScore(i, j) = sum(sum(signifMatrix));
    end
end
