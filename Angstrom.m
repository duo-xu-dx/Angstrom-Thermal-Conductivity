% Thermal conductivity calculation using Angstrom method

% Input: 1D tempreature variation in .csv format in the working directory.
% Each row is the temperature spatial distribution at a given frame.
% Output: a .csv file with "results_" prefix, summarizing the thermal
% diffusivity ad conductivity measured at each pixel relative to the first.

% Duo Xu 10/09/2023

%% Clear the workspace
clear all; close all; clc;

%% Parameters
% Settings
heatFreq = 0.01; % Hz
pixelLength = 0.17; % mm
cameraFreq = 9; % fps
timestep = 1/cameraFreq; % seconds
powerCriterion = 0.2;
sampleStartPixle = 3;
loadingStartFrame = 2;

% Material constants
density = 875; % kg/m^3
cp = 2200; % J/kg/K

%% Find all .csv files in the current directory
files = dir('*.csv');
csvfiles = {files.name};

%% Processing
fprintf(['Processing starts.\n',repmat('-',1,54),'\n']);
for i = 1:length(csvfiles)
    filename = csvfiles{i};

    % Skip files that start with "Result_"
    if startsWith(filename, 'results_')
        continue;
    end
    fprintf([repmat('-',1,54),'\n']);
    fprintf('Processing file %s\n', filename);
    dataArray = [];

    % Extracting data
    fid = fopen(filename, 'r');
    for j = 1:16
        fgetl(fid);
    end
    line = fgetl(fid);
    while ischar(line) && numel(strsplit(line,',')) > 2
        parts = strsplit(line, ',');
        timestamp = str2double(parts{1});
        dataArray = [dataArray; timestamp,...
            cellfun(@str2double,parts(2:end-1))];
        line = fgetl(fid);
    end
    fclose(fid);
    timeArray = dataArray(loadingStartFrame:end, 1);
    dataArray = dataArray(loadingStartFrame:end, sampleStartPixle:end);

    % Remove outliers in IR camera reading
    samplestd = std(dataArray,0,'all');
    for j = 2:size(dataArray, 1)
        for k = 1:size(dataArray, 2)
            if abs(dataArray(j,k)-dataArray(j-1,k))>samplestd
                dataArray(j, k) = dataArray(j-1,k);
            end
        end
    end

    % Find correct discretized FFT length
    mathlcm = lcm(1/heatFreq, cameraFreq);
    dataPoints = floor(length(timeArray)/mathlcm)*mathlcm;
    timeArray = timeArray(1:dataPoints);
    dataArray = dataArray(1:dataPoints,:);
    totalTime = timeArray(end)-timeArray(1);

    % Check power spectrum
    pixelIdxList = [];
    deltaTArray = [];
    amplitudeArray = [];
    phaseArray = [];
    freq = (0:length(timeArray)-1)/(timestep*length(timeArray));
    freqIDX = find(abs(freq-heatFreq) == min(abs(freq-heatFreq)),1);
    meanTempArray = [];
    for j = 1:size(dataArray,2)
        deltaT = dataArray(:,j)-mean(dataArray(:,j));
        deltaTfft = fft(deltaT);
        powerSpec = abs(deltaTfft).^2;
        if powerSpec(freqIDX)/sum(powerSpec)>powerCriterion
            pixelIdxList = [pixelIdxList,j];
            deltaTArray = [deltaTArray, deltaT];
            meanTempArray = [meanTempArray, mean(dataArray(:,j))];
            amplitudeArray = [amplitudeArray,abs(deltaTfft(freqIDX))...
                * 2/length(timeArray)];
            phasefft = -angle(deltaTfft(freqIDX))/(2*pi*heatFreq);
            if phasefft < 0
                phasefft = phasefft+1/heatFreq;
            end
            if phasefft > 1/heatFreq
                phasefft = phasefft-1/heatFreq;
            end
            if ~isempty(phaseArray) &&...
                    phasefft-phaseArray(end) < -1/(2*heatFreq)
                phasefft = phasefft+1/heatFreq;
            end
            phaseArray = [phaseArray,phasefft];
        else
            break;
        end
    end

    fprintf('\tUsable pixel number %d\n',length(pixelIdxList));
    amplitudeArray = transpose(amplitudeArray);
    phaseArray = transpose(phaseArray);
    distanceArray = transpose(pixelLength*(pixelIdxList-1));
    meanTempArray = meanTempArray';
    deltaTArray = deltaTArray';

    dt = phaseArray - phaseArray(1);
    l = distanceArray;
    alpha = (l.^2) ./ (2.* dt.*log(amplitudeArray(1)./amplitudeArray));
    % Convert to m^2/s
    alphanew = alpha*1e-6;
    k = alphanew.*density.*cp;

    % Create result table and save to CSV
    result = table(distanceArray,amplitudeArray,phaseArray,alpha,k,...
        repmat(powerCriterion,size(distanceArray)),...
        repmat(sampleStartPixle,size(distanceArray)),...
        repmat(density,size(distanceArray)),...
        repmat(cp, size(distanceArray)), ...
        meanTempArray,...
        'VariableNames',{'distance_mm','amplitude_K','phase',...
        'Diffusivity_mm2_s','conductivity_W_mK','Power_Criterion',...
        'Start_Pixel','Density_kg_m3','Specific_Heat_J_kg_K',...
        'pixelMeanTemperature'});
    writetable(result, ['results_',filename]);
end
fprintf([repmat('-',1,54),'\n',repmat('-',1,54),'\n',...
    'Processing complete.\n']);