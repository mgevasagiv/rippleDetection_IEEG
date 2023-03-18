classdef RippleDetector_class < handle
    
    properties
        
        samplingRate = 1000;
        dataFilePrefix = 'CSC';
        
        %all the parameters are from Staresina et al (and Zhang et al is identical)
        minFreq = 80;
        maxFreq = 100;
        
        % ripple detection on MICRO channels - parameters as in Le Van
        % Quyen et al 2008
        % minFreq = 80;
        % maxFreq = 200;
        
        RMSWindowDuration = 20; %ms
        rippleThreshPercentile = 99;
        minDurationAboveThresh = 38; %ms
        minNumOfExtreme = 3; %min number of peaks / troughs in a ripple event
        NPointsToAverage = 3;
        minPercNaNAllowed = 0.1;
        minDistBetweenRipples = 20; %ms
        
        %IIS removal constants
        windowAroundIIS = 500; %ms
        
        %STIMULATION removal constants
        windowAroundSTIM = 200; %ms
      
        %sleep scoring parameters
        scoringEpochDuration = 0.001; % How many seconds represented by one individual value in the scoring vector [scalar].
        sleepEpochs = [1]; % all the values in the scoring vector which represent sleep stages for which we want to perform the analysis (like NREM/REM/transitions) [1D vector].
        
        %filtering constants
        defaultFilterOrder = 3;
        nanWarning = 0.01;
        
        %plot params
        secondBefAfter = 0.5; %seconds
        subplotSizeY = 4;
        subplotSizeX = 3;
        nInPlotMicro = 4;
        
        spikeMultiUnits = true;
        
        %control params for spike rate around ripples
        minDistControl = 300; %ms
        maxDistControl = 1300; %ms
        
        firingRateWinSize = 10; %ms
        windowSpikeRateAroundRip = 500; %ms
        windowForSignificance = 100; %ms
        
        %params for spike rate around stimulations
        windowSpikeRateAroundStim = 500;
        windowSpikeRateForComparison = 500;  %ms - for comparing between stim and control
        controlDistForStim = 1000; %ms
        
        avgRippleBeforeAfter = 1; %second
        freqoiForAvgSpec = [0:0.5:10];
        freqRangeForAvgSpec = [1:250];
        timeBeforeAfterEventRipSpec = 1; %second
        timeForBaselineRip = 1; %second
        minNCycles = 5;
        minWinSizeSpec = 100; %ms
        minNripples = 5;
        
        %micro constants
        ripplesDistMicrChanForMerge = 15; %ms
        minNripplesForAnalysis_perChannel = 20; % dropping channels with low ripple rate
        
        %plotting constants
        nBinsHist = 10;
        maxLinesInFigureRipSpike = 4;
        maxColumnsInFigureDataMicro = 4;
        
        winFromLastSpike = 1000; %ms
        shortTimeRangeAfterStim = 3;%seconds
        midTimeRangeAfterStim = 60; %seconds
        stimulusDuration = 50; %ms
        minSpikeRateToIncludeUnit = 1;
        
        freqRangeForShowingSpindles = [5:30];
        xaxisForRipSp = [-500:500];
        specStartPointRipSp = 500;
        freqRangeSpRip = [50:150];
        nBinsPolar = 18;
    end
    
    methods
        
        function [rippleTimes, rippleStartEnd] = detectRipple(obj, data, sleepScoring, IIStimes, stim_times)
            
            % The method finds ripples based on Staresina et al (2015)
            % Input:
            % data - the raw data on which detection is required
            % sleepScoring - a vector in which each element represents the
            % sleep stage for an epoch of length obj.sleepEpochs (e.g. 1
            % second). The values which represent the required sleep stages
            % are kept in obj.sleepEpochs
            % Output:
            % rippleTimes - times of ripple peaks
            
            %convert sizes according to the sampling rate
            RMSWindowDuration = obj.RMSWindowDuration*obj.samplingRate/1000;
            minDurationAboveThresh = obj.minDurationAboveThresh*obj.samplingRate/1000;
            minDistBetweenRipples = obj.minDistBetweenRipples*obj.samplingRate/1000;
            
            if nargin < 3
                sleepScoring = [];
            end
            
            removeIIS = true;
            if nargin < 4 || isempty(IIStimes)
                removeIIS = false;
            end
            
            removeSTIM_artifacts = true;
            if nargin < 5 || isempty(stim_times)
                removeSTIM_artifacts = false;
            end
            
            
            %filter data to required range
            filteredData = obj.bandpass(data, obj.minFreq, obj.maxFreq);
            
            % if sleepScoring is nonempty: leave only the
            % segments in which there was sleep at the desired stage for
            % the analysis
            if ~isempty(sleepScoring)
                segLength = obj.scoringEpochDuration*obj.samplingRate;
                isSleep = zeros(1,length(sleepScoring)*segLength);
                for iEpoch = 1:length(sleepScoring)
                    if ismember(sleepScoring(iEpoch),obj.sleepEpochs)
                        isSleep((iEpoch-1)*segLength+1:iEpoch*segLength) = ones(1,segLength);
                    end
                end
                %match the length of data with the length of sleepScoring -
                %might get rid of some data points at the end if required, assuming it's
                %negligible
                if length(isSleep)>length(filteredData)
                    isSleep = isSleep(1:length(filteredData));
                else if length(isSleep)<length(filteredData)
                        filteredData = filteredData(1:length(isSleep));
                        data = data(1:length(isSleep));
                    end
                end
                %only leave segments of "real" sleep in data
                filteredData(~isSleep) = nan;
                data(~isSleep) = nan;
            end
            
            %remove windowAroundSTIM ms before and after every stimulation as
            %provided as input parameter
            if removeSTIM_artifacts
                winAroundIIS = obj.windowAroundSTIM*obj.samplingRate/1000;
                for iTime = 1:length(stim_times)
                    pointsBefore = min(stim_times(iTime),winAroundIIS);
                    pointsAfter = min(length(filteredData)-stim_times(iTime),winAroundIIS);
                    filteredData(stim_times(iTime)-pointsBefore+1:stim_times(iTime)+pointsAfter) = nan;
                    data(stim_times(iTime)-pointsBefore+1:stim_times(iTime)+pointsAfter) = nan;
                end
            end
            
            
            %remove windowAroundIIS ms before and after every IIS as
            %provided as input parameter
            if removeIIS
                winAroundIIS = obj.windowAroundIIS*obj.samplingRate/1000;
                for iTime = 1:length(IIStimes)
                    pointsBefore = min(IIStimes(iTime),winAroundIIS);
                    pointsAfter = min(length(filteredData)-IIStimes(iTime),winAroundIIS);
                    filteredData(IIStimes(iTime)-pointsBefore+1:IIStimes(iTime)+pointsAfter) = nan;
                    data(IIStimes(iTime)-pointsBefore+1:IIStimes(iTime)+pointsAfter) = nan;
                end
            end
            
            %calculate the root mean squared signal for windows of length RMSWindowDuration
            rmsSignal = zeros(1,length(filteredData)-RMSWindowDuration+1);
            for iPoint = 1:length(filteredData)-RMSWindowDuration+1
                rmsSignal(iPoint) = rms(filteredData(iPoint:iPoint+RMSWindowDuration-1));
            end
            
            %calculate the threshold as the rippleThreshPercentile
            %percentile of the rms signal
            rippleThresh = prctile(rmsSignal,obj.rippleThreshPercentile);
            
            %find windows that pass the thresh
            didPassThresh = rmsSignal>=rippleThresh;
            
            %find segments that pass the threshold for a duration longer than the threshold:
            %minDurationAboveThresh milliseconds
            rippleSegs = [];
            ind = 1;
            indRipple = 1;
            while ind <= length(didPassThresh)-minDurationAboveThresh+1
                if all(didPassThresh(ind:ind+minDurationAboveThresh-1))
                    rippleSegs(indRipple,1) = ind;
                    endSeg = ind+find(didPassThresh(ind:end)==0,1) - 2;
                    rippleSegs(indRipple,2) = endSeg + RMSWindowDuration - 1;
                    indRipple = indRipple+1;
                    ind = endSeg+2;
                else
                    ind = ind+1;
                end
            end
            
            %merge ripples who are close
            rippleSegsMerged = [];
            if ~isempty(rippleSegs)
                rippleDiffsSmall = (rippleSegs(2:end,1)-rippleSegs(1:end-1,2))<minDistBetweenRipples;
                
                indOld = 1;
                indNew = 0;
                
                while indOld < size(rippleSegs,1)
                    indNew = indNew+1;
                    if rippleDiffsSmall(indOld)==0
                        rippleSegsMerged(indNew,:) = rippleSegs(indOld,:);
                        indOld = indOld+1;
                    else
                        nextMerge = find(rippleDiffsSmall(indOld+1:end)==0,1)+indOld;
                        if isempty(nextMerge)
                            nextMerge = size(rippleSegs,1);
                        end
                        rippleSegsMerged(indNew,:) = [rippleSegs(indOld,1) rippleSegs(nextMerge,2)];
                        indOld = nextMerge+1;
                    end
                end
                if sum(rippleDiffsSmall)==0 || nextMerge<size(rippleSegs,1)
                    rippleSegsMerged(end+1,:) = rippleSegs(end,:);
                end
                
                rippleSegs = rippleSegsMerged;
                %go over ripples and leave only those with minNumOfExtreme
                %extremes in the raw data
                
                %smooth the data by averaging 3 consecutive points
                rippleTimes = [];
                rippleStartEnd = [];
                isRipple = false;
                for iRipple = 1:size(rippleSegs,1)
                    currRipple = data(rippleSegs(iRipple,1):rippleSegs(iRipple,2));
                    %average ripple
                    if isnan(currRipple)/length(currRipple)>=obj.minPercNaNAllowed
                        continue;
                    end
                    currRipple = movmean(currRipple, obj.NPointsToAverage);
                    localMax = islocalmax(currRipple);
                    if sum(localMax) >= obj.minNumOfExtreme
                        isRipple = true;
                    else
                        localMin = islocalmin(currRipple);
                        if sum(localMin) >= obj.minNumOfExtreme
                            isRipple = true;
                        end
                    end
                    
                    if isRipple
                        %find the location of the largest peak
                        currFilteredRipple = filteredData(rippleSegs(iRipple,1):rippleSegs(iRipple,2));
                        [~,maxLoc] = max(currFilteredRipple);
                        %                     localMaxInds = find(localMax);
                        %                     absMaxInd = localMaxInds(maxLoc)+rippleSegs(iRipple,1)-1;
                        absMaxInd = maxLoc+rippleSegs(iRipple,1)-1;
                        %add it to the ripple indices list
                        rippleTimes(end+1) = absMaxInd;
                        rippleStartEnd(end+1,:) = rippleSegs(iRipple,:);
                    end
                    
                    isRipple = false;
                end
            else
                
                rippleTimes = [];
                rippleStartEnd = [];
                
            end
        end
        
        %         function [rippleTimes, rippleStartEnd] = filterRipplesByRef(obj, rippleTimes, rippleStartEnd, refRipples)
        %         end
        
        function updateRippleTimes(obj, runData)
            %go over ripple times in folders and change max point - a one
            %time thing
            nPatients = length(runData);
            for iPatient = 1:nPatients
                filesInFolder = dir([runData(iPatient).RipplesFileNames,'*']);
                filesInFolder = {filesInFolder.name};
                for iFile = 1:length(filesInFolder)
                    disp(['updating ',filesInFolder{iFile}]);
                    
                    numRipple = str2num(filesInFolder{iFile}(12:end-4));
                    
                    try
                        data = [runData(iPatient).DataFolder '\CSC' num2str(numRipple) '.mat'];
                        data = load(data);
                        data = data.data;
                    catch
                        disp([runData(iPatient).DataFolder '\CSC' num2str(numRipple) '.mat doesn''t exist']);
                        continue;
                    end
                    
                    %filter data to required range
                    filteredData = obj.bandpass(data, obj.minFreq, obj.maxFreq);
                    
                    %load ripples
                    try
                        ripplesData = [runData(iPatient).RipplesFileNames num2str(numRipple) '.mat'];
                        ripplesData = load(ripplesData);
                        ripplesTimes = ripplesData.ripplesTimes;
                        ripplesStartEnd = ripplesData.ripplesStartEnd;
                    catch
                        disp([runData(iPatient).RipplesFileNames num2str(numRipple) '.mat doesn''t exist']);
                        ripplesTimes = [];
                        ripplesStartEnd = [];
                    end
                    
                    nRipples = length(ripplesTimes);
                    newRipplesTimes = zeros(1,nRipples);
                    for iRipple = 1:nRipples
                        [~,maxLoc] = max(filteredData(ripplesStartEnd(iRipple,1):ripplesStartEnd(iRipple,2)));
                        newRipplesTimes(iRipple) = maxLoc+ripplesStartEnd(iRipple,1)-1;
                    end
                    ripplesTimes = newRipplesTimes;
                    save([runData(iPatient).RipplesFileNames num2str(numRipple) '.mat'],'ripplesTimes','ripplesStartEnd');
                    disp(['re-saved ',filesInFolder{iFile}]);
                end
            end
        end
        
        function BP = bandpass(obj, timecourse, lowLimit, highLimit, filterOrder)
            
            
            %bandpass code - from Maya
            
            if (nargin < 5)
                filterOrder = obj.defaultFilterOrder;
            end
            
            % Maya GS - handle NAN values
            indices = find(isnan(timecourse));
            if length(indices) > obj.nanWarning*length(timecourse)
                warning('many NaN values in filtered signal')
            end
            timecourse(indices) = 0;
            %
            
            [b, a] = butter(filterOrder, [(lowLimit/obj.samplingRate)*2 (highLimit/obj.samplingRate)*2]);
            BP = filtfilt(b, a, timecourse );
            BP(indices) = NaN;
        end
        
        function plotRipples(obj, data, rippleTimes, folderToSave)
            
            %The method plots the ripples in rippleTimes with obj.secondBefAfter
            %before and after each ripple peak.
            %The subplot properties are set by class properties subplotSizeY,
            %subplotSizeX
            
            if nargin < 3 || isempty(folderToSave)
                toSave = false;
            else
                toSave = true;
            end
            
            %convert window size to sampling points
            secondBefAfter = obj.secondBefAfter*obj.samplingRate;
            
            %filter data to required range
            filteredData = obj.bandpass(data, obj.minFreq, obj.maxFreq);
            
            nRipples = length(rippleTimes);
            nInPlot = obj.subplotSizeX*obj.subplotSizeY;
            nPlots = ceil(nRipples/nInPlot);
            
            indRipple = 1;
            figInd = 1;
            for iPlot = 1:nPlots-1
                f = figure;
                
                for iInPlot = 1:nInPlot
                    subplot(obj.subplotSizeY,obj.subplotSizeX,iInPlot);
                    minInd = max(rippleTimes(indRipple)-secondBefAfter,1);
                    maxInd = min(rippleTimes(indRipple)+secondBefAfter,length(data));
                    
                    plot([-secondBefAfter:secondBefAfter]/obj.samplingRate,data(minInd:maxInd));
                    hold all;
                    plot([-secondBefAfter:secondBefAfter]/obj.samplingRate,filteredData(minInd:maxInd));
                    %                     xlim([1 secondBefAfter*2]);
                    
                    title(['Ripple time = ', num2str(rippleTimes(indRipple)/obj.samplingRate/60),' mins']);
                    indRipple = indRipple+1;
                end
                if toSave
                    set(f, 'Position', get(0, 'Screensize'));
                    saveas(f, [folderToSave,'\all_ripples_' num2str(figInd),'.jpg']);
                    close(f);
                else
                    pause;
                end
                figInd = figInd+1;
            end
            

        end
        
        
            end % methods
        
        
end
        