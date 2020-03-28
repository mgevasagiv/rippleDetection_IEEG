classdef RippleDetector < handle
    
    properties
        
        samplingRate = 1000;
        dataFilePrefix = 'CSC';
        
        %all the parameters are from Staresina et al (and Zhang et al is identical)
        minFreq = 80;
        maxFreq = 100;
        
        RMSWindowDuration = 20; %ms
        rippleThreshPercentile = 99;
        minDurationAboveThresh = 38; %ms
        minNumOfExtreme = 3; %min number of peaks / troughs in a ripple event
        NPointsToAverage = 3;
        minPercNaNAllowed = 0.1;
        minDistBetweenRipples = 20; %ms
        
        %IIS removal constants
        windowAroundIIS = 500; %ms
        
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
        
        firingRateWinSize = 100; %ms
        windowSpikeRateAroundRip = 500; %ms
        windowForSignificance = 100; %ms
        
        %params for spike rate around stimulations
        windowSpikeRateAroundStim = 500;
        controlDistForStim = 1000; %ms
        
        avgRippleBeforeAfter = 1; %second
        freqoiForAvgSpec = [0:0.5:10];
        freqRangeForAvgSpec = [1:200];
        timeBeforeAfterEventRipSpec = 1; %second
        timeForBaselineRip = 1; %second
        minNCycles = 5;
        minWinSizeSpec = 100; %ms
        
        %micro constants
        ripplesDistMicrChanForMerge = 15; %ms
        
        %plotting constants
        nBinsHist = 10;
        maxLinesInFigureRipSpike = 4;
        
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
        
        function [rippleTimes, rippleStartEnd] = detectRipple(obj, data, sleepScoring, IIStimes)
            
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
            
%             remPlots = rem(nRipples,nInPlot);
%             for iInPlot = 1:remPlots
%                 subplot(obj.subplotSizeY,obj.subplotSizeX,iInPlot);
%                 minInd = max(rippleTimes(indRipple)-secondBefAfter,1);
%                 maxInd = min(rippleTimes(indRipple)+secondBefAfter,length(data));
%                 hold off;
%                 plot(data(minInd:maxInd));
%                 hold all;
%                 plot(filteredData(minInd:maxInd));
%                 title(['Ripple time = ', num2str(rippleTimes(indRipple))]);
%                 indRipple = indRipple+1;
%             end
        end
        
        function [fireRateRip, fireRateControl] = checkSpikeRateAtRip(obj, rippleTimes, rippleStartEnd, spikeTimes, dataDuration)
            
            % a method which receives the ripple and spike times and returns the spike rate around ripples and around 
            % controls
            % The input is:
            % rippleTimes – array of ripple times
            % rippleStartEnd – matrix where each row is a ripple and index <I,1> is the start point of ripple I and index 
            % <I,2> is the end point of ripple i.
            % spikeTimes – array of spike times, can also be a cell array of arrays of spike times, in which case the 
            % analysis will be performed on the average spike rate.
            % dataDuration (optional) – the length of the session in ms.
            % The output is:
            % fireRateRip – matrix with the number of rows and the number of ripples where each row is the average spike 
            % rate around the ripples.
            % fireRateControl - matrix with the number of rows and the number of ripples where each row is the average 
            % spike rate around the controls.
            
            
            if nargin < 5
                dataDuration = max(rippleTimes) + obj.windowSpikeRateAroundRip+obj.maxDistControl;
            end
            
            %possible - add a check whether control falls on IIS (not very
            %important as control is very close to ripple and ripples are
            %not detected in the vicinity of IIS). same goes for sleep
            %scoring
            %             removeIIS = true;
            %             if nargin < 4 || isempty(IIStimes)
            %                 removeIIS = false;
            %             end
            %
            
            nRipples = size(rippleStartEnd,1);
            %convert indices to ms
            rippleStartEnd = round(rippleStartEnd*1000/obj.samplingRate);
            rippleTimes = round(rippleTimes*1000/obj.samplingRate);
            
            %mark all ripple segments (for choosing control points later)
            ripplesSegsLog = zeros(1,dataDuration);
            for iRipple = 1:nRipples
                ripplesSegsLog(rippleStartEnd(iRipple,1):rippleStartEnd(iRipple,2)) = 1;
            end
            
            if ~iscell(spikeTimes)
                spikeTimes = {spikeTimes};
            end
            
            nUnits = length(spikeTimes);
            
            spikes = cell(1,nUnits);
            spikeRate = cell(1,nUnits);
            
            %convert spike times to smoothed spike rate function
            for iUnit = 1:nUnits
                spikes{iUnit} = zeros(1, dataDuration);
                spikes{iUnit}(round(spikeTimes{iUnit})) = 1;
                spikeRate{iUnit} = movsum(spikes{iUnit},obj.firingRateWinSize,'Endpoints','fill')/(obj.firingRateWinSize/1000);
            end
            
            fireRateRip = [];
            fireRateControl = [];
            
            for iRipple=1:nRipples

                %the control is a random point in the vicinity (as determined by min/maxDistControl) of the ripple
                %which doesn't fall on another ripple
                
                possInds = [rippleStartEnd(iRipple,1)-obj.maxDistControl:rippleStartEnd(iRipple,1)-obj.minDistControl rippleStartEnd(iRipple,2)+obj.minDistControl:rippleStartEnd(iRipple,2)+obj.maxDistControl];
                indsWithNoRipples = find(ripplesSegsLog(possInds)==0);
                %if after taking into account other ripples no indices are
                %available for control, enlarge the window (theoretically
                %this could fail also, but not likely)
                if isempty(indsWithNoRipples)
                    possInds = [rippleStartEnd(iRipple,1)-2*obj.maxDistControl:rippleStartEnd(iRipple,1)-obj.minDistControl rippleStartEnd(iRipple,2)+obj.minDistControl:rippleStartEnd(iRipple,2)+obj.maxDistControl*2];
                    indsWithNoRipples = find(ripplesSegsLog(possInds)==0);
                end
                controlIndForRate = possInds(indsWithNoRipples(randi(length(indsWithNoRipples))));
                     
                %add to spike rate averaged around ripple 
                if rippleTimes(iRipple)-obj.windowSpikeRateAroundRip>=1 && rippleTimes(iRipple)+obj.windowSpikeRateAroundRip<=dataDuration
                    currUnitsAvg = zeros(1,2*obj.windowSpikeRateAroundRip+1);
                    currUnitsAvgControl = zeros(1,2*obj.windowSpikeRateAroundRip+1);
                    for iUnit = 1:nUnits
                        currUnitsAvg = currUnitsAvg+spikeRate{iUnit}(rippleTimes(iRipple)-obj.windowSpikeRateAroundRip:rippleTimes(iRipple)+obj.windowSpikeRateAroundRip);
                        currUnitsAvgControl = currUnitsAvgControl+spikeRate{iUnit}(controlIndForRate-obj.windowSpikeRateAroundRip:controlIndForRate+obj.windowSpikeRateAroundRip);
                    end
                    fireRateRip(end+1,:) = currUnitsAvg/nUnits;
                    fireRateControl(end+1,:) = currUnitsAvgControl/nUnits;
                end

            end
            
        end
        
        function [spikeRateSession, rateAroundStim, rateAroundControl] = getSpikeRateAtStim(obj, stimTimes, spikeTimes, dataDuration)
            
            % receives the spike times and the stimulation times and calculates spike rates around stimulations and around 
            % control, and the spike rate through the session.
            % The controls for stimulations are the time points controlDistForStim (by default 1000ms) before the 
            % stimulations.
            % The spike rates are calculated by using movsum with a window size firingRateWinSize (by default 100 ms).
            % Input:
            % stimTimes – an array of stimulation times
            % spikeTimes – an array of spike times
            % dataDuration (optional) – length of session (ms)
            % output:
            % spikeRateSession – the spike rate function through the entire session
            % rateAroundStim – a matrix where each row is the spike rate around a stimulation
            % rateAroundControl – a matrix where each row is the spike rate around the control

            if nargin < 4
                dataDuration = max(stimTimes) + obj.shortTimeRangeAfterStim*obj.samplingRate;
            end
            
            %if data duration of interest is shorter than spike times
            spikeTimes = spikeTimes(spikeTimes<=dataDuration);
            
            nStim = length(stimTimes);
            
            %convert spike times to a spike rate function
            spikes = zeros(1, dataDuration);
            spikes(round(spikeTimes)) = 1;
            spikeRateSession = movsum(spikes,obj.firingRateWinSize,'Endpoints','fill')/(obj.firingRateWinSize/1000);
                        
            rateAroundStim = zeros(nStim,obj.windowSpikeRateAroundStim*2+1);
            rateAroundControl = zeros(nStim,obj.windowSpikeRateAroundStim*2+1);
            
            %go over the stimulations and build the rate around
            %stimulations and rate around control matrices
            for iStim = 1:nStim
                
                currStimTime = stimTimes(iStim);
                currControlTime = currStimTime-obj.controlDistForStim;
                
                rateAroundStim(iStim,:) = spikeRateSession(currStimTime-obj.windowSpikeRateAroundStim:currStimTime+obj.windowSpikeRateAroundStim);
                rateAroundControl(iStim,:) = spikeRateSession(currControlTime-obj.windowSpikeRateAroundStim:currControlTime+obj.windowSpikeRateAroundStim);
                
            end
            
            %the rate function through the entire session (by default up
            %until the last stimulation)
            spikeRateSession = movsum(spikes,10000,'Endpoints','fill')/(10000/1000);
            
        end
        
        function results = runRippleData(obj, runData, fileNameResults)
            
            % The method produces information about the ripples in a channel that includes:
            % A. Average ripple – before stimulation and during stimulations (short effect).
            % B. Spectrum of average ripple – before stimulation and during stimulations (short effect).
            % C. Average of TFR around ripples – before stimulation and during stimulations (short effect).
            % D. Average of TFR around spindles – before stimatulation.
            % E. Polar histogram of synchronization index between spindles and ripples range.
            %
            % The input runData is a struct in the length of number of patients (for which the analysis is required). 
            % In addition it receives the input parameter fileNameResults which includes the file name into which the results 
            % will be saved (optional).
            % Each element (=patient) in runData should include the fields:
            % patientName
            % channelsToRunOn – list of channel indices for which to perform the analysis.
            % DataFolder – The folder in which the raw data files are saved (the method assumes the prefix for the files is 
            % CSC, can be changed by the property dataFilePrefix).
            % macroMontageFileName - the file name (including path) of the macromontage.
            % RipplesFileNames - name (including path) of the ripple mat files in which the ripple times for the macro 
            % channels are saved (the method assumes the name of the file is RipplesFileNames <#channel index>
            % SpindlesFileNames - name (including path) of the spindle mat files in which the spindle times for the macro 
            % channels are saved (the method assumes the name of the file is SpindlesFileNames <#channel index>
            % SpikesFileNames - name (including path) of the spikes mat files in which the spikes times for the macro channels 
            % are saved (the method assumes the name of the file is SpikesFileNames <#channel index>). If not provided spikes 
            % will not be removed from the data for the analysis.
            % sleepScoringFileName – file name (including path) of the sleep scoring mat file. If not provided all the data will be used.
            %
            % The output struct results includes all the results of the analysis, which can then be plotted using 
            % plotResultsRipplesData. The output struct is a struct with the length of the number of patients (=the length 
            % of runData), where each element includes:
            % patientName
            % resultsPerChan – a struct in the length of the number of channels required for the analysis per the patient. 
            % Each element in resultsPerChan includes the fields:
            % channelNum
            % area
            % nRipplesBefore, nRipplesStim - number of ripples before stimulations and after stimulations (short effect) 
            % respectively
            % avgBefore, avgStim – average ripples before stimulations and after stimulations (short effect) respectively
            % stdBefore, stdStim – std of ripples before stimulations and after stimulations (short effect) respectively
            % specBefore, specStim – spectrum of ripples average before stimulations and after stimulations (short effect) 
            % respectively
            % meanTFRRipBefore, meanTFRRipStim – mean of ripple triggered TFR before stimulations and after stimulations 
            % (short effect) respectively
            % SIanglesSpRip – Synchronization Indices of spindles-ripples before stimaulations (an array with the length as 
            % number of spindles)
            % R – results of the r-test for the polar histogram (of SIanglesSpRip)
            % V – results of the v-test for the polar histogram (of SIanglesSpRip)
            % meanSpecs – mean spindle-triggered TFR before stimulations.
            % meanEpochs – mean spindle for the spindles for which meanSpecs were calculated.
            % nEpochs – number of spindles in the spindle-ripple analyses (all before stimatulions).

            
            if nargin < 3
                fileNameResults = '';
            end
            
            shortTimeRangeAfterStim = obj.shortTimeRangeAfterStim*obj.samplingRate;
            midTimeRangeAfterStim = obj.midTimeRangeAfterStim*obj.samplingRate;
            stimulusDuration = obj.stimulusDuration*obj.samplingRate/1000;
            avgRippleBeforeAfter = obj.avgRippleBeforeAfter*obj.samplingRate;
            timeWin = min((1./obj.freqRangeForAvgSpec)*obj.minNCycles,ones(size(obj.freqRangeForAvgSpec))*obj.minWinSizeSpec);
            
            nPatients = length(runData);
            
            %go over all required patients
            for iPatient = 1:nPatients
                disp(['Patient ',runData(iPatient).patientName,' ',num2str(iPatient),'/',num2str(nPatients)]);
                results(iPatient).patientName = runData(iPatient).patientName;
                
                %load exp data for stimulation timings
                expData = load(runData(iPatient).ExpDataFileName);
                expData = expData.EXP_DATA;
                stimTimes = expData.stimTiming.validatedTTL_NLX;
                firstStim = stimTimes(1);
                
                %load sleep scoring
                if isfield(runData(iPatient), 'sleepScoringFileName') && ~isempty(runData(iPatient).sleepScoringFileName)
                    try
                        sleepScoring = load(runData(iPatient).sleepScoringFileName);
                        sleepScoring = sleepScoring.sleep_score_vec;
                    catch
                        disp([runData(iPatient).sleepScoringFileName '.mat doesn''t exist']);
                        sleepScoring = [];
                    end
                end
                
                %load macro montage
                macroMontage = load(runData(iPatient).macroMontageFileName);
                macroMontage = macroMontage.MacroMontage;
                
                %go over required channels per patient
                resultsPerChan = [];
                nChans = length(runData(iPatient).channelsToRunOn);
                for iChan = 1:nChans
                    currChan = runData(iPatient).channelsToRunOn(iChan);
                    resultsPerChan(iChan).channelNum = currChan;
                    
                    disp(['channel ',num2str(currChan),' ',num2str(iChan)','/',num2str(nChans)]);
                    
                    currArea = macroMontage(resultsPerChan(iChan).channelNum).Area;
                    resultsPerChan(iChan).area = currArea;
                    
                    %load the data
                    try
                        data = [runData(iPatient).DataFolder '\CSC' num2str(currChan) '.mat'];
                        data = load(data);
                        data = data.data;
                    catch
                        disp([runData(iPatient).DataFolder '\CSC' num2str(currChan) '.mat doesn''t exist']);
                        continue;
                    end
                    
                    %load IIS times
                    if isfield(runData(iPatient), 'SpikesFileNames') && ~isempty(runData(iPatient).SpikesFileNames)
                        try
                            IIStimes = [runData(iPatient).SpikesFileNames num2str(currChan) '.mat'];
                            IIStimes = load(IIStimes);
                            IIStimes = IIStimes.peakTimes;
                        catch
                            disp([runData(iPatient).SpikesFileNames num2str(currChan) '.mat doesn''t exist']);
                            IIStimes = [];
                        end
                    else
                        IIStimes = [];
                    end
                    
                    %load ripples
                    if isfield(runData(iPatient), 'RipplesFileNames') && ~isempty(runData(iPatient).RipplesFileNames)
                        try
                            ripplesData = [runData(iPatient).RipplesFileNames num2str(currChan) '.mat'];
                            ripplesData = load(ripplesData);
                            ripplesTimes = ripplesData.ripplesTimes;
                        catch
                            disp([runData(iPatient).RipplesFileNames num2str(currChan) '.mat doesn''t exist']);
                            ripplesTimes = [];
                        end
                    else
                        ripplesTimes = [];
                    end
                    
                    %load spindles
                    if isfield(runData(iPatient), 'SpindlesFileNames') && ~isempty(runData(iPatient).SpindlesFileNames)
                        try
                            spindlesTimes = [runData(iPatient).SpindlesFileNames num2str(currChan) '.mat'];
                            spindlesTimes = load(spindlesTimes);
                            spindlesTimes = spindlesTimes.spindlesTimes;
                        catch
                            disp([runData(iPatient).SpindlesFileNames num2str(currChan) '.mat doesn''t exist']);
                            spindlesTimes = [];
                        end
                    end
                    
                    %indices of ripples before stimulations
                    ripplesBeforeInds = ripplesTimes<firstStim-obj.windowSpikeRateAroundRip;
                    %get inds of ripples that are short and mid effect and not too
                    %close to the stimulus
                    dataDuration = stimTimes(end)+midTimeRangeAfterStim;
                    stimInds = zeros(1,dataDuration);
                    %short effect
                    for iStim = 1:length(stimTimes)
                        stimInds(stimTimes(iStim)+stimulusDuration:stimTimes(iStim)+shortTimeRangeAfterStim) = 1;
                    end
                    %add also mid effect
                    stimDiffs = diff(stimTimes);
                    stimIndsWithMidPauseAfter = [find(stimDiffs >= midTimeRangeAfterStim) length(stimTimes)];
                    stimIndsWithMidPauseAfterTimes = stimTimes(stimIndsWithMidPauseAfter);
                    for iStim = 1:length(stimIndsWithMidPauseAfter)
                        stimInds(stimIndsWithMidPauseAfterTimes(iStim)+shortTimeRangeAfterStim:stimIndsWithMidPauseAfterTimes(iStim)+midTimeRangeAfterStim) = 1;
                    end
                    ripplesIndsLog = zeros(1,dataDuration);
                    ripplesIndsLog(ripplesTimes(ripplesTimes<=dataDuration)) = 1;
                    %inds of short and mid effect ripples
                    ripplesDuringStimInds = ismember(ripplesTimes,find(ripplesIndsLog & stimInds));
                    
                    
                    %calculate average ripples and TFRs
                    ripplesTimesBefore = ripplesTimes(ripplesBeforeInds);
                    nRipplesBefore = length(ripplesTimesBefore);
                    
                    ripplesTimesStim = ripplesTimes(ripplesDuringStimInds);
                    nRipplesStim = length(ripplesTimesStim);
                    
                    %calculate average of ripples
                    ripplesBefore = zeros(nRipplesBefore,length([-avgRippleBeforeAfter:avgRippleBeforeAfter]));
                    ripplesStim = zeros(nRipplesStim,length([-avgRippleBeforeAfter:avgRippleBeforeAfter]));
                    
                    %before stimulations
                    for iRipple = 1:nRipplesBefore
                        ripplesBefore(iRipple,:) = data(ripplesTimesBefore(iRipple)-avgRippleBeforeAfter:ripplesTimesBefore(iRipple)+avgRippleBeforeAfter);
                    end
                    if ~isempty(ripplesBefore)
                        avgBefore = nanmean(ripplesBefore);
                        stdBefore = nanstd(ripplesBefore);
                        %spectrum of average
                        specBefore = ft_specest_mtmfft(avgBefore,[1:length(avgBefore)]/obj.samplingRate,'freqoi',obj.freqoiForAvgSpec,'taper','hanning');
                        specBefore = abs(squeeze(specBefore(1,1,:,:)));
                    else
                        avgBefore = nan(1,size(ripplesBefore,2));
                        stdBefore = nan(1,size(ripplesBefore,2));
                        specBefore = nan(length(obj.freqoiForAvgSpec),1);
                    end
                        
                    %after stimulations
                    for iRipple = 1:nRipplesStim
                        ripplesStim(iRipple,:) = data(ripplesTimesStim(iRipple)-avgRippleBeforeAfter:ripplesTimesStim(iRipple)+avgRippleBeforeAfter);
                    end
                    if ~isempty(ripplesStim)
                        avgStim = nanmean(ripplesStim);
                        stdStim = nanstd(ripplesStim);
                        specStim = ft_specest_mtmfft(avgStim,[1:length(avgStim)]/obj.samplingRate,'freqoi',obj.freqoiForAvgSpec,'taper','hanning');
                        specStim = abs(squeeze(specStim(1,1,:,:)));
                    else
                        avgStim = nan(1,size(ripplesStim,2));
                        stdStim = nan(1,size(ripplesStim,2));
                        specStim = nan(1,length(obj.freqoiForAvgSpec));
                    end
                    
                    %get ripple centered TFRs (implemented in
                    %PACCalculator)
                    pacCalc = PACCalculator;
                    pacCalc.freqRange = obj.freqRangeForAvgSpec;
                    pacCalc.timeBeforeAfterEvent = obj.timeBeforeAfterEventRipSpec; %seconds
                    pacCalc.timeForBaseline = obj.timeForBaselineRip; %seconds, from Starestina et al
                    pacCalc.minNCycles = obj.minNCycles;
                    
                    meanTFRRipBefore = pacCalc.plotAvgSpecDiff(data, ripplesTimesBefore);
                    meanTFRRipStim = pacCalc.plotAvgSpecDiff(data, ripplesTimesStim);
                    
                    resultsPerChan(iChan).nRipplesBefore = nRipplesBefore;
                    resultsPerChan(iChan).nRipplesStim = nRipplesStim;
                    
                    resultsPerChan(iChan).avgBefore = avgBefore;
                    resultsPerChan(iChan).avgStim = avgStim;
                    resultsPerChan(iChan).stdBefore = stdBefore;
                    resultsPerChan(iChan).stdStim = stdStim;
                    resultsPerChan(iChan).specBefore = specBefore;
                    resultsPerChan(iChan).specStim = specStim;
                    resultsPerChan(iChan).meanTFRRipBefore = meanTFRRipBefore;
                    resultsPerChan(iChan).meanTFRRipStim = meanTFRRipStim;
                    
                    %get spindle-ripple coupling (SI angles and spectrum) for before the stimluations -
                    %implemented in AnalyzeCoupling
                    ac = AnalyzeCoupling;
                    spindlesBefore = spindlesTimes(spindlesTimes<firstStim-ac.timeBeforeAfterEventSpRip*obj.samplingRate);
                    detections.slowWavesTimes = [];
                    detections.spindlesTimes = spindlesBefore;
                    detections.ripplesTimes = [];
                    ac.freqRangeSpRip = obj.freqRangeSpRip;
                    [~, SIanglesSpRip, rs, vs, meanSpecs, meanEpochs ~, nEpochs] = ac.analyzeUsingSpectrumsAndSIIndex(data, data, {IIStimes,IIStimes}, sleepScoring, detections, [], true, false, true);
                    
                    resultsPerChan(iChan).SIanglesSpRip = SIanglesSpRip;
                    resultsPerChan(iChan).r = rs(2);
                    resultsPerChan(iChan).v = vs(2);
                    resultsPerChan(iChan).meanSpecs = meanSpecs{2};
                    resultsPerChan(iChan).meanEpochs = meanEpochs{2};
                    resultsPerChan(iChan).nEpochs = nEpochs{2};
                    
                end
                
                results(iPatient).resultsPerChan = resultsPerChan;
            end
            
            if ~isempty(fileNameResults)
                save(fileNameResults,'results');
            end
            
        end
        
        function results = runRipSpikes(obj, runData, fileNameResults)
            % The method produces information about correlation between ripples and spikes.
            % The analysis is performed per requested macro channel and examines the relationship between ripples in that 
            % macro channel and spikes in all the micro channels in the same area (based on the macro and micro montages). 
            % Analysis is performed for the macro channel vs spike rate in each single/multi unit in the area and also vs 
            % the average spike rate of the single/multi units. The property spikeMultiUnits sets whether the analysis is 
            % on single or multi units (by default it’s set to true – i.e. multi units).
            % The main result is the firing rate in a time window around ripples and around controls for ripples before 
            % the stimulations and for ripples during the stimulations (short+mid effect). The time window is set by 
            % windowSpikeRateAroundRip. The controls are random points which are: a. at least minDistControl 
            % before or after the ripple, b. no more than maxDistControl before or after the ripple, c. do not coincide 
            % with a different ripple.
            % Single/multi units in which the spike rate is below minSpikeRateToIncludeUnit (by default 1Hz) are not 
            % included in the analysis.
            % The method also returns stimulation triggered spike rate and the spike rate through the entire session 
            % (before and during stimulations). The stimulation triggered spike rate is presented vs a control, which is 
            % the time point controlDistForStim before the stimulation.
            % The spike rates are calculated by using movsum with a window size firingRateWinSize.
            %
            % The input runData is a struct in the length of number of patients (for which the analysis is required). 
            % In addition it receives the input parameter fileNameResults which includes the file name into which the 
            % results will be saved (optional).
            % Each element (=patient) in runData should include the fields:
            % patientName
            % channelsToRunOn – list of channel indices for which to perform the analysis.
            % ExpDataFileName – name (including path) of the EXP_DATA for the patient.
            % spikeData – name (including path) of the file which includes the spike data for the patient (i.e. the file 
            % which usually has the name <patientName>_spike_timestamps_post_processing.mat).
            % macroMontageFileName - the file name (including path) of the macromontage.
            % RipplesFileNames - name (including path) of the ripple mat files in which the ripple times for the macro 
            % channels are saved (the method assumes the name of the file is RipplesFileNames <#channel index>
            %
            % The output struct results includes all the results of the analysis, which can then be plotted using 
            % plotResultsSpikes. The output struct is a struct with the length of the number of patients (=the length of runData), where each element includes:
            % patientName
            % resultsPerChan – a struct in the length of the number of channels required for the analysis per the patient. 
            % Each element in resultsPerChan includes the fields:
            % channelNum
            % area
            % unitInds – a cell array with the single/multi unit indices for this area. If using multi units each element 
            % in the cell array is an array including all the single units per multi unit.
            % spikeTimes – a cell array with the length as the number of units where each element is the spike times for 
            % that unit.
            % fireRateRipBefore, fireRateRipStim – cell array with the length as the number of units where each element is 
            % the average spike rate around ripples for that unit for all the ripples before the stimulation and during the 
            % stimulations (short+mid effect together) respectively.
            % fireRateRipAvgUnitsBefore, fireRateRipAvgUnitsStim – the average spike rate around ripples averaged over all 
            % units for all the ripples before the stimulation and during the stimulations (short+mid effect together) 
            % respectively.
            % fireRateControlBefore, fireRateControlStim - cell array with the length as the number of units where each 
            % element is the average spike rate around controls for that unit for all the ripples before the stimulation 
            % and during the stimulations (short+mid effect together) respectively.
            % fireRateControlAvgUnitsBefore, fireRateControlAvgUnitsStim - the average spike rate around controls for all 
            % the ripples averaged over all units before the stimulation and during the stimulations (short+mid effect 
            % together) respectively.
            % allSessionFireRates - cell array with the length as the number of units where each element is spike rate 
            % across the session (before and during stimulations).
            % stimTriggeredFireRates - cell array with the length as the number of units where each element is a matrix 
            % where each row is the spike rate function around a stimulation.
            % controlForStimTriggered - cell array with the length as the number of units where each element is a matrix 
            % where each row is the spike rate function around the control per stimulation.
            % firstStim – time of first stimulation.
            % lastStim – time of last stimulation.

            if nargin < 3
                fileNameResults = '';
            end
            
            shortTimeRangeAfterStim = obj.shortTimeRangeAfterStim*obj.samplingRate;
            midTimeRangeAfterStim = obj.midTimeRangeAfterStim*obj.samplingRate;
            stimulusDuration = obj.stimulusDuration*obj.samplingRate/1000;
            
            %go over the patients
            nPatients = length(runData);
            for iPatient = 1:nPatients
                disp(['Patient ',runData(iPatient).patientName]);
                results(iPatient).patientName = runData(iPatient).patientName;
                
                %load exp data for stimulation timings
                expData = load(runData(iPatient).ExpDataFileName);
                expData = expData.EXP_DATA;
                stimTimes = expData.stimTiming.validatedTTL_NLX;
                firstStim = stimTimes(1);
                lastStim = stimTimes(end);
                
                %load macro montage
                macroMontage = load(runData(iPatient).macroMontageFileName);
                macroMontage = macroMontage.MacroMontage;
                
                %load spike data
                try
                    spikeData = load(runData(iPatient).spikeData);
                catch
                    disp([runData(iPatient).spikeData ' doesn''t exist, continuing']);
                    continue;
                end
                
                %go over channels for analysis
                resultsPerChan = [];
                nChans = length(runData(iPatient).channelsToRunOn);
                for iChan = 1:nChans
                    currChan = runData(iPatient).channelsToRunOn(iChan);
                    resultsPerChan(iChan).channelNum = currChan;
                    
                    disp(['channel ',num2str(currChan),' ',num2str(iChan)','/',num2str(nChans)]);
                    
                    currArea = macroMontage(resultsPerChan(iChan).channelNum).Area;
                    resultsPerChan(iChan).area = currArea;
                    
                    %load the ripples per channel
                    if isfield(runData(iPatient), 'RipplesFileNames') && ~isempty(runData(iPatient).RipplesFileNames)
                        try
                            ripplesData = [runData(iPatient).RipplesFileNames num2str(currChan) '.mat'];
                            ripplesData = load(ripplesData);
                            ripplesTimes = ripplesData.ripplesTimes;
                            ripplesStartEnd = ripplesData.ripplesStartEnd;
                        catch
                            disp([runData(iPatient).RipplesFileNames num2str(currChan) '.mat doesn''t exist']);
                            ripplesTimes = [];
                            ripplesStartEnd = [];
                        end
                    else
                        ripplesTimes = [];
                        ripplesStartEnd = [];
                    end
                    
                    %option to add ripple detection - requires data and
                    %IIStimes which are currently not loaded
                    %                     if isempty(rippleTimes)
                    %                         [ripplesTimes, ripplesStartEnd] = obj.detectRipple(data, sleepScoring, IIStimes);
                    %                     end
                    
                    %ripple times before the stimulations
                    ripplesBeforeStimInds = ripplesTimes<(firstStim-obj.windowSpikeRateAroundRip);

                    %ripple times during the stimulations - only ripples that are short and mid effect and not too
                    %close to the stimulus are counted
                    dataDuration = stimTimes(end)+midTimeRangeAfterStim;
                    stimInds = zeros(1,dataDuration);
                    %short effect
                    for iStim = 1:length(stimTimes)
                        stimInds(stimTimes(iStim)+stimulusDuration:stimTimes(iStim)+shortTimeRangeAfterStim) = 1;
                    end
                    %add mid effect
                    stimDiffs = diff(stimTimes);
                    stimIndsWithMidPauseAfter = [find(stimDiffs >= midTimeRangeAfterStim) length(stimTimes)];
                    stimIndsWithMidPauseAfterTimes = stimTimes(stimIndsWithMidPauseAfter);
                    for iStim = 1:length(stimIndsWithMidPauseAfter)
                        stimInds(stimIndsWithMidPauseAfterTimes(iStim)+shortTimeRangeAfterStim:stimIndsWithMidPauseAfterTimes(iStim)+midTimeRangeAfterStim) = 1;
                    end
                    ripplesIndsLog = zeros(1,dataDuration);
                    ripplesIndsLog(ripplesTimes(ripplesTimes<=dataDuration)) = 1;
                    ripplesDuringStimInds = ismember(ripplesTimes,find(ripplesIndsLog & stimInds));
                    
                    %find all units recorded in the area using the spike
                    %data
                    try
                        unitInds = find(strcmp(extractfield(spikeData.micro_channels_spike_summary.unit_list_xls,'Location'),currArea));
                    catch
                        unitInds = [];
                    end
                    
                    %if obj.spikeMultiUnits is true regard multi units for
                    %the analysis, otherwise regard single units
                    if obj.spikeMultiUnits
                        if ~isempty(unitInds)
                            currChannels = [spikeData.micro_channels_spike_summary.unit_list_xls(unitInds).Channel];
                        else
                            currChannels = [];
                        end
                        %find all the unique channels for the single units
                        %(i.e. multiunits channels)
                        uChannels = unique(currChannels);
                        
                        %a cell array that will include unit indices - each
                        %element has the indices for the multiunit in the
                        %area
                        newUnitInds = {};
                        
                        spikeTimes = {};
                        
                        for iChannel = 1:length(uChannels)
                            currUinds = unitInds(currChannels==uChannels(iChannel));
                            %merge spike data to create multi unit
                            currSpikeTimes = sort(cat(1,spikeData.micro_channels_spike_summary.spike_timestamps{currUinds}));
                            currSpikeTimes = currSpikeTimes(currSpikeTimes <= dataDuration);
                            %only use the multiunit if the firing rate is above  a threshold 
                            if length(currSpikeTimes)/(dataDuration/obj.samplingRate) >= obj.minSpikeRateToIncludeUnit
                                spikeTimes{end+1} = currSpikeTimes;
                                newUnitInds{end+1} = currUinds;
                            end
                        end
                        
                        unitInds = newUnitInds;
                        nUnits = length(unitInds);
                    else %single units
                        %remove units with too low spike rate
                        newUnitInds = {};
                        spikeTimes = {};
                        for iUnit = 1:length(unitInds)
                            currSpikeTimes = spikeData.micro_channels_spike_summary.spike_timestamps{unitInds(iUnit)};
                            currSpikeTimes = currSpikeTimes(currSpikeTimes <= dataDuration);
                            if length(currSpikeTimes)/(dataDuration/obj.samplingRate) >= obj.minSpikeRateToIncludeUnit
                                newUnitInds{end+1} = unitInds(iUnit);
                                spikeTimes{end+1} = currSpikeTimes;
                            end
                        end
                        unitInds = newUnitInds;
                        nUnits = length(unitInds);
                    end
                    
                    %for each unit - find the spike rate for ripples and controls before
                    %and during the stimulations
                    fireRateRipBefore = cell(1,nUnits);
                    fireRateControlBefore = cell(1,nUnits);
                    
                    fireRateRipStim = cell(1,nUnits);
                    fireRateControlStim = cell(1,nUnits);
                    
                    %and the stimulations triggered spike rate and their
                    %controls, and the fire rate through the entire session
                    stimTriggeredFireRates = cell(1,nUnits);
                    controlForStimTriggered = cell(1,nUnits);
                    allSessionFireRates = cell(1,nUnits);
                    
                    %go over the units (single/multi)
                    for iUnit = 1:nUnits
                        [fireRateRipBefore{iUnit}, fireRateControlBefore{iUnit}] = obj.checkSpikeRateAtRip(ripplesTimes(ripplesBeforeStimInds), ripplesStartEnd(ripplesBeforeStimInds,:), spikeTimes{iUnit});
                        [fireRateRipStim{iUnit}, fireRateControlStim{iUnit}] = obj.checkSpikeRateAtRip(ripplesTimes(ripplesDuringStimInds), ripplesStartEnd(ripplesDuringStimInds,:), spikeTimes{iUnit});
                        
                        [allSessionFireRates{iUnit}, stimTriggeredFireRates{iUnit}, controlForStimTriggered{iUnit}] = obj.getSpikeRateAtStim(stimTimes, spikeTimes{iUnit});
                    end
                                        
                    %find the spike rate for ripples and controls before
                    %and during the stimulations for the average spike rate
                    %over the units
                    [fireRateRipAvgUnitsBefore, fireRateControlAvgUnitsBefore] = obj.checkSpikeRateAtRip(ripplesTimes(ripplesBeforeStimInds), ripplesStartEnd(ripplesBeforeStimInds,:), spikeTimes);
                    [fireRateRipAvgUnitsStim, fireRateControlAvgUnitsStim] = obj.checkSpikeRateAtRip(ripplesTimes(ripplesDuringStimInds), ripplesStartEnd(ripplesDuringStimInds,:), spikeTimes);
                    
                    resultsPerChan(iChan).unitInds = unitInds;
                    resultsPerChan(iChan).spikeTimes = spikeTimes;
                    
                    resultsPerChan(iChan).fireRateRipBefore = fireRateRipBefore;
                    resultsPerChan(iChan).fireRateRipAvgUnitsBefore = fireRateRipAvgUnitsBefore;
                    resultsPerChan(iChan).fireRateControlBefore = fireRateControlBefore;
                    resultsPerChan(iChan).fireRateControlAvgUnitsBefore = fireRateControlAvgUnitsBefore;
                    
                    resultsPerChan(iChan).fireRateRipStim = fireRateRipStim;
                    resultsPerChan(iChan).fireRateRipAvgUnitsStim = fireRateRipAvgUnitsStim;
                    resultsPerChan(iChan).fireRateControlStim = fireRateControlStim;
                    resultsPerChan(iChan).fireRateControlAvgUnitsStim = fireRateControlAvgUnitsStim;
                    
                    resultsPerChan(iChan).allSessionFireRates = allSessionFireRates;
                    resultsPerChan(iChan).stimTriggeredFireRates = stimTriggeredFireRates;
                    resultsPerChan(iChan).controlForStimTriggered = controlForStimTriggered;
                    resultsPerChan(iChan).firstStim = firstStim;
                    resultsPerChan(iChan).lastStim = lastStim;
                    
                end
                
                results(iPatient).resultsPerChan = resultsPerChan;
            end
            
            if ~isempty(fileNameResults)
                save(fileNameResults,'results');
            end
            
        end
        
        function [rippleTimesMerged,rippleStartEndMerged] = getRipplesFromAllMicroElectrodesInArea(obj, ripplesTimes, ripplesStartEnd)
            % The ripples detection on micro channels is performed per area – i.e. the detection takes into account all the 
            % channels in that area. The method getRipplesFromMicroElectrodesInArea receives a cell array of ripple detections 
            % for the channels in a given area (i.e. 8 channels if the area has 8 micro channels) and a cell array of 
            % start-end indices per ripple.
            % Ripple will be considered “legitimate” if it appears in at least two channels in the area. Ripple in channel X 
            % and ripple in channel Y will be considered the same ripple (and thus appear in both areas) if they are less than 
            % ripplesDistMicrChanForMerge ms apart.
            %
            % The method returns the merged ripple times and ripple start-end times.

            
            nChans = length(ripplesTimes);
            maxRippleTime = -inf;
            for iChan = 1:nChans
                maxRippleTime = max(maxRippleTime,max(ripplesTimes{iChan}));
            end
            dataDuration = maxRippleTime+obj.ripplesDistMicrChanForMerge;
            
            %build ripples matrix
            ripplesMat = zeros(nChans, dataDuration);
            for iChan = 1:nChans
                ripplesMat(iChan,ripplesTimes{iChan}) = 1;
            end
            
            rippleTimesMerged = [];
            rippleStartEndMerged = [];
            for iChan = 1:nChans-1
                currRipTimes = find(ripplesMat(iChan,:));
                for iRipple = 1:length(currRipTimes)
                    minInd = max(1, currRipTimes(iRipple)-obj.ripplesDistMicrChanForMerge);
                    maxInd = min(currRipTimes(iRipple)+obj.ripplesDistMicrChanForMerge,size(ripplesMat,2));
                    otherChansClose = ripplesMat(iChan+1:end,minInd:maxInd);
                    %if the ripple also appears on another channel
                    if any(otherChansClose(:))
                        rippleTimesMerged(end+1) = currRipTimes(iRipple);
                        rippleStartEndMerged(end+1,:) = ripplesStartEnd{iChan}(ripplesTimes{iChan}==currRipTimes(iRipple),:);
                        ripplesMat(iChan+1:end,currRipTimes(iRipple)-obj.ripplesDistMicrChanForMerge:currRipTimes(iRipple)+obj.ripplesDistMicrChanForMerge) = 0;
                    end
                end
            end
        end
        
        function saveRipplesDetectionsMicro(obj, runData, runByChannel, useExistingRipples)
            
            % A wrapper for running ripples detection on micro channels. Ripple detection on micro channels is performed per area, where a ripple will be 
            % considered “legitimate” if it appears in at least two channels in the area. Ripple in channel X and ripple in channel Y will be considered the 
            % same ripple (and thus appear in both areas) if they are less than ripplesDistMicrChanForMerge ms apart (by default – 15 ms).
            % By default, the wrapper detects and saves ripples per area, which means: A. Running and saving ripples for all of the micro channels in that area, 
            % B. Running the ripples merge method (getRipplesFromMicroElectrodesInArea) and saving the merge ripples for the area.
            % The list of areas on which the wrapper will run on is provided as part of the runData input struct, if left empty then the wrapper will run on all 
            % the areas for that patient (based on the micro montage).
            % Another option is to run ripple detection per micro channels, without the merging per area step. This is possible if the input parameter 
            % runByChannel is set to true (by default it’s false) and the list of channels to run on is provided as part of the runData input struct.
            % By default when requested to run on an area the method will first run the detections on all the channels in that area and save them. If 
            % useExistingRipples is true, then before running the detection on a micro channel it will first check whether a ripples file for that channel 
            % already exists and if it does will load it instead of running the detections – this is useful if the detections for some or all of the micro 
            % channels in the area were already run and saved and we only want to run the missing ones + the merging step or only want to run the merging step 
            % if all of them were already saved.
            %
            % The wrapper receives as input runData which is a struct array with the length as the number of patients. Each element (=patient) should have the 
            % fields:
            % PatientName
            % areasToRunOn – a list of areas on which the ripples detection per channel and per area will be run. If left empty the method will run on all the 
            % areas for that patient. If the input parameter runByChannel is true then this field is ignored and the run is by channel and not by area (by 
            % default it’s false).
            % channelsToRunOn – a list of channels for which the ripples detection will be performed and saved. This field is only relevant if the input 
            % parameter runByChannel is true (by default, false), otherwise this field is ignored.
            % microMontageFileName – the file name (including path) of the micromontage.
            % MicroDataFolder – The folder in which the raw micro data is stored. The property fileNamePrefix of the class includes the prefix for the data 
            % filenames (by default: ‘CSC’).
            % MicroRipplesFileNames – name (including path) of the ripple mat files in which the ripple times for the micro channels should be saved. The 
            % method will save ripples per channel as MicroRipplesFileNames<#channel index> and ripples per area as MicroRipplesFileNames<area name>.
            % MicroSpikesFileNames – The filenames (including path) of the mat files that include the spike times. The method assumes the filename format is 
            % SpikesFileNames<#area_name>. That means that for the detection of all the ripples on micro channels in the RAH area (for example) the same spikes 
            % file is loaded with the name <MicroSpikesFileNames>RAH.mat.
            % sleepScoringFileName – file name (including path) of the sleep scoring mat file. If not provided all the data will be used.

            
            if nargin < 3 || isempty(runByChannel)
                runByChannel = false;
            end
            
            if nargin < 4 || isempty(useExistingRipples)
                useExistingRipples = false;
            end
            
            rd = RippleDetector;
             
            %go over the patients
            nPatients = length(runData);
            for iPatient=1:nPatients
                disp(['patient ',runData(iPatient).patientName]);
                %default - run ripples by area
                
                if isfield(runData(iPatient), 'sleepScoringFileName') && ~isempty(runData(iPatient).sleepScoringFileName)
                    try
                        sleepScoring = load(runData(iPatient).sleepScoringFileName);
                        sleepScoring = sleepScoring.sleep_score_vec;
                    catch
                        disp([runData(iPatient).sleepScoringFileName '.mat doesn''t exist']);
                        sleepScoring = [];
                    end
                end
                
                %load micro montage
                microMontage = load(runData(iPatient).microMontageFileName);
                microMontage = microMontage.Montage;
                
                if ~runByChannel
                    allAreas = {microMontage.Area};
                    allChans = [microMontage.Channel];
                    
                    %areas can be provided in runData, or otherwise all
                    %areas are loaded from the micro montage
                    if isfield(runData(iPatient), 'areasToRunOn') && ~isempty(runData(iPatient).areasToRunOn)
                        areasToRunOn = runData(iPatient).areasToRunOn;
                    else %all areas from the micromontage
                        areasToRunOn = unique(allAreas);
                        areasToRunOn = areasToRunOn(cellfun(@(x)~isempty(x),areasToRunOn));
                    end
                    
                    nAreas = length(areasToRunOn);
                    for iArea = 1:nAreas
                        currArea = areasToRunOn{iArea};                        
                        disp(['area ',currArea,' ',num2str(iArea)','/',num2str(nAreas)]);
                        
                        %find micro channels indices
                        currChanInds = allChans(strcmp(allAreas,currArea));
                        nChans = length(currChanInds);
                        
                        %loading spikes for area
                        disp(['loading spikes']);
                        try
                            peakTimes = load([runData(iPatient).MicroSpikesFileNames,currArea,'.mat']);
                            peakTimes = peakTimes.peakTimes;
                        catch
                            disp([runData(iPatient).MicroSpikesFileNames,currArea,'.mat doesn''t exist']);
                            peakTimes = [];
                        end
                        disp('running ripples per channel');
                        ripplesTimesAll = cell(1,nChans);
                        ripplesStartEndAll = cell(1,nChans);
                        for iChan = 1:nChans
                            
                            currChan = currChanInds(iChan);
                            %if useExistingRipples is true, first check
                            %whether a ripples file already exists for this
                            %channel
                            if useExistingRipples
                                try
                                    ripplesData = [runData.MicroRipplesFileNames num2str(currChan) '.mat'];
                                    ripplesData = load(ripplesData);
                                    ripplesTimesAll{iChan} = ripplesData.ripplesTimes;
                                    ripplesStartEndAll{iChan} = ripplesData.ripplesStartEnd;
                                    continue;
                                catch
                                    ripplesTimesAll{iChan} = [];
                                    ripplesStartEndAll{iChan} = [];
                                end
                            end
                 
                            try
                                currData = load([runData(iPatient).MicroDataFolder,'\',obj.dataFilePrefix ,num2str(currChan),'.mat']);
                                currData = currData.data;
                            catch
                                disp([runData(iPatient).MicroDataFolder,'\',obj.dataFilePrefix,num2str(currChan),'.mat doesnt exist']);
                                continue;
                            end
                            
                            [ripplesTimes, ripplesStartEnd] = rd.detectRipple(currData, sleepScoring, peakTimes);
                            save([runData(iPatient).MicroRipplesFileNames,num2str(currChan),'.mat'],'ripplesTimes','ripplesStartEnd');
                            ripplesTimesAll{iChan} = ripplesTimes;
                            ripplesStartEndAll{iChan} = ripplesStartEnd;
                        end
                        [ripplesTimes,ripplesStartEnd] = obj.getRipplesFromAllMicroElectrodesInArea(ripplesTimesAll, ripplesStartEndAll);
                        
                        save([runData(iPatient).MicroRipplesFileNames,currArea,'.mat'],'ripplesTimes','ripplesStartEnd');
                    end
                else %if we want to run for specific channels and not by area, no merging will be performed in this case
                    nChannels = length(runData(iPatient).channelsToRunOn);
                    for iChan = 1:nChans
                        currChan = runData(iPatient).channelsToRunOn(iChan);
                        currArea = microMontage(currChan).Area;
                        disp(['Channel ',num2str(currChan)]);
                        
                        try
                            currData = load([runData(iPatient).MicroDataFolder,'\',obj.dataFilePrefix ,num2str(currChan),'.mat']);
                            currData = currData.data;
                        catch
                            disp([runData(iPatient).MicroDataFolder,'\',obj.dataFilePrefix,num2str(currChan),'.mat doesnt exist']);
                            continue;
                        end
                        
                        %loading spikes for area
                        disp(['loading spikes']);
                        try
                            peakTimes = load([runData(iPatient).MicroSpikesFileNames,currArea,'.mat']);
                            peakTimes = peakTimes.peakTimes;
                        catch
                            disp([runData(iPatient).MicroSpikesFileNames,currArea,'.mat doesn''t exist']);
                            peakTimes = [];
                        end
                        
                        [ripplesTimes, ripplesStartEnd] = rd.detectRipple(currData, sleepScoring, peakTimes);
                        save([runData(iPatient).MicroRipplesFileNames,num2str(currChan),'.mat'],'ripplesTimes','ripplesStartEnd');
                    end
                end
            end
        end
        
        function plotRipplesMicro(obj, runData, areaName, folderToSave)
            % The method plots ripples detected on micro channels. Ripple for single micro channels should be detected and 
            % saved in advance. The method loads ripples for single micro channels in the area, merges them and plots all 
            % the single ripples in figures where each column is a channel and each row is a ripple. Red circles appear in 
            % all the channels for which a ripple was detected – i.e. each row should have at least two red circles as only 
            % ripples that are detected in at least two channels are considered legitimate.
            % It receives as input:
            % A. runData – a struct containing the fields:
            % patientName
            % microMontageFileName – the file name (including path) of the micromontage.
            % MicroDataFolder – folder from which to read the micro data (the method assumes the prefix for the 
            % files is CSC, can be changed by the property dataFilePrefix)
            % MicroRipplesFileNames – name (including path) of the ripple mat files in which the ripple times for the micro 
            % channels are saved (the method assumes the name of the file is MicroRipplesFileNames<#channel index>
            % B. areaName – name of the area (as appears in the micro montage) to plot the ripples for.
            % C. folderToSave (optional) – folder into which to save the figures.

            
            if nargin < 4 || isempty(folderToSave)
                toSave = false;
            else
                toSave = true;
            end
            
            secondBefAfter = obj.secondBefAfter*obj.samplingRate;
            
            %load macro montage
            microMontage = load(runData.microMontageFileName);
            microMontage = microMontage.Montage;
            allAreas = {microMontage.Area};
            allChans = [microMontage.Channel];
            
            currChanInds = allChans(strcmp(allAreas,areaName));
            nChans = length(currChanInds);
            
            ripplesTimes = cell(1,nChans);
            rippleTimesLog = cell(1,nChans);
            ripplesStartEnd = cell(1,nChans);
            datas = cell(1,nChans);
            
            for iChan = 1:nChans
                currChan = currChanInds(iChan);
                
                try
                    ripplesData = [runData.MicroRipplesFileNames num2str(currChan) '.mat'];
                    ripplesData = load(ripplesData);
                    ripplesTimes{iChan} = ripplesData.ripplesTimes;
                    ripplesStartEnd{iChan} = ripplesData.ripplesStartEnd;
                catch
                    disp([runData.MicroRipplesFileNames num2str(currChan) '.mat doesn''t exist']);
                    ripplesTimes = [];
                end
                
                rippleTimesLog{iChan} = zeros(1,max(ripplesTimes{iChan}+obj.ripplesDistMicrChanForMerge));
                rippleTimesLog{iChan}(ripplesTimes{iChan}) = 1;
                
                try
                    datas{iChan} = load([runData.MicroDataFolder,'\',obj.dataFilePrefix ,num2str(currChan),'.mat']);
                    datas{iChan} = datas{iChan}.data;
                catch
                    disp([runData.MicroDataFolder,'\',obj.dataFilePrefix,num2str(currChan),'.mat doesnt exist']);
                    continue;
                end
            end
            
            [rippleTimesMerged,~] = obj.getRipplesFromAllMicroElectrodesInArea(ripplesTimes, ripplesStartEnd);

            nRipples = length(rippleTimesMerged);
            nPlots = ceil(nRipples/obj.nInPlotMicro);
            
            indRipple = 1;
            figInd = 1;
            xax = [-secondBefAfter:secondBefAfter];
            for iPlot = 1:nPlots-1
                f = figure;
                
                for iInPlot = 1:obj.nInPlotMicro
                    for iChan = 1:nChans
                        subplot(obj.nInPlotMicro,nChans,(iInPlot-1)*nChans+iChan);
                        minInd = max(rippleTimesMerged(indRipple)-secondBefAfter,1);
                        maxInd = min(rippleTimesMerged(indRipple)+secondBefAfter,length(datas{iChan}));
                        
                        maxInd2 = min(maxInd,length(rippleTimesLog{iChan}));
                        if minInd<length(rippleTimesLog{iChan})
                            currRippleTime = find(rippleTimesLog{iChan}(minInd:maxInd2));
                            currRippleTime = currRippleTime(currRippleTime>=secondBefAfter-obj.ripplesDistMicrChanForMerge & currRippleTime<=secondBefAfter+obj.ripplesDistMicrChanForMerge);
                        else
                            currRippleTime = [];
                        end
                        
                        plot(xax/obj.samplingRate,datas{iChan}(minInd:maxInd));
                        hold all;
                        if ~isempty(currRippleTime)
                            plot(xax(currRippleTime)/obj.samplingRate,datas{iChan}(minInd+currRippleTime-1),'O','color','r','markersize',10);
                        end
                        
                        if iChan==1
                            title(['Ripple time = ', num2str(rippleTimesMerged(indRipple)/obj.samplingRate/60),' mins']);
                        end
                    end
                    indRipple = indRipple+1;
                    if indRipple >= nRipples
                        break;
                    end
                end
                if toSave
                        set(f, 'Position', get(0, 'Screensize'));
                        saveas(f, [folderToSave,'\',runData.patientName,'_',areaName,'_all_ripples_' num2str(figInd),'.jpg']);
                        close(f);
                else
                    pause;
                end
                figInd = figInd+1;
            end
            
            
        end
        
        function results = runRipSpikesMicro(obj, runData, fileNameResults)
            % This method is the same as runRipSpikes except it uses ripple times as detected on micro channels, and not 
            % macro channels. See documentation for runRipSpikes for more details.
            % Ripple detection on micro channels is performed for an entire area rather than just one channel – see the 
            % documentation on getRipplesFromMicroElectrodesInArea for more details. Because the detection is per area and 
            % not per channel, the difference from runRipSpikes is in that the struct runData instead of a field 
            % channelsToRunOn there should be a field ‘areasToRunOn’. If this field doesn’t exist or is empty the method 
            % will try to run on all areas for that patient according to the micro montages.
            % Two more differences in the input runData struct, it should have the fields:
            % microMontageFileName - the file name (including path) of the micromontage (this is tinstead of 
            % macroMontageFileName).
            % MicroRipplesFileNames - name (including path) of the ripple mat files in which the ripple times for the 
            % micro channels are saved (the method assumes the name of the file is MicroRipplesFileNames <#channel index>). 
            % This is instead of RipplesFileNames.

            
            if nargin < 3
                fileNameResults = '';
            end
            
            shortTimeRangeAfterStim = obj.shortTimeRangeAfterStim*obj.samplingRate;
            midTimeRangeAfterStim = obj.midTimeRangeAfterStim*obj.samplingRate;
            stimulusDuration = obj.stimulusDuration*obj.samplingRate/1000;
            
            nPatients = length(runData);
            for iPatient = 1:nPatients
                disp(['Patient ',runData(iPatient).patientName]);
                results(iPatient).patientName = runData(iPatient).patientName;
                
                %load exp data for stimulation timings
                expData = load(runData(iPatient).ExpDataFileName);
                expData = expData.EXP_DATA;
                stimTimes = expData.stimTiming.validatedTTL_NLX;
                firstStim = stimTimes(1);
                lastStim = stimTimes(end);
                dataDuration = lastStim+midTimeRangeAfterStim;
                
                %load micro montage
                microMontage = load(runData(iPatient).microMontageFileName);
                microMontage = microMontage.Montage;
                allAreas = {microMontage.Area};
                allChans = [microMontage.Channel];
                
                %load spike data
                try
                    spikeData = load(runData(iPatient).spikeData);
                catch
                    disp([runData(iPatient).spikeData ' doesn''t exist, continuing']);
                    continue;
                end
                
                %go over areas
                resultsPerChan = [];
                if isfield(runData(iPatient), 'areasToRunOn') && ~isempty(runData(iPatient).areasToRunOn)
                    areasToRunOn = runData(iPatient).areasToRunOn;
                else %all areas from the micromontage
                    areasToRunOn = unique(allAreas);
                    areasToRunOn = areasToRunOn(cellfun(@(x)~isempty(x),areasToRunOn));
                end
                
                nAreas = length(areasToRunOn);
                for iArea = 1:nAreas
                    currArea = areasToRunOn{iArea};
                    resultsPerChan(iArea).area = currArea;
                    
                    disp(['area ',currArea,' ',num2str(iArea)','/',num2str(nAreas)]);
                    
                    %find micro channels indices 
                    currChanInds = allChans(strcmp(allAreas,currArea));
                    resultsPerChan(iArea).channelNum = currChanInds;
                    nChans = length(currChanInds);
                                        
                    if isfield(runData(iPatient),'MicroRipplesFileNames') && ~isempty(runData(iPatient).MicroRipplesFileNames)
                        
                        ripplesTimes = cell(1,nChans);
                        ripplesStartEnd = cell(1,nChans);
                        
                        for iChan = 1:nChans
                            try
                                ripplesData = [runData(iPatient).MicroRipplesFileNames num2str(currChanInds(iChan)) '.mat'];
                                ripplesData = load(ripplesData);
                                ripplesTimes{iChan} = ripplesData.ripplesTimes;
                                ripInds = ripplesTimes{iChan} <= dataDuration;
                                ripplesTimes{iChan} = ripplesTimes{iChan}(ripInds);
                                ripplesStartEnd{iChan} = ripplesData.ripplesStartEnd(ripInds,:);
                            catch
                                disp([runData(iPatient).MicroRipplesFileNames num2str(currChanInds(iChan)) '.mat doesn''t exist']);
                                ripplesTimes = [];
                                ripplesStartEnd = [];
                            end
                        end
                    else
                        ripplesTimes = {};
                        ripplesStartEnd = {};
                    end
                    
                    [rippleTimesMerged,rippleStartEndMerged] = obj.getRipplesFromAllMicroElectrodesInArea(ripplesTimes, ripplesStartEnd);
                    
                    ripplesTimes = rippleTimesMerged;
                    ripplesStartEnd = rippleStartEndMerged;
                    
                    ripplesBeforeStimInds = ripplesTimes<firstStim-obj.windowSpikeRateAroundRip;
                    %                     tmp = ripplesTimes>=firstStim&ripplesTimes<=stimTimes(end);
                    %only ripples that are short and mid effect and not too
                    %close to the stimulus
                    stimInds = zeros(1,dataDuration);
                    %short effect
                    for iStim = 1:length(stimTimes)
                        stimInds(stimTimes(iStim)+stimulusDuration:stimTimes(iStim)+shortTimeRangeAfterStim) = 1;
                    end
                    stimDiffs = diff(stimTimes);
                    stimIndsWithMidPauseAfter = [find(stimDiffs >= midTimeRangeAfterStim) length(stimTimes)];
                    stimIndsWithMidPauseAfterTimes = stimTimes(stimIndsWithMidPauseAfter);
                    for iStim = 1:length(stimIndsWithMidPauseAfter)
                        stimInds(stimIndsWithMidPauseAfterTimes(iStim)+shortTimeRangeAfterStim:stimIndsWithMidPauseAfterTimes(iStim)+midTimeRangeAfterStim) = 1;
                    end
                    ripplesIndsLog = zeros(1,dataDuration);
                    ripplesIndsLog(ripplesTimes(ripplesTimes<=dataDuration)) = 1;
                    ripplesDuringStimInds = ismember(ripplesTimes,find(ripplesIndsLog & stimInds));
                    
                    %find all units in spike data recorded in area
                    unitInds = find(strcmp(extractfield(spikeData.micro_channels_spike_summary.unit_list_xls,'Location'),currArea));
                    
                    if obj.spikeMultiUnits
                        currChannels = [spikeData.micro_channels_spike_summary.unit_list_xls(unitInds).Channel];
                        uChannels = unique(currChannels);
                        
                        newUnitInds = {};
                        
                        spikeTimes = {};
                        
                        for iChannel = 1:length(uChannels)
                            currUinds = unitInds(currChannels==uChannels(iChannel));
                            %merge spike data to create multi unit
                            currSpikeTimes = sort(cat(1,spikeData.micro_channels_spike_summary.spike_timestamps{currUinds}));
                            currSpikeTimes = currSpikeTimes(currSpikeTimes <= dataDuration);
                            if length(currSpikeTimes)/(dataDuration/obj.samplingRate) >= obj.minSpikeRateToIncludeUnit
                                spikeTimes{end+1} = currSpikeTimes;
                                newUnitInds{end+1} = currUinds;
                            end
                        end
                        
                        unitInds = newUnitInds;
                        nUnits = length(unitInds);
                    else %single units
                        %remove units with too low spike rate
                        newUnitInds = {};
                        spikeTimes = {};
                        for iUnit = 1:length(unitInds)
                            currSpikeTimes = spikeData.micro_channels_spike_summary.spike_timestamps{unitInds(iUnit)};
                            currSpikeTimes = currSpikeTimes(currSpikeTimes <= dataDuration);
                            if length(currSpikeTimes)/(dataDuration/obj.samplingRate) >= obj.minSpikeRateToIncludeUnit
                                newUnitInds{end+1} = unitInds(iUnit);
                                spikeTimes{end+1} = currSpikeTimes;
                            end
                        end
                        unitInds = newUnitInds;
                        nUnits = length(unitInds);
                    end
                    
                    %for each unit - run the spikes-rip correlation before
                    %stim and during stim
                    ripRatesBefore = cell(1,nUnits);
                    controlRatesBefore = cell(1,nUnits);
                    fireRateRipBefore = cell(1,nUnits);
                    fireRateControlBefore = cell(1,nUnits);
                    
                    fireRateRipStim = cell(1,nUnits);
                    fireRateControlStim = cell(1,nUnits);
                    
                    stimTriggeredFireRates = cell(1,nUnits);
                    controlForStimTriggered = cell(1,nUnits);
                    allSessionFireRates = cell(1,nUnits);
                    
                    
                    %                     maxSpikeTime = 0;
                    for iUnit = 1:nUnits
                        [fireRateRipBefore{iUnit}, fireRateControlBefore{iUnit}] = obj.checkSpikeRateAtRip(ripplesTimes(ripplesBeforeStimInds), ripplesStartEnd(ripplesBeforeStimInds,:), spikeTimes{iUnit});
                        [fireRateRipStim{iUnit}, fireRateControlStim{iUnit}] = obj.checkSpikeRateAtRip(ripplesTimes(ripplesDuringStimInds), ripplesStartEnd(ripplesDuringStimInds,:), spikeTimes{iUnit});
                        
                        [allSessionFireRates{iUnit}, stimTriggeredFireRates{iUnit}, controlForStimTriggered{iUnit}] = obj.getSpikeRateAtStim(stimTimes, spikeTimes{iUnit});
                    end
                    
                    %                                         dataDuration = maxSpikeTime+obj.winFromLastSpike;
                    
                    %average of units
                    [fireRateRipAvgUnitsBefore, fireRateControlAvgUnitsBefore] = obj.checkSpikeRateAtRip(ripplesTimes(ripplesBeforeStimInds), ripplesStartEnd(ripplesBeforeStimInds,:), spikeTimes);
                    [fireRateRipAvgUnitsStim, fireRateControlAvgUnitsStim] = obj.checkSpikeRateAtRip(ripplesTimes(ripplesDuringStimInds), ripplesStartEnd(ripplesDuringStimInds,:), spikeTimes);
                    
                    resultsPerChan(iArea).unitInds = unitInds;
                    resultsPerChan(iArea).spikeTimes = spikeTimes;
                    
                    resultsPerChan(iArea).fireRateRipBefore = fireRateRipBefore;
                    resultsPerChan(iArea).fireRateRipAvgUnitsBefore = fireRateRipAvgUnitsBefore;
                    resultsPerChan(iArea).fireRateControlBefore = fireRateControlBefore;
                    resultsPerChan(iArea).fireRateControlAvgUnitsBefore = fireRateControlAvgUnitsBefore;
                    
                    resultsPerChan(iArea).fireRateRipStim = fireRateRipStim;
                    resultsPerChan(iArea).fireRateRipAvgUnitsStim = fireRateRipAvgUnitsStim;
                    resultsPerChan(iArea).fireRateControlStim = fireRateControlStim;
                    resultsPerChan(iArea).fireRateControlAvgUnitsStim = fireRateControlAvgUnitsStim;
                    
                    resultsPerChan(iArea).allSessionFireRates = allSessionFireRates;
                    resultsPerChan(iArea).stimTriggeredFireRates = stimTriggeredFireRates;
                    resultsPerChan(iArea).controlForStimTriggered = controlForStimTriggered;
                    resultsPerChan(iArea).firstStim = firstStim;
                    resultsPerChan(iArea).lastStim = lastStim;
                    
                end
                
                results(iPatient).resultsPerChan = resultsPerChan;
            end
            
            if ~isempty(fileNameResults)
                save(fileNameResults,'results');
            end
            
        end
        
        function plotResultsSpikes(obj, results, folderToSave)
            % The output of runRipSpikes can be plotted using this method.
            % It receives the results struct (the output of runRipSpikes) and folderToSave (optional) for saving the figures.
            
            % The figures per channel show the spike rate for ripple vs control before and during stimulations per 
            % single/multi unit and for the average of units. They also show the stimulation triggered spike rate per unit 
            % and the spike rate for the entire session per unit. This method also calculates the pvalue of the comparison 
            % between the presented conditions (ripple vs control or before vs stimulations), using the average over the 
            % event-triggered window, the window size is defined by windowForSignificance.
            % The maximal number of units shown in the single units figure is set by maxLinesInFigureRipSpike. 
            % The purple in the spike rate figures represent stimulation spochs.

            
            if nargin < 3 || isempty(folderToSave)
                toSave = false;
            else
                toSave = true;
            end
            
            windowForSignificance = round(obj.windowForSignificance*obj.samplingRate/1000);
            indsForSign = obj.windowSpikeRateAroundRip+1-windowForSignificance:obj.windowSpikeRateAroundRip+windowForSignificance;
            avgRippleBeforeAfter = obj.avgRippleBeforeAfter*obj.samplingRate;
            segDuration = (windowForSignificance*2+1)/obj.samplingRate;
            
            nPatients = length(results);
            
            for iPatient = 1:nPatients
                nChans = length(results(iPatient).resultsPerChan);
                
                for iChan = 1:nChans
                    
%                     ymin = inf;
%                     ymax = -inf;
                    nUnits = length(results(iPatient).resultsPerChan(iChan).spikeTimes);
                    if nUnits==0
                        continue;
                    end
                    nFiguresForSingle = ceil(nUnits/obj.maxLinesInFigureRipSpike);
                    if nFiguresForSingle>1
                        nInFigure = obj.maxLinesInFigureRipSpike;
                    else
                        nInFigure = nUnits;
                    end
                    fs = cell(1,nFiguresForSingle+1);
                    for iFigure = 1:nFiguresForSingle
                        fs{iFigure} = figure('Name',['Ripple triggered spike rate single units ',results(iPatient).patientName,' chan ',num2str(results(iPatient).resultsPerChan(iChan).channelNum),' Area ',results(iPatient).resultsPerChan(iChan).area]);
                        if iFigure*obj.maxLinesInFigureRipSpike > nUnits
                            nUnitsInFigure = nUnits - (iFigure-1)*obj.maxLinesInFigureRipSpike;
                        else
                            nUnitsInFigure = obj.maxLinesInFigureRipSpike;
                        end
                        
                        for iUnit = 1:nUnitsInFigure
                            
                            ymin = inf;
                            ymax = -inf;
                            
                            currUnit = (iFigure-1)*obj.maxLinesInFigureRipSpike+iUnit;
                            
                            %                             currRip = results(iPatient).resultsPerChan(iChan).ripRatesBefore{currUnit};
                            %                             currControl = results(iPatient).resultsPerChan(iChan).controlRatesBefore{currUnit};
                            %                             [h1,inds1] = hist(currRip,obj.nBinsHist);
                            %                             h1 = h1/sum(h1);
                            %                             h2 = hist(currControl,inds1);
                            %                             h2 = h2/sum(h2);
                            %                             bar(inds1,[h1' h2']);
                            %                             [~,pval] = ttest(currRip,currControl);
                            %                             title(['Before, unit ',num2str(results(iPatient).resultsPerChan(iChan).unitInds(currUnit)),' pval = ',num2str(pval)]);
                            %                             if iUnit == 1
                            %                                 legend({'Ripples','Control'});
                            %                                 xlabel('Spike rate');
                            %                                 ylabel('Probability');
                            %                             end
                            
                            
                            %                             subplot(nInFigure,4,(iUnit-1)*4+3);
                            %                             currRip = results(iPatient).resultsPerChan(iChan).ripRatesStim{currUnit};
                            %                             currControl = results(iPatient).resultsPerChan(iChan).controlRatesStim{currUnit};
                            %                             [h1,inds1] = hist(currRip,obj.nBinsHist);
                            %                             h1 = h1/sum(h1);
                            %                             h2 = hist(currControl,inds1);
                            %                             h2 = h2/sum(h2);
                            %                             bar(inds1,[h1' h2']);
                            %                             [~,pval] = ttest(currRip,currControl);
                            %                             title(['Stimulations, unit ',num2str(results(iPatient).resultsPerChan(iChan).unitInds(currUnit)),' pval = ',num2str(pval)]);
                            %
                            
                            %first plot - before ripple vs control
                            subplot(nInFigure*2,4,[(iUnit-1)*8+1 (iUnit-1)*8+5]);
                            if ~isempty(results(iPatient).resultsPerChan(iChan).fireRateRipBefore{currUnit})
                                resRipB = results(iPatient).resultsPerChan(iChan).fireRateRipBefore{currUnit};
                                resCont = results(iPatient).resultsPerChan(iChan).fireRateControlBefore{currUnit};
                                nRipplesBefore = size(resRipB,1);
                                
                                shadedErrorBar([-obj.windowSpikeRateAroundRip:obj.windowSpikeRateAroundRip], mean(resRipB), std(resRipB)/sqrt(nRipplesBefore),'lineprops','-g');
                                hold all;
                                shadedErrorBar([-obj.windowSpikeRateAroundRip:obj.windowSpikeRateAroundRip], mean(resCont), std(resCont)/sqrt(nRipplesBefore),'lineprops','-b');
                                hold off;
                                
                                aucRipB = sum(resRipB(:,indsForSign),2);
                                aucCont = sum(resCont(:,indsForSign),2);
                                frB = mean(mean(resRipB(:,indsForSign),2));
                                frC = mean(mean(resCont(:,indsForSign),2));
                                
                                [~,p] = ttest(aucRipB,aucCont,'tail','right');
                                currylim = ylim;
                                ymin = min(currylim(1),ymin);
                                ymax = max(currylim(2),ymax);
                                
                                if iUnit ==1
                                    legend({'Ripple','Control'});
                                end
                                xlabel('Time (ms)');
                                ylabel('Spike rate');
                                if iUnit==1
                                    title({['Control vs Ripple, Before, nRipples = ',num2str(nRipplesBefore)],['Unit ',num2str(results(iPatient).resultsPerChan(iChan).unitInds{currUnit}),' firing rate ripple=',num2str(frB,'%0.2g'),' control=',num2str(frC,'%0.2g'),' p=',num2str(p)]});
                                else
                                    title(['Units ',num2str(results(iPatient).resultsPerChan(iChan).unitInds{currUnit}),' firing rate ripple=',num2str(frB,'%0.2g'),' control=',num2str(frC,'%0.2g'),' p=',num2str(p)]);
                                end
                            else
                                nRipplesBefore = 0;
                                resRipB = [];
                                title('No data');
                            end
                            
                            %second plot - stimulation ripple vs control
                            subplot(nInFigure*2,4,[(iUnit-1)*8+2 (iUnit-1)*8+6]);
                            
                            if ~isempty(results(iPatient).resultsPerChan(iChan).fireRateRipStim{currUnit})
                                resRipS = results(iPatient).resultsPerChan(iChan).fireRateRipStim{currUnit};
                                resCont = results(iPatient).resultsPerChan(iChan).fireRateControlStim{currUnit};
                                nRipplesStim = size(resRipS,1);
                                
                                shadedErrorBar([-obj.windowSpikeRateAroundRip:obj.windowSpikeRateAroundRip], mean(resRipS), std(resRipS)/sqrt(nRipplesStim),'lineprops','-r');
                                hold all;
                                shadedErrorBar([-obj.windowSpikeRateAroundRip:obj.windowSpikeRateAroundRip], mean(resCont), std(resCont)/sqrt(nRipplesStim),'lineprops','-b');
                                hold off;
                                
                                aucRipS = sum(resRipS(:,indsForSign),2);
                                aucCont = sum(resCont(:,indsForSign),2);
                                frS = mean(mean(resRipS(:,indsForSign),2));
                                frC = mean(mean(resCont(:,indsForSign),2));
                                
                                [~,p] = ttest(aucRipS,aucCont,'tail','right');
                                currylim = ylim;
                                ymin = min(currylim(1),ymin);
                                ymax = max(currylim(2),ymax);
                                
                                if iUnit ==1
                                    legend({'Ripple','Control'});
                                end
                                xlabel('Time (ms)');
                                ylabel('Spike rate');
                                if iUnit==1
                                    title({['Control vs Ripple, Stimulations, nRipples = ',num2str(nRipplesStim)],['firing rate ripple=',num2str(frS,'%0.2g'),' control=',num2str(frC,'%0.2g'),' p=',num2str(p)]});
                                else
                                    title(['firing rate ripple=',num2str(frS,'%0.2g'),' control=',num2str(frC,'%0.2g'),' p=',num2str(p)]);
                                end
                            else
                                resRipS = [];
                                nRipplesStim = 0;
                                title('No data');
                            end
                            
                            %third plot - before vs stim
                            subplot(nInFigure*2,4,[(iUnit-1)*8+3 (iUnit-1)*8+7]);
                            if ~isempty(resRipB) && ~isempty(resRipS)
                                shadedErrorBar([-obj.windowSpikeRateAroundRip:obj.windowSpikeRateAroundRip], mean(resRipB), std(resRipB)/sqrt(nRipplesBefore),'lineprops','-g');
                                hold all;
                                shadedErrorBar([-obj.windowSpikeRateAroundRip:obj.windowSpikeRateAroundRip], mean(resRipS), std(resRipS)/sqrt(nRipplesStim),'lineprops','-r');
                                hold off;
                                
                                [~,p] = ttest2(aucRipB,aucRipS);
                                
                                currylim = ylim;
                                ymin = min(currylim(1),ymin);
                                ymax = max(currylim(2),ymax);
                                if iUnit==1
                                    legend({'Before','Stimulations'});
                                end
                                xlabel('Time (ms)');
                                ylabel('Spike rate');
                                if iUnit==1
                                    title({['Before vs Stimulations'], ['firing rate before=',num2str(frB,'%0.2g'),' stim=',num2str(frS,'%0.2g'),' p(two tailed)=',num2str(p)]});
                                else
                                    title(['firing rate before=',num2str(frB,'%0.2g'),' stim=',num2str(frS,'%0.2g'),' p(two tailed)=',num2str(p)]);
                                end
                                %make all ylimits the same
                                for iPlot =  1:3
                                    subplot(nInFigure*2,4,[(iUnit-1)*8+iPlot (iUnit-1)*8+iPlot+4]);
                                    ylim([ymin ymax]);
                                end
                            else
                                title('No data');
                            end
                            
                            subplot(nInFigure*2,4,(iUnit-1)*8+4);
                            
                            if ~isempty(results(iPatient).resultsPerChan(iChan).stimTriggeredFireRates{currUnit})
                                resStim = results(iPatient).resultsPerChan(iChan).stimTriggeredFireRates{currUnit};
                                resStimControl = results(iPatient).resultsPerChan(iChan).controlForStimTriggered{currUnit};
                                nStim = size(resStim,1);
                                shadedErrorBar([-obj.windowSpikeRateAroundStim:obj.windowSpikeRateAroundStim]/obj.samplingRate, mean(resStim), std(resStim)/sqrt(nStim),'lineprops','-m');
                                hold all;
                                shadedErrorBar([-obj.windowSpikeRateAroundStim:obj.windowSpikeRateAroundStim]/obj.samplingRate, mean(resStimControl), std(resStimControl)/sqrt(nStim),'lineprops','-k');
                                if iUnit == 1
                                    xlabel('Time (sec)');
                                    ylabel('Spike rate');
                                    %                             legend({'Stimulation','Control'});
                                    title('Stimulation triggered spike rate');
                                end
                            else
                                title('No data');
                                nStim = 0;
                            end
                            
                            if ~isempty(results(iPatient).resultsPerChan(iChan).allSessionFireRates{currUnit})
                                subplot(nInFigure*2,4,(iUnit-1)*8+8);
                                plot(results(iPatient).resultsPerChan(iChan).allSessionFireRates{currUnit},'color','k');
                                hold all;
                                stimRange = [results(iPatient).resultsPerChan(iChan).firstStim:results(iPatient).resultsPerChan(iChan).lastStim+obj.shortTimeRangeAfterStim*obj.samplingRate];
                                plot(stimRange,results(iPatient).resultsPerChan(iChan).allSessionFireRates{currUnit}(stimRange),'color','m');
                                if iUnit ==1
                                    xlabel('Time (ms)');
                                    ylabel('Spike rate');
                                    title('Spike rate for entire session');
                                end
                            else
                                title('No data');
                            end
                            
                        end
                        suptitle(['Ripples-spikes per unit ',results(iPatient).patientName,' chan ',num2str(results(iPatient).resultsPerChan(iChan).channelNum),' Area ',results(iPatient).resultsPerChan(iChan).area]);
                    end
                    
                    ymin = inf;
                    ymax = -inf;
                    
                    %average over units
                    fs{nFiguresForSingle+1} = figure('Name',['Ripple triggered spike rate - average over units ',results(iPatient).patientName,' chan ',num2str(results(iPatient).resultsPerChan(iChan).channelNum),' Area ',results(iPatient).resultsPerChan(iChan).area]);
                    
                    %first plot - before ripple vs control
                    subplot(1,3,1);
                    if ~isempty(results(iPatient).resultsPerChan(iChan).fireRateRipAvgUnitsBefore)
                        resRipB = results(iPatient).resultsPerChan(iChan).fireRateRipAvgUnitsBefore;
                        resCont = results(iPatient).resultsPerChan(iChan).fireRateControlAvgUnitsBefore;
                        nRipplesBefore = size(resRipB,1);
                        shadedErrorBar([-obj.windowSpikeRateAroundRip:obj.windowSpikeRateAroundRip], mean(resRipB), std(resRipB)/sqrt(nRipplesBefore),'lineprops','-g');
                        hold all;
                        shadedErrorBar([-obj.windowSpikeRateAroundRip:obj.windowSpikeRateAroundRip], mean(resCont), std(resCont)/sqrt(nRipplesBefore),'lineprops','-b');
                        hold off;
                        %ttest around ripple
                        
                        aucRipB = sum(resRipB(:,indsForSign),2);
                        aucCont = sum(resCont(:,indsForSign),2);
                        frB = mean(mean(resRipB(:,indsForSign),2));
                        frC = mean(mean(resCont(:,indsForSign),2));
                        
                        [~,p] = ttest(aucRipB,aucCont,'tail','right');
                        currylim = ylim;
                        ymin = min(currylim(1),ymin);
                        ymax = max(currylim(2),ymax);
                        
                        legend({'Ripple','Control'});
                        xlabel('Time (ms)');
                        ylabel('Spike rate');
                        title({['Control vs Ripple, Before, nRipples = ',num2str(nRipplesBefore)],['firing rate ripple=',num2str(frB,'%0.2g'),' control=',num2str(frC,'%0.2g'),' p=',num2str(p)]});
                    else
                        title('No Data');
                        nRipplesBefore = 0;
                    end
                    
                    %second plot - stimulation ripple vs control
                    subplot(1,3,2);
                    if ~isempty(results(iPatient).resultsPerChan(iChan).fireRateRipAvgUnitsStim)
                        resRipS = results(iPatient).resultsPerChan(iChan).fireRateRipAvgUnitsStim;
                        resCont = results(iPatient).resultsPerChan(iChan).fireRateControlAvgUnitsStim;
                        nRipplesStim = size(resRipS,1);
                        shadedErrorBar([-obj.windowSpikeRateAroundRip:obj.windowSpikeRateAroundRip], mean(resRipS), std(resRipS)/sqrt(nRipplesStim),'lineprops','-r');
                        hold all;
                        shadedErrorBar([-obj.windowSpikeRateAroundRip:obj.windowSpikeRateAroundRip], mean(resCont), std(resCont)/sqrt(nRipplesStim),'lineprops','-b');
                        hold off;
                        
                        aucRipS = sum(resRipS(:,indsForSign),2);
                        aucCont = sum(resCont(:,indsForSign),2);
                        frS = mean(mean(resRipS(:,indsForSign),2));
                        frC = mean(mean(resCont(:,indsForSign),2));
                        
                        [~,p] = ttest(aucRipS,aucCont,'tail','right');
                        currylim = ylim;
                        ymin = min(currylim(1),ymin);
                        ymax = max(currylim(2),ymax);
                        
                        legend({'Ripple','Control'});
                        xlabel('Time (ms)');
                        ylabel('Spike rate');
                        title({['Control vs Ripple, Stimulations, nRipples = ',num2str(nRipplesStim)],['firing rate ripple=',num2str(frS,'%0.2g'),' control=',num2str(frC,'%0.2g'),' p=',num2str(p)]});
                    else
                        title('No Data');
                        nRipplesStim = 0;
                    end
                    
                    %third plot - before vs stim
                    subplot(1,3,3);
                    if ~isempty(results(iPatient).resultsPerChan(iChan).fireRateRipAvgUnitsBefore) && ~isempty(results(iPatient).resultsPerChan(iChan).fireRateRipAvgUnitsStim)
                        shadedErrorBar([-obj.windowSpikeRateAroundRip:obj.windowSpikeRateAroundRip], mean(results(iPatient).resultsPerChan(iChan).fireRateRipAvgUnitsBefore), std(results(iPatient).resultsPerChan(iChan).fireRateRipAvgUnitsBefore)/sqrt(nRipplesBefore),'lineprops','-g');
                        hold all;
                        shadedErrorBar([-obj.windowSpikeRateAroundRip:obj.windowSpikeRateAroundRip], mean(results(iPatient).resultsPerChan(iChan).fireRateRipAvgUnitsStim), std(results(iPatient).resultsPerChan(iChan).fireRateRipAvgUnitsStim)/sqrt(nRipplesStim),'lineprops','-r');
                        hold off;
                        legend({'Before','Stimulations'});
                        xlabel('Time (ms)');
                        ylabel('Spike rate');
                        
                        [~,p] = ttest2(aucRipB,aucRipS);
                        currylim = ylim;
                        ymin = min(currylim(1),ymin);
                        ymax = max(currylim(2),ymax);
                        title({['Before vs Stimulations'], ['firing rate before=',num2str(frB,'%0.2g'),' stim=',num2str(frS,'%0.2g'),' p(two tailed)=',num2str(p)]});
                        
                        for iPlot =  1:3
                            subplot(1,3,iPlot);
                            ylim([ymin ymax]);
                        end
                    else
                        title('No Data');
                    end
                    
                    %                     subplot(2,2,1);
                    %                     currRip = results(iPatient).resultsPerChan(iChan).ripRatesAvgUnitsBefore;
                    %                     currControl = results(iPatient).resultsPerChan(iChan).controlRatesAvgUnitsBefore;
                    %                     [h1,inds1] = hist(currRip,obj.nBinsHist);
                    %                     h1 = h1/sum(h1);
                    %                     h2 = hist(currControl,inds1);
                    %                     h2 = h2/sum(h2);
                    %                     bar(inds1,[h1' h2']);
                    %                     [~,pval] = ttest(currRip,currControl);
                    %                     title(['Before, pval = ',num2str(pval)]);
                    %                     legend({'Ripples','Control'});
                    %                     xlabel('Spike rate');
                    %                     ylabel('Probability');
                    
                    %                     subplot(2,2,2);
                    %                     shadedErrorBar([-obj.windowSpikeRateAroundRip:obj.windowSpikeRateAroundRip], results(iPatient).resultsPerChan(iChan).fireRateAvgRipAvgUnitsBefore, results(iPatient).resultsPerChan(iChan).fireRateErrRipAvgUnitsBefore);
                    %                     title('Ripples-centered spikes rate');
                    %                     xlabel('Miliseconds');
                    %                     ylabel('Spike rate');
                    %
                    %                     subplot(2,2,3);
                    %                     currRip = results(iPatient).resultsPerChan(iChan).ripRatesAvgUnitsStim;
                    %                     currControl = results(iPatient).resultsPerChan(iChan).controlRatesAvgUnitsStim;
                    %                     [h1,inds1] = hist(currRip,obj.nBinsHist);
                    %                     h1 = h1/sum(h1);
                    %                     h2 = hist(currControl,inds1);
                    %                     h2 = h2/sum(h2);
                    %                     bar(inds1,[h1' h2']);
                    %                     [~,pval] = ttest(currRip,currControl);
                    %                     title(['Stimulations, pval = ',num2str(pval)]);
                    %                     legend({'Ripples','Control'});
                    %                     xlabel('Spike rate');
                    %                     ylabel('Probability');
                    %
                    %                     subplot(2,2,4);
                    %                     shadedErrorBar([-obj.windowSpikeRateAroundRip:obj.windowSpikeRateAroundRip], results(iPatient).resultsPerChan(iChan).fireRateAvgRipAvgUnitsStim, results(iPatient).resultsPerChan(iChan).fireRateErrRipAvgUnitsStim);
                    %                     title('Ripples-centered spikes rate');
                    %                     xlabel('Miliseconds');
                    %                     ylabel('Spike rate');
                    
                    suptitle(['Average over units, ',results(iPatient).patientName,' chan ',num2str(results(iPatient).resultsPerChan(iChan).channelNum),' Area ',results(iPatient).resultsPerChan(iChan).area]);
                    
                    %set y limits to be the same in all figures
                    %                     for iFigure = 1:nFiguresForSingle+1
                    %                         nChild = numel(fs{iFigure}.Children);
                    %                         for iChild = 1:nChild
                    %                             if strcmp(fs{iFigure}.Children(iChild).Type,'axes')
                    %                                 fs{iFigure}.Children(iChild).YLim = [ymin ymax];
                    %                             end
                    %                         end
                    %                     end
                    %
                    if toSave
                        
                        if obj.spikeMultiUnits
                            filename = '\ripple_spikes_multi_units_';
                        else
                            filename = '\ripple_spikes_single_units_';
                        end
                        
                        for iFigure = 1:nFiguresForSingle
                            set(fs{iFigure}, 'Position', get(0, 'Screensize'));
                            saveas(fs{iFigure},[folderToSave,'\',results(iPatient).patientName,filename, num2str(results(iPatient).resultsPerChan(iChan).channelNum) '_',num2str(iFigure),'.jpg']);
                            close(fs{iFigure});
                        end
                        
                        set(fs{end}, 'Position', get(0, 'Screensize'));
                        saveas(fs{end},[folderToSave,'\',results(iPatient).patientName,'\ripple_spikes_avg_units_' num2str(results(iPatient).resultsPerChan(iChan).channelNum),'.jpg']);
                        close(fs{end});
                    end
                    
                    %                     ylim(fs, [ymin ymax]);
                end
            end
        end
        
        function plotResultsRipplesData(obj, results, folderToSave)
            
            %The method receives as input a struct with the format of the output of runRipplesData and produces figures. 
            %The input folderToSave (optional) sets the folder into which the figures will be saved.
            
            if nargin < 3 || isempty(folderToSave)
                toSave = false;
            else
                toSave = true;
            end
            
            avgRippleBeforeAfter = obj.avgRippleBeforeAfter*obj.samplingRate;
            
            nPatients = length(results);
            
            for iPatient = 1:nPatients
                nChans = length(results(iPatient).resultsPerChan);
                
                for iChan = 1:nChans
                    %ripple summary figure - average ripple and spectogram
                    %before and during stimulation
                    f1 = figure('Name',['Ripples Data patient ',results(iPatient).patientName,' channel ',num2str(results(iPatient).resultsPerChan(iChan).channelNum),' Area ',results(iPatient).resultsPerChan(iChan).area]);
                    
                    subplot(4, 3, [1 4]);
                    meanSpecs = results(iPatient).resultsPerChan(iChan).meanSpecs;
                    meanEpochs = results(iPatient).resultsPerChan(iChan).meanEpochs;
                    imagesc(obj.xaxisForRipSp/obj.samplingRate, obj.freqRangeSpRip, meanSpecs(:,obj.specStartPointRipSp+1:end-obj.specStartPointRipSp));
                    set(gca, 'YDir','normal');
                    xlabel('Time (sec)');
                    ylabel('Frequency (Hz)');
                    colorbar;
                    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
                    
                    yyaxis right;
                    if length(meanEpochs)>1
                        plot(obj.xaxisForRipSp/obj.samplingRate,meanEpochs(obj.specStartPointRipSp+1:end-obj.specStartPointRipSp),'color','k','linewidth',1);
                    end
                    title({['Power spectrum around spindles, Before'],['nSpindles=',num2str(results(iPatient).resultsPerChan(iChan).nEpochs)]});
                    ylabel('uV');
                    
                    subplot(4,3,[7 10]);
                    currAngles = results(iPatient).resultsPerChan(iChan).SIanglesSpRip;
                    if ~isempty(currAngles)
                        circm = circ_mean(currAngles(~isnan(currAngles))');
                        tmpplot = polarhistogram(currAngles(~isnan(currAngles)),obj.nBinsPolar,'Normalization','probability');
                        hold all;
                        polarplot([circm circm],[0 max(tmpplot.Values)/2],'color','r','linewidth',2);
                        hold off;
                        nSpindles = sum(~isnan(currAngles));
                        title({['pval r = ',num2str(results(iPatient).resultsPerChan(iChan).r)], ['v = ', num2str(results(iPatient).resultsPerChan(iChan).v)], ['nSpindles = ',num2str(nSpindles)]});
                    else
                        title('No spindles');
                    end
                    
                    subplot(4,3,2);
                    shadedErrorBar([-avgRippleBeforeAfter:avgRippleBeforeAfter]/obj.samplingRate,results(iPatient).resultsPerChan(iChan).avgBefore,results(iPatient).resultsPerChan(iChan).stdBefore/sqrt(results(iPatient).resultsPerChan(iChan).nRipplesBefore));
                    xlabel('Time (sec)');
                    ylabel('uV');
                    title(['Average of detected ripple - before, nRipples=',num2str(results(iPatient).resultsPerChan(iChan).nRipplesBefore)]);
                    
                    subplot(4,3,3);
                    shadedErrorBar([-avgRippleBeforeAfter:avgRippleBeforeAfter]/obj.samplingRate,results(iPatient).resultsPerChan(iChan).avgStim,results(iPatient).resultsPerChan(iChan).stdStim/sqrt(results(iPatient).resultsPerChan(iChan).nRipplesStim));
                    xlabel('Time (sec)');
                    ylabel('uV');
                    title(['Average of detected ripple - stimulation, nRipples=',num2str(results(iPatient).resultsPerChan(iChan).nRipplesStim)]);
                    
                    subplot(4,3,5);
                    plot(obj.freqoiForAvgSpec,results(iPatient).resultsPerChan(iChan).specBefore);
                    xlabel('Frequency');
                    title('Spectrum of average 0-10 Hz');
                    
                    subplot(4,3,6);
                    plot(obj.freqoiForAvgSpec,results(iPatient).resultsPerChan(iChan).specStim);
                    xlabel('Frequency');
                    title('Spectrum of average 0-10 Hz');
                    
                    
                    cmin = min(min(results(iPatient).resultsPerChan(iChan).meanTFRRipBefore(:)),min(results(iPatient).resultsPerChan(iChan).meanTFRRipStim(:)));
                    cmax = max(max(results(iPatient).resultsPerChan(iChan).meanTFRRipBefore(:)),max(results(iPatient).resultsPerChan(iChan).meanTFRRipStim(:)));
                    subplot(4,3,8);
                    if size(results(iPatient).resultsPerChan(iChan).meanTFRRipBefore,1)>1
                        imagesc([-obj.timeBeforeAfterEventRipSpec*obj.samplingRate:obj.timeBeforeAfterEventRipSpec*obj.samplingRate]/obj.samplingRate,obj.freqRangeForAvgSpec, results(iPatient).resultsPerChan(iChan).meanTFRRipBefore,[cmin cmax]);
                        set(gca, 'YDir','normal');
                        xlabel('Time (sec)');
                        ylabel('frequency');
                        title('Average of ripple triggerd TFR - before');
                        colorbar;
                    else
                        title('No ripples (before)');
                    end
                    
                    subplot(4,3,9);
                    if size(results(iPatient).resultsPerChan(iChan).meanTFRRipStim,1)>1
                        imagesc([-obj.timeBeforeAfterEventRipSpec*obj.samplingRate:obj.timeBeforeAfterEventRipSpec*obj.samplingRate]/obj.samplingRate,obj.freqRangeForAvgSpec, results(iPatient).resultsPerChan(iChan).meanTFRRipStim,[cmin cmax]);
                        set(gca, 'YDir','normal');
                        xlabel('Time (sec)');
                        ylabel('frequency');
                        title('Average of ripple triggerd TFR - stimulation');
                        colorbar;
                    else
                        title('No ripples (stim)');
                    end
                    
                    freqIndsSpindlesRange = find(obj.freqRangeForAvgSpec==obj.freqRangeForShowingSpindles(1)):find(obj.freqRangeForAvgSpec==obj.freqRangeForShowingSpindles(end));
                    if size(results(iPatient).resultsPerChan(iChan).meanTFRRipBefore,1)>1
                        TFRSmallBefore = results(iPatient).resultsPerChan(iChan).meanTFRRipBefore(freqIndsSpindlesRange,:);
                    else
                        TFRSmallBefore = nan;
                    end
                    if size(results(iPatient).resultsPerChan(iChan).meanTFRRipStim,1)>1
                        TFRSmallStim = results(iPatient).resultsPerChan(iChan).meanTFRRipStim(freqIndsSpindlesRange,:);
                    else
                        TFRSmallStim = nan;
                    end
                    
                    cmin = min(min(TFRSmallBefore(:)),min(TFRSmallStim(:)));
                    cmax = max(max(TFRSmallBefore(:)),max(TFRSmallStim(:)));
                    subplot(4,3,11);
                    if size(TFRSmallBefore,1)>1
                        imagesc([-obj.timeBeforeAfterEventRipSpec*obj.samplingRate:obj.timeBeforeAfterEventRipSpec*obj.samplingRate]/obj.samplingRate,obj.freqRangeForShowingSpindles, TFRSmallBefore,[cmin cmax]);
                        set(gca, 'YDir','normal');
                        xlabel('Time (sec)');
                        ylabel('frequency');
                        title('Average of ripple triggerd TFR zoomed - before');
                        colorbar;
                    else
                        title('No ripples');
                    end
                    
                    subplot(4,3,12);
                    if size(TFRSmallStim,1)>1
                        imagesc([-obj.timeBeforeAfterEventRipSpec*obj.samplingRate:obj.timeBeforeAfterEventRipSpec*obj.samplingRate]/obj.samplingRate,obj.freqRangeForShowingSpindles, TFRSmallStim,[cmin cmax]);
                        set(gca, 'YDir','normal');
                        xlabel('Time (sec)');
                        ylabel('frequency');
                        title('Average of ripple triggerd TFR zoomed - stimulation');
                        colorbar;
                    else
                        title('No ripples');
                    end
                    
                    suptitle(['Ripples Channel Data ',results(iPatient).patientName,' ',num2str(num2str(results(iPatient).resultsPerChan(iChan).channelNum)),' Area ',results(iPatient).resultsPerChan(iChan).area]);
                    
                    if toSave
                        set(f1, 'Position', get(0, 'Screensize'));
                        saveas(f1,[folderToSave,'\',results(iPatient).patientName,'\ripple_channel_data_' num2str(results(iPatient).resultsPerChan(iChan).channelNum),'.jpg']);
                        close(f1);
                    end
                end
            end
        end
        
        function plotRipplesBiPolar(obj, runData, folderToSave)
            %
            
            nPatients = length(runData);
            
            for iPatient = 1:nPatients
                ss = load(runData(iPatient).sleepScoringFileName);
                ss = ss.sleep_score_vec;
                
                nChannelCouples = size(runData(iPatient).biPolarCouples,1);
                
                for iCouple = 1:nChannelCouples
                    
                    currChan = runData(iPatient).biPolarCouples(iCouple,1);
                    refChan = runData(iPatient).biPolarCouples(iCouple,2);
                    
                    disp('loading data');
                    try
                        currData = load([runData(iPatient).DataFolder,'\',obj.dataFilePrefix ,num2str(currChan),'.mat']);
                        currData = currData.data;
                    catch
                        disp([runData(iPatient).DataFolder,'\',obj.dataFilePrefix,num2str(currChan),'.mat doesnt exist']);
                        continue;
                    end
                    try
                        refData = load([runData(iPatient).DataFolder,'\',obj.dataFilePrefix ,num2str(refChan),'.mat']);
                        refData = refData.data;
                    catch
                        disp([runData(iPatient).DataFolder,'\',obj.dataFilePrefix,num2str(refChan),'.mat doesnt exist']);
                        continue;
                    end
                    
                    disp(['loading spikes']);
                    try
                        peakTimes = load([runData(iPatient).SpikesFileNames,num2str(currChan),'.mat']);
                        IISChan = peakTimes.peakTimes;
                    catch
                        disp([runData(iPatient).SpikesFileNames,num2str(currChan),'.mat doesn''t exist']);
                        IISChan = [];
                    end
                    
                    try
                        peakTimes = load([runData(iPatient).SpikesFileNames,num2str(refChan),'.mat']);
                        IISRef = peakTimes.peakTimes;
                    catch
                        disp([runData(iPatient).SpikesFileNames,num2str(refChan),'.mat doesn''t exist']);
                        IISRef = [];
                    end
                    
                    disp('detecting ripples');
                    rippleTimes = obj.detectRipple(currData-refData,ss,[IISChan IISRef]);
                    
                    disp('plotting and saving ripples');
                    obj.plotRipples(currData,rippleTimes,[folderToSave,'\',runData(iPatient).patientName,'_',num2str(currChan),'_',num2str(refChan)]);
                end
            end
        end
    end
end