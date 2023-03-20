classdef AnalyzeSleepOsc < handle
    properties
        samplingRate = 1000;
        
        windowAroundIIS = 500; %ms

        %filtering constants
        defaultFilterOrder = 1;
        nanWarning = 0.01;

        %detections
        dataFilePrefix = 'CSC';
        bipolarDataFilePrefix = 'CSC_biPolar_';
        
       
        %sleep scoring parameters
        scoringEpochDuration = 0.001; % How many seconds represented by one individual value in the scoring vector [scalar].
        sleepEpochs = 1; % all the values in the scoring vector which represent sleep stages for which we want to perform the analysis (like NREM/REM/transitions) [1D vector].
                
    end
    
    methods
        
        function saveDetectionResults(obj, runData, whatToRun)
            % Saves detections of the various types of events, a wrapper for the various detectors.
            % Receives runData with information regarding for which patients to run. In addition receives whatToRun, a struct
            % that sets which detections will be run, with the fields:
            % runSpikes – if true will run IIS detection
            % runRipples – run RippleDetector.detectRipples (Staresina’s method).
            % runRipplesBiPolar – references each MTL channel to the provided reference channel (see documentation of runData)
            % and then runs on it RippleDetector.detectRipples.
            % runSpindles – runs SpindleDetector.detectSpindles (i.e. the method of Andrillon, Nir et al).
            % runSpindlesStaresina - runs SpindleDetector.detectSpindlesStaresina (i.e. the method of Staresina et al).
            % runSWStaresina – runs SlowWavesDetector.findSlowWavesStaresina, i.e. Staresina’s method.
            % runSWMaingret - runs SlowWavesDetector.findSlowWavesMaingret, i.e. Maingret’s method.
            %
            % If whatToRun is not provided, the default is to run ripples, SWStaresina and SpindlesStaresina.
            %
            % The input runData is a struct in the length of number of patients (for which the detection is required). Each
            % element (=patient) in runData should include the fields:
            % patientName
            % channelsToRunOn – an array of channels for which the detections should be run (can be empty only if the only
            % detection required is bipolar ripples, in which case the field biPolarCouples should be provided).
            % biPolarCouples – if bipolar ripples detection is required, should be provided as a matrix where each row has a
            % couple of channel indices – the first is the channel for which ripples should be detected, the second is its
            % reference channel (i.e. ripples will be detected on channel#1-channel#2).
            % DataFolder – The folder in which the raw data files are saved (the method assumes the prefix for the files is CSC,
            % can be changed by the property dataFilePrefix).
            % sleepScoringFileName – file name (including path) of the sleep scoring mat file. If not provided all the data will
            % be used.
            % SpikesFileNames - name (including path) of the spikes mat files in which the spikes times for the macro channels
            % are saved (the method assumes the name of the file is SpikesFileNames <#channel index>). This file name will be
            % used to save the detection results if the user sets whatToRun.runSpikes = true, and it will be used to load spikes
            % before other detections if whatToRun.runSpikes = false. If not provided spikes will not be removed from the data
            % for the analysis.
            % The following fields are required depending on the detection that the user wishes to run:
            % RipplesFileNames – if ripple detection (not bipolar) is required, name (including path) of the ripple mat files in
            % which the ripple times for the macro channels will be saved (will save for each channel as
            % RipplesFileNames<#channel index>).
            % RipplesBipolarFileNames – if bipolar ripple detection is required, name (including path) of the bipolar ripple mat
            % files in which the ripple times for the macro channels will be saved (will save for each channel as
            % RipplesBipolarFileNames <#channel index>).
            % SpindlesFileNames - if spindles (not Staresina’s method) detection is required, name (including path) of the
            % spindles mat files in which the spindles times for the macro channels will be saved (will save for each channel
            % as SpindlesFileNames <#channel index>).
            % SpindlesStaresinaFileNames - if spindles (Staresina’s method) detection is required, name (including path) of
            % the spindles mat files in which the spindles times for the macro channels will be saved (will save for each
            % channel as SpindlesStaresinaFileNames <#channel index>).
            % SWStaresinaFileName - if slow waves (Staresina’s method) detection is required, name (including path) of the
            % slow waves mat files in which the slow waves times for the macro channels will be saved (will save for each
            % channel as SWStaresinaFileName <#channel index>).
            % SWMaingretFileName - if slow waves (Maingret’s method) detection is required, name (including path) of the slow
            % waves mat files in which the slow waves times for the macro channels will be saved (will save for each channel
            % as SWMaingretFileName <#channel index>).
            
            %go over the patients
            nPatients = length(runData);
            for iPatient=1:nPatients
                disp(['patient ',runData(iPatient).patientName]);
                nChannels = length(runData(iPatient).channelsToRunOn);
                
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
                
                %go over requried channels
                for iChannel = 1:nChannels
                    
                    disp(['channel ',num2str(runData(iPatient).channelsToRunOn(iChannel))]);
                    currChan = runData(iPatient).channelsToRunOn(iChannel);
                    
                    %load the data
                    try
                        currData = load([runData(iPatient).DataFolder,'\',obj.dataFilePrefix ,num2str(currChan),'.mat']);
                        currData = currData.data;
                    catch
                        disp([runData(iPatient).DataFolder,'\',obj.dataFilePrefix,num2str(currChan),'.mat doesnt exist']);
                        continue;
                    end
                    
                    %run detections as required
                    if whatToRun.runSpikes
                        disp('running spikes');
                        IIS_det = SpikeWaveDetector;
                        if ~exist('sleepScoring','var')
                            disp('sleep scoring missing')
                            sleepScoring = ones(1,length(currData));
                        end
                        [peakTimes, peakStats] = IIS_det.detectTimes(currData(:)', true, sleepScoring);
                        save([runData(iPatient).SpikesFileNames,num2str(currChan),'.mat'],'peakTimes','peakStats','IIS_det');
                    else
                        disp(['loading spikes']);
                        try
                            peakTimes = load([runData(iPatient).SpikesFileNames,num2str(currChan),'.mat']);
                            peakTimes = peakTimes.peakTimes;
                        catch
                            disp([runData(iPatient).SpikesFileNames,num2str(currChan),'.mat doesn''t exist']);
                            peakTimes = [];
                        end
                    end
                    
                    if whatToRun.runRipples
                        disp('running ripples');
                        rd = RippleDetector_class;
                        if ~isfield(runData(iPatient),'stimulation_times')
                            runData(iPatient).stimulation_times = [];
                        end
                        [ripplesTimes, ripplesStartEnd] = rd.detectRipple(currData, sleepScoring, peakTimes, runData(iPatient).stimulation_times);
                        save([runData(iPatient).RipplesFileNames,num2str(currChan),'.mat'],'ripplesTimes','ripplesStartEnd');
                    end
                    
                    if whatToRun.runSpindles
                        clear LL
                        if ~isempty(dir([runData(iPatient).SpindlesFileNames,num2str(currChan),'.mat']))
                            LL = load([runData(iPatient).SpindlesFileNames,num2str(currChan),'.mat'],'isVerified');
                        end
                        if exist('LL','var') & isfield(LL,'isVerified')
                            disp(['ch',num2str(currChan),'spindle file exists'])
                        else
                            disp('running spindles');
                            sd = SpindleDetector;
                            
                            isVerified = sd.verifyChannelStep1(currData,sleepScoring,peakTimes);
                            if isVerified
                                returnStats = 1;
                                [spindlesTimes,spindleStats,spindlesStartEndTimes] = sd.detectSpindles(currData, sleepScoring, peakTimes, returnStats);
                            else
                                spindlesTimes = []; spindlesStartEndTimes = [];spindleStats = [];
                            end
                            save([runData(iPatient).SpindlesFileNames,num2str(currChan),'.mat'],'isVerified','spindlesTimes','spindlesStartEndTimes','spindleStats');
                        end
                    end
                    
                    
                     if whatToRun.HighFreqSpindles
                        
                        if ~isempty(dir([runData(iPatient).HighFreqSpindlesFileNames,num2str(currChan),'.mat']))
                            LL = load([runData(iPatient).HighFreqSpindlesFileNames,num2str(currChan),'.mat'],'isVerified');
                        end
                        if exist('LL','var') & isfield(LL,'isVerified')
                            disp(['ch',num2str(currChan),'spindle file exists'])
                        else
                            disp('running spindles');
                            sd = SpindleDetector;
                            
                            isVerified = sd.verifyChannelStep1(currData,sleepScoring,peakTimes);
                            if isVerified
                                returnStats = 0;
                                sd.spindleRangeMin = 11;
                                [spindlesTimes,spindleStats,spindlesStartEndTimes] = sd.detectSpindles(currData, sleepScoring, peakTimes);
                            else
                                spindlesTimes = []; spindlesStartEndTimes = [];spindleStats = [];
                            end
                            save([runData(iPatient).HighFreqSpindlesFileNames,num2str(currChan),'.mat'],'isVerified','spindlesTimes','spindlesStartEndTimes','spindleStats');
                        end
                     end
                    
                     
                    if whatToRun.runSpindlesStaresina
                        disp('running spindles');
                        sd = SpindleDetector;
                        [spindlesTimes,spindlesStartEndTimes] = sd.detectSpindlesStaresina(currData, sleepScoring, peakTimes);
                        save([runData(iPatient).SpindlesStaresinaFileNames,num2str(currChan),'.mat'],'spindlesTimes','spindlesStartEndTimes');
                    end
                    
                    if whatToRun.runSWStaresina
                        disp('running slow waves staresina');
                        swd = SlowWavesDetector;
                        slowWavesTimes = swd.findSlowWavesStaresina(currData, sleepScoring, peakTimes);
                        save([runData(iPatient).SWStaresinaFileName,num2str(currChan),'.mat'],'slowWavesTimes');
                    end
                    
                    if whatToRun.runSWMaingret
                        disp('running slow waves maingret');
                        swd = SlowWavesDetector;
                        slowWavesTimes = swd.findSlowWavesMaingret(currData, sleepScoring, peakTimes);
                        save([runData(iPatient).SWMaingretFileName,num2str(currChan),'.mat'],'slowWavesTimes');
                    end
                    
                end
                
                %ripples bipolar running is performed according to the
                %biPolarCouples field (and not channelsToRunOn field)
                if whatToRun.runRipplesBiPolar
                    disp('running ripples bi polar');
                    nCouples = size(runData(iPatient).biPolarCouples,1);
                    for iCouple = 1:nCouples
                        disp(['channel ',num2str(runData(iPatient).biPolarCouples(iCouple,1)),' ref ',num2str(runData(iPatient).biPolarCouples(iCouple,2))]);
                        currChan = runData(iPatient).biPolarCouples(iCouple,1);
                        refChan = runData(iPatient).biPolarCouples(iCouple,2);
                        try
                            currData = load([runData(iPatient).DataFolder,'\',obj.dataFilePrefix ,num2str(currChan),'.mat']);
                            currData = currData.data;
                        catch
                            disp([runData(iPatient).DataFolder,'\',obj.dataFilePrefix,num2str(currChan),'.mat doesnt exist']);
                            continue;
                        end
                        %load the reference channel
                        try
                            currRef = load([runData(iPatient).DataFolder,'\',obj.dataFilePrefix ,num2str(refChan),'.mat']);
                            currRef = currRef.data;
                        catch
                            disp([runData(iPatient).DataFolder,'\',obj.dataFilePrefix,num2str(refChan),'.mat doesnt exist']);
                            continue;
                        end
                        
                        %load spikes from both channels (data and
                        %reference)
                        disp(['loading spikes']);
                        try
                            peakTimesChan = load([runData(iPatient).SpikesFileNames,num2str(currChan),'.mat']);
                            peakTimesChan = peakTimesChan.peakTimes;
                        catch
                            disp([runData(iPatient).SpikesFileNames,num2str(currChan),'.mat doesn''t exist']);
                            peakTimesChan = [];
                        end
                        try
                            peakTimesRef = load([runData(iPatient).SpikesFileNames,num2str(refChan),'.mat']);
                            peakTimesRef = peakTimesRef.peakTimes;
                        catch
                            disp([runData(iPatient).SpikesFileNames,num2str(refChan),'.mat doesn''t exist']);
                            peakTimesRef = [];
                        end
                        
                        rd = RippleDetector;
                        rd.samplingRate = obj.samplingRate;
                        %run detection on data-reference
                        [ripplesTimes, ripplesStartEnd] = rd.detectRipple(currData-currRef, sleepScoring, [peakTimesChan peakTimesRef]);
                        save([runData(iPatient).RipplesBipolarFileNames,num2str(currChan),'.mat'],'ripplesTimes','ripplesStartEnd');
                    end
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
        
    end
end
