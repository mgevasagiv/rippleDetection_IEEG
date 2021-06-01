%the struct runData holds data about patients and where the different event
%types are stored
DETECTION_RUN = 0;

data_p_path = 'E:\Data_p\';

patients = {'p488','p487','p489','p486','p499','p485','p490','p496','p498','p510','p505','p515','p510','p520','p497'};
expNames = {'EXP4','EXP3','EXP3','EXP8','EXP8','EXP8','EXP3','EXP8','EXP3','EXP7','EXP1','EXP9','EXP1','EXP3','EXP2'};
sleepScoreFileName = {'sleepScore_manualValidated_p488_4_LPHG2',...
    'sleepScore_manualValidated_p487_3_ROF2',...
    'sleepScore_manualValidated_p489_3_RPHG3',...
    'sleepScore_manualValidated_p486_8_RA1',...
    'sleepScore_manualValidated_p499_8_RPT2',...
    'sleepScore_manualValidated_p485_8_ROF7_0016',...
    'sleepScore_manualValidated_p490_3_RSTG1',...
    'sleepScore_manualValidated_p496_8_RSO4',...
    'sleepScore_manualValidated_p498_3_ROF1',...
    'sleepScore_manualValidated_p510_7_RPSM2',...
    'sleepScore_manualValidated_p505_1_RSTG1',...
    'sleepScore_manualValidated_p515_9_LOF3',...
    'sleepScore_manualValidated_p510_1_RPSM2',...
    'sleepScore_manualValidated_p520_3_CSC26',...
    'sleepScore_manualValidated_p497_2_RPHG4'};


%note macro montage is not necessary to run AnalyzeCoupling
% macroMontageNames = {'MM488','MM487','MM489','MM486','MM499','MM485','MM490','MM496','MM498','MM510','MM505','MM515','MM510_1'};

%channels on which detections will be performed (just an example)
channelsPerPatient = {[78  1  23  ],... % p488
    [16  10  39 31 60],... % p487
    [1  8  36  29  44 ],... % p489
    [8 46],... % p486
    [46 38 24],... % p499
    [52 ],... % p485
    [1],... % 490 - RANH - excluded from analysis
    [1],...% 496
    [32 40 48],... % 498
    [9 47],... % p510
    [2],... % p505
    [9 39 65],... %p515
    [9 47],...% p510
    [17 9 ],...% 520 
    [9 16 1 40 48]};  % 497


%for bipolar ripple detection - in every row the first index is the channel in which ripple
%detection is required and the second is the reference channel
% patients = {'p488','p487','p489','p486','p499','p485','p490','p496','p498','p510','p505','p515','p510','p520','p497'};
biPolarCouplesPerPatient = {[78 81;1 4;23 26],... % 488
                            [16 18; 10 11; 39 41; 31 32; 60 62],... % 487
                            [1 4; 8 11; 36 39; 29 33; 44 46],...% 489
                            [8 11; 46 48],... % p486
                            [46 49; 38 40; 24 25],... % p499
                            [52 54],... % p485
                            [1 3; 8 10; 50 53],... % 490 - RANH - excluded from analysis
                            [1 3],...% 496
                            [32 35; 39 42; 40 42; 48 50],...%498
                            [9 13; 47 51],... %510
                            [2 4],... 505
                            [9 12; 39 42; 40 42; 65 67; 66 67],... % 515
                            [9 13; 47 51],... %510
                            [17 18; 9 10],...% 520
                            [9 10; 16 19; 1 2; 40 41; 48 49]};  % 497

%building the run data, not that for all file names of detections the
%methods assume the name is <provided name according to runData><channel
%num>, e.g. if runData(iPatient).SWStaresinaFileName='c:\slow_wave', then
%the slow waves file for channel 1 is 'c:\slow_wave1.mat'.
runData = [];
nPatients = length(patients);
for iPatient = 1:nPatients 
    runData(iPatient).patientName = patients{iPatient};
    %The folder where the raw data is stored - you will need to change it
    runData(iPatient).DataFolder = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MACRO'];
       
    %The folder+filename from which spindles are going to be loaded
    %The folder+filename from which spindles are going to be loaded (should be the same
    %as the line above if the detections are first run and saved into SpindlesStaresinaFileNames
    spindleFolder = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MACRO\spindleResults'];
    runData(iPatient).SpindlesFileNames = fullfile(spindleFolder,'spindleTimes');
    if isempty(dir(spindleFolder))
        mkdir(spindleFolder)
    end
    
    %The folder+filename into which the bipolar ripples detections results is going
    %to be stored
    rippleFolder = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MACRO\rippleBipolarResults'];
    runData(iPatient).RipplesBipolarFileNames = fullfile(rippleFolder,'rippleTimes');
    if isempty(dir(rippleFolder))
        mkdir(rippleFolder)
    end
    
    %The folder+filename from which ripples are going to be loaded (should be the same
    %as the line above if the bipolar detections are first run and saved into
    %RipplesBipolarFileNames)
    rippleFolder = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MACRO\rippleBipolarResults'];
    runData(iPatient).RipplesFileNames = fullfile(rippleFolder,'rippleTimes');
    if isempty(dir(rippleFolder))
        mkdir(rippleFolder)
    end
    
    %list of couples for bipolar ripple detection - where each row has the channel
    %in which the detection is performed in the first index, and the
    %reference channel in the second
    runData(iPatient).biPolarCouples = biPolarCouplesPerPatient{iPatient};
    %The folder+filename into which the spikes results is going to be stored or is
    %already stored if the spikes detection was already run (the folder should
    %exist)
    runData(iPatient).SpikesFileNames = fullfile(runData(iPatient).DataFolder, ...
        sprintf('MacroInterictalSpikeTimesFor_%s_%s_',patients{iPatient},expNames{iPatient}));
   
    %name of the EXP data for the patient
    runData(iPatient).ExpDataFileName = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\',patients{iPatient},'_',expNames{iPatient},'_dataset.mat'];
    % Extract stimulation times
    mm = matfile(runData(iPatient).ExpDataFileName); EXP_DATA = mm.EXP_DATA;
    runData(iPatient).stimulation_times = EXP_DATA.stimTiming.validatedTTL_NLX;
    clear mm EXP_DATA
    
    %name of the sleep scoring mat file for the patient
    runData(iPatient).sleepScoringFileName = [runData(iPatient).DataFolder,'\',sleepScoreFileName{iPatient},'.mat'];
    %channels that the ripples analyses will be performed on
    
    % For ripple detection
    if DETECTION_RUN
        runData(iPatient).channelsToRunOn = biPolarCouplesPerPatient{iPatient}(:);   
    else
    % For analysis
        runData(iPatient).channelsToRunOn = biPolarCouplesPerPatient{iPatient}(:,1);   
    end
    %macromontage file name per patient
    runData(iPatient).macroMontageFileName = [data_p_path,'MACRO_MONTAGE','\',patients{iPatient},'\',expNames{iPatient},'\MacroMontage.mat'];
    %The micro montage filename is required for micro channels analysis
    %which is not performed in this example
    %runData(iPatient).microMontageFileName = ['C:\Users\user\google drive\MayaProject\montages\',patients{iPatient},'\',expNames{iPatient},'\montage.mat'];
    
    %name of file with single units info
    runData(iPatient).spikeData = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\averagedRef\',patients{iPatient},'_spike_timestamps_post_processing.mat'];
end


%% an example for saving ripples using the wrapper AnalyzeCoupling.saveDetectionResults
ac = AnalyzeCoupling;

%setting which detections to run - in this example, bipolar ripples will be
%detected for channelsToRunOn per patient
whatToRun.runSpikes = true;
whatToRun.runRipples = false;
whatToRun.runRipplesBiPolar = true;
whatToRun.runSpindles = false;
whatToRun.runSpindlesStaresina = false;
whatToRun.runSWStaresina = false;
whatToRun.runSWMaingret = false;

% Prepare a biPolar version of MTL channels
ac.saveBipolarDataset(runData);

%saving detections (in this example, bipolar ripples detection)
parfor ii = 1:length(runData)
    ac.saveDetectionResults(runData(ii), whatToRun);
end

%% an example for detecting ripples directly using RippleDetector (it's the same thing the wrapper does inside)
rd = RippleDetector;

%an example of using the ripple detection directly and not with the wrapper
%(on the first channel of the first patient for this example)
currChan = runData(iPatient).channelsToRunOn(1);

%loading - sleep scoring, IIS, data
sleepScoring = load(runData(1).sleepScoringFileName);
sleepScoring = sleepScoring.sleep_score_vec;
peakTimes = load([runData(iPatient).SpikesFileNames,num2str(currChan),'.mat']);
peakTimes = peakTimes.peakTimes;
currData = load([runData(iPatient).DataFolder,'\CSC',num2str(currChan),'.mat']);
currData = currData.data;
%detecting the ripples
[ripplesTimes, ripplesStartEnd] = rd.detectRipple(currData, sleepScoring, peakTimes);
%plotting the single ripples and saving the figures
rd.plotRipples(currData,ripplesTimes,'D:\singleRipples');

%% ripple related analyses

% Plot ripple summary figures per patient
for ii = 1:length(runData)
    resultsFile = fullfile('E:\Data_p\ClosedLoopDataset\rippleDetResults\MACRO\',...
        sprintf('resultsRipplesMACRO_%s_%s.mat',runData(ii).patientName,expNames{ii} ));
    runData(ii).channelsToRunOn = biPolarCouplesPerPatient{ii}(:,1);
    if isempty(dir(resultsFile))
        resultsData = rd.runRippleData(runData(ii),resultsFile);
    else
        disp(['plotting results ',runData(ii).patientName])
        mm = matfile(resultsFile);
        resultsData = mm.results;
        figureFolder = 'E:\Data_p\ClosedLoopDataset\rippleDetResults\MACROrippleFigures_biPolar';
        
        % Compare to unipolar detections without a REF :
        % figureFolder = 'E:\Data_p\ClosedLoopDataset\rippleDetResults\MACROrippleFigures_uniPolar';
        rd.plotResultsRipplesData(resultsData,figureFolder);
    end
end

% collecting results in one struct
for ii = 1:length(runData)
    resultsFile = fullfile('E:\Data_p\ClosedLoopDataset\rippleDetResults\MACRO\',...
        sprintf('resultsRipplesMACRO_%s_%s.mat',runData(ii).patientName,expNames{ii} ));
    if isempty(dir(resultsFile))
        continue
    else
        disp(['uploading results ',runData(ii).patientName])
        mm = matfile(resultsFile);
        resultsData(ii) = mm.results;
    end
end


% Alternatively - run for all patients 
%runs the "ripple overview" analysis for the channels specified in
%channelsToRunOn
resultsData = rd.runRippleData(runData,...
    'E:\Data_p\ClosedLoopDataset\rippleDetResults\MACRO\resultsRipplesMACRO.mat');
%plots the results
rd.plotResultsRipplesData(resultsData,'E:\Data_p\ClosedLoopDataset\rippleDetResults\\ripplesFigures');

%runs the ripples-spikes correlation analysis for the channels specified in
%channelsToRunOn
spikeResultsFilename = 'E:\Data_p\ClosedLoopDataset\rippleDetResults\MACRO\resultsRipplesRS.mat';
resultsRS = rd.runRipSpikes(runData,spikeResultsFilename);
%plots the results
mm = matfile(spikeResultsFilename);
resultsRS = mm.results;
rd.plotResultsSpikes(resultsRS,'E:\Data_p\ClosedLoopDataset\rippleDetResults\MACROrippleFigures_biPolar\ripplesSpikesFigures');

%% Population figure for all patients
% using bi-polar detections 



nChan = 0;
nRipThreshold = 20;
nRip = [];
ptNum_v = []; chNum_v = [];
vtest_v = []; rtest_v = [];
chMTL_area = cell(1,1);
avgBefore = zeros(1,2001);
stdBefore = zeros(1,2001);
nSU = 0;
clear avgBefore stdBefore
for iiP = 1:length(runData)
    resultsFile = fullfile('E:\Data_p\ClosedLoopDataset\rippleDetResults\MACRO\',...
        sprintf('resultsRipplesMACRO_%s_%s.mat',runData(iiP).patientName,expNames{iiP} ));
    resultsData = load(resultsFile, 'results');
    resultsData = resultsData.results;
    ptNum = str2num(resultsData.patientName(2:end));
    
    for iiC = 1:length(resultsData.resultsPerChan)
        if resultsData.resultsPerChan(iiC).nRipplesBefore > nRipThreshold
            if mean(resultsData.resultsPerChan(iiC).stdBefore) > 500
                continue
            end
                
            A = resultsData.resultsPerChan(iiC).avgBefore(995:1005);
            [pks,locs,w,p]  = findpeaks(A,'SortStr','descend','NPeaks',1);
            
            if isempty(locs)
                continue
            end
            nChan = nChan + 1;

            D1 = resultsData.resultsPerChan(iiC).avgBefore;
            D2 = resultsData.resultsPerChan(iiC).stdBefore;
            L = length(D1);
            if locs == 7
                avgBefore(nChan,:) = D1;
                stdBefore(nChan,:) = D2;
            elseif locs == 8
                avgBefore(nChan,:) = [D1(2:L),0];
                stdBefore(nChan,:) = [D2(2:L),0];
            elseif locs == 9
                avgBefore(nChan,:) = [D1(3:L),0,0];
                stdBefore(nChan,:) = [D2(3:L),0,0];
            elseif locs == 6
                avgBefore(nChan,:) = [0,D1(1:L-1)];
                stdBefore(nChan,:) = [0,D2(1:L-1)];
            else
                disp(['maxima ', num2str(locs)])
                disp(['ch #', num2str(nChan)])
            end
            
            
            ptNum_v = [ptNum_v, ptNum];
            chNum_v = [chNum_v, resultsData.resultsPerChan(iiC).channelNum];
            area = classifyArea(resultsData.resultsPerChan(iiC).area);
            if ~area.isMTL; disp(ptNum); disp(resultsData.resultsPerChan(iiC).area); warning('wrong area included in ripple detection');end
            chMTL_area{nChan} = area;
            vtest_v = [vtest_v, resultsData.resultsPerChan(iiC).v];
            rtest_v = [rtest_v, resultsData.resultsPerChan(iiC).r];
            nRip(nChan) = resultsData.resultsPerChan(iiC).nRipplesBefore;

            meanTFRRipBefore(nChan,:,:) = resultsData.resultsPerChan(iiC).meanTFRRipBefore;
            
            % Adding spike data where it exists
            if ~isempty(resultsRS(iiP).resultsPerChan)
                for iiU = 1:length(resultsRS(iiP).resultsPerChan(iiC).unitInds)
                    nSU = nSU + 1;
                    SU_mat(nSU,:) = [nChan, nSU];
                    meanfireRateRipBefore(nSU,:) = nanmean(resultsRS(iiP).resultsPerChan(iiC).fireRateRipBefore{iiU},1);
                    meanfireRateControlBefore(nSU,:) = nanmean(resultsRS(iiP).resultsPerChan(iiC).fireRateControlBefore{iiU},1);
                end
            end
            
        end
    end 
end

rippleInfo.chNum_v = chNum_v;
rippleInfo.ptNum_v = ptNum_v;
rippleInfo.area = chMTL_area;
rippleInfo.rtest_v = rtest_v;
rippleInfo.meanTFRRipBefore = meanTFRRipBefore;
rippleInfo.avgBefore = avgBefore;
rippleInfo.stdBefore = stdBefore;
rippleInfo.nRip = nRip;
rippleInfo.meanfireRateRipBefore = meanfireRateRipBefore;
rippleInfo.meanfireRateControlBefore= meanfireRateControlBefore;
rippleInfo.SU_mat = SU_mat;
rippleInfo.nSU = nSU;
rippleInfo.artifactInd = [6 12 34 35];

find(rippleInfo.ptNum_v == 487 & rippleInfo.chNum_v == 39);

save(fullfile('E:\Data_p\ClosedLoopDataset\rippleDetResults\MACRO\','rippleInfo'),'rippleInfo')

%% Figure 4 (and sup figure 4) ripple panels
rd = RippleDetector;

mm = matfile(fullfile('E:\Data_p\ClosedLoopDataset\rippleDetResults\MACRO\','rippleInfo'));
rippleInfo = mm.rippleInfo; 
chMTL_area = rippleInfo.area;
nRip = rippleInfo.nRip;
avgBefore = rippleInfo.avgBefore;
stdBefore = rippleInfo.stdBefore;
meanTFRRipBefore = rippleInfo.meanTFRRipBefore;
SU_mat = rippleInfo.SU_mat;
meanfireRateRipBefore = rippleInfo.meanfireRateRipBefore;
meanfireRateControlBefore = rippleInfo.meanfireRateControlBefore;
nSU = rippleInfo.nSU;
artifactInd = rippleInfo.artifactInd;
nRipThreshold = 20;

nChan = size(avgBefore,1);


for ii_a = 1:4
    if ii_a == 1
        hipInd = zeros(1,length(chMTL_area));
        for ii = 1:length(chMTL_area)
            area = chMTL_area{ii};
            if area.isHip
                hipInd(ii) = 1;
            end
        end
        ind = find(hipInd);
        figName = 'MACROrippleGrandAverage_hip';
    elseif ii_a == 2
        hipInd = zeros(1,length(chMTL_area));
        for ii = 1:length(chMTL_area)
            area = chMTL_area{ii};
            if area.isEC
                hipInd(ii) = 1;
            end
        end
        ind = find(hipInd);
        figName = 'MACROrippleGrandAverage_EC';
    elseif ii_a == 3
        hipInd = zeros(1,length(chMTL_area));
        for ii = 1:length(chMTL_area)
            area = chMTL_area{ii};
            if area.isPHG
                hipInd(ii) = 1;
            end
        end
        ind = find(hipInd);
        figName = 'MACROrippleGrandAverage_PHG';
    elseif ii_a == 4
        for ii = 1:length(chMTL_area)
            area = chMTL_area{ii};
            if area.isMTL
                hipInd(ii) = 1;
            end
        end
        ind = find(hipInd);
        figName = 'MACROrippleGrandAverage';
    end
    
    ind(ismember(ind,artifactInd)) = [];
    ind(ismember(ind,nRip < nRipThreshold)) = [];
    
    NRip = sum(nRip(ind));

    % stats
    disp(ii_a)
    disp(sprintf('%d electrodes, %d pts, %d ripples',length(ind),length(unique(rippleInfo.ptNum_v(ind))),NRip))
    

    avgWeighted = zeros(1,length(avgBefore)); stdWeighted = zeros(1,length(avgBefore));
    meanTFRRipWeighted = zeros(250,length(avgBefore));
    for iiC = ind
        avgWeighted = avgWeighted + (nRip(iiC)/NRip) * avgBefore(iiC,:);
        stdWeighted = stdWeighted + (sqrt( (nRip(iiC)-1)* stdBefore(iiC,:).^2 )/(NRip-iiC));
        meanTFRRipWeighted(:,:) = meanTFRRipWeighted + (nRip(iiC)/NRip) * squeeze(meanTFRRipBefore(iiC,:,:)) ;
    end
    
    NRip = sum(nRip(ind));
    avgWeighted = zeros(1,length(avgBefore)); stdWeighted = zeros(1,length(avgBefore));
    meanTFRRipWeighted = zeros(250,length(avgBefore));
    for iiC = ind
        avgWeighted = avgWeighted + (nRip(iiC)/NRip) * avgBefore(iiC,:);
        stdWeighted = stdWeighted + (sqrt( (nRip(iiC)-1)* stdBefore(iiC,:).^2 )/(NRip-iiC));
        meanTFRRipWeighted(:,:) = meanTFRRipWeighted + (nRip(iiC)/NRip) * squeeze(meanTFRRipBefore(iiC,:,:)) ;
    end
    
    
    f0 = figure('Name', figName,'NumberTitle','off');
    set(gcf,'DefaultAxesFontSize',8);
    set(gcf,'DefaultAxesFontName','arial');
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0.2 0.2 15 7]); % this size is the maximal to fit on an A4 paper when printing to PDF
    set(gcf,'PaperOrientation','portrait');
    set(gcf,'Units','centimeters','Position', get(gcf,'paperPosition')+[1 1 0 0]);
    colormap('jet');
    
    % Ripple waveform average
    axes('position',[0.1,0.58,0.4,0.3])
    shadedErrorBar(1:length(avgBefore),avgWeighted,stdWeighted/sqrt(nChan))
    axis tight
    set(gca,'xlim',[700 1300]);
    XLIM = get(gca,'xlim');
    
    if ii_a == 4
        line([XLIM(1),XLIM(1)]+10,[15 45],'color','k','linewidth',7)
    elseif ismember(ii_a,[2,3])
        line([XLIM(1),XLIM(1)]+10,[10 20],'color','k','linewidth',7)
    else
        line([XLIM(1),XLIM(1)]+10,[15 45],'color','k','linewidth',7)
    end
    
    line([XLIM(1), XLIM(1)+100] ,[-20 -20],'color','k','linewidth',7)
    axis off
    
    % Spike rate triggered by ripples
    axes('position',[0.62,0.1,0.3,0.3])
    suInd = find(ismember(SU_mat(:,1),ind));
    meanM = nanmean(meanfireRateRipBefore(suInd,:));
    stdM = nanstd(meanfireRateRipBefore(suInd,:));
    shadedErrorBar(1:length(meanM),meanM,stdM/sqrt(nSU))
    hold on
    meanM = nanmean(meanfireRateControlBefore(suInd,:));
    stdM = nanstd(meanfireRateControlBefore(suInd,:));
    shadedErrorBar(1:length(meanM),meanM,stdM/sqrt(length(suInd)),'lineProps',{'color',[0.2,.7,0.2]})
    axis tight
    set(gca,'xtick',[0 500 1000], 'xticklabels',[-0.5,0,0.5]);
    if ii_a == 1
        set(gca,'ylim',[0 30], 'yticklabels',[0:10:30]);
    else
        set(gca,'ylim',[0 20], 'yticklabels',[0,20]);
    end
    XLIM = get(gca,'xlim');
  
    
    % Average TFR 
    axes('position',[0.12,0.1,0.4,0.4])
    meanTFRRipWeighted(isnan(meanTFRRipWeighted)) = 0;
    M = ceil(max(max(meanTFRRipWeighted)));
    h = imagesc(meanTFRRipWeighted, [-M,M]);
    set(gca,'xtick',[1,1000,2000],'xticklabels',[-1,0,1])    
    set(gca,'ytick',[5,50:50:250])    
    axis xy
    
    cc = colorbar;
    cc.Ticks = [-M,0,M];
    cc.TickLabels = [-M,0,M]*100;
    
    [spectrum,ntaper,freqoi]  = ft_specest_mtmfft(avgBefore(:,500:1500),[1:length(avgBefore)]/rd.samplingRate,'freqoi',[0.5:0.1:20],'taper','dpss','tapsmofrq',8);
    specBefore = 10*log10(abs(squeeze(spectrum(1,1,:)).^2));
    % specBefore = specBefore/max(specBefore);
    
    axes('position',[0.62,0.58,0.15,0.4])
    plot(freqoi, specBefore', 'k','linewidth',2)
    set(gca,'ylim',[-15 25],'ytick',[-15 0 25])    
    box off
        
    outputFigureFolder = 'E:\Dropbox\ILLUSTRATOR\Figure4\Links\';
    res =  600;
    a = gcf;
    eval(['print ', [outputFigureFolder,'\',a.Name], ' -f', num2str(a.Number),sprintf(' -dtiff  -r%d',res), '-cmyk' ]); % adding r600 slows down this process significantly!
end

