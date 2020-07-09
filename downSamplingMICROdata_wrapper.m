mc = microChannels;
ii = 1;


% 485 - to do - need to re-extract MICROs
%this struct can be longer than 1, where each element represents a patient
runDataForMicroSaving(ii).microChannelsFolderToLoad = 'E:\Data_p\p485\EXP8\spikeSorting';
runDataForMicroSaving(ii).microChannelsFolderToSave = 'E:\Data_p\p485\EXP8\Denoised_Downsampled_InMicroVolt\MICRO';

% 486 - DONE (2018)
runDataForMicroSaving(1).microChannelsFolderToSave = 'E:\Data_p\p486\EXP8\Denoised_Downsampled_InMicroVolt\MICRO';

% 487 - DONE (2018)
runDataForMicroSaving(1).microChannelsFolderToSave = 'E:\Data_p\p487\EXP3\Denoised_Downsampled_InMicroVolt\MICRO';

% 488 - reDONE
ii = ii+1;
%this struct can be longer than 1, where each element represents a patient
runDataForMicroSaving(ii).microChannelsFolderToLoad = 'E:\Data_p\p488\EXP4\spikeSorting';
runDataForMicroSaving(ii).microChannelsFolderToSave = 'E:\Data_p\p488\EXP4\Denoised_Downsampled_InMicroVolt\MICRO';
runDataForMicroSaving(ii).microChannelsToRunOn = 1:9;
mc.oldSamplingRate = 40000;
mc.Nfiles = 7;
mc.subFilesExist = 1;

% runDataForMicroSaving(1).microChannelsToRunOn = []; % Optional - if not
% provided will run on all the channels (files in the form: CSC#_#) in the
% folder microChannelsFolderToLoad.

% 498 - redone
ii = ii+1;
%this struct can be longer than 1, where each element represents a patient
runDataForMicroSaving(ii).microChannelsFolderToLoad = 'E:\Data_p\p498\EXP3\spikeSorting\STIM_ONLY';
runDataForMicroSaving(ii).microChannelsFolderToSave = 'E:\Data_p\p498\EXP3\Denoised_Downsampled_InMicroVolt\MICRO';
runDataForMicroSaving(ii).microChannelsToRunOn = 1;
mc.oldSamplingRate = 32000;
mc.subFilesExist = 0;
mc.fileNumberLength = 0;

% 499
ii = ii+1;
%this struct can be longer than 1, where each element represents a patient
runDataForMicroSaving(ii).microChannelsFolderToLoad = 'E:\Data_p\p499\EXP8\spikeSorting\STIM_ONLY';
runDataForMicroSaving(ii).microChannelsFolderToSave = 'E:\Data_p\p499\EXP8\Denoised_Downsampled_InMicroVolt\MICRO';
mc.oldSamplingRate = 32000;
mc.subFilesExist = 0;
mc.fileNumberLength = 0;
mc.fileSuffixLength = 4;

% 520 
ii = ii+1;
%this struct can be longer than 1, where each element represents a patient
runDataForMicroSaving(ii).microChannelsFolderToLoad = 'E:\Data_p\p520\EXP4\spikeSorting';
runDataForMicroSaving(ii).microChannelsFolderToSave = 'E:\Data_p\p520\EXP4\Denoised_Downsampled_InMicroVolt\MICRO';

%Note the code assumes the format of the filename of the parts is
%CSC<chan_number>_<part number in 3 digits exactly>. if it's different the
%properties microChannels.fileNamePrefix (currently 'CSC') and
%microChannels.fileNumberLength (currently 3) should be changed accordingly

%These are the defaults, change them if it's not true
mc.samplingRate = 1000; % required sampling rate
mc.mergeAndDownsample(runDataForMicroSaving(ii));
