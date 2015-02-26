function [val,Fs,gain,description,ecg_idx,bp_idx] = mvqrs_loadrecord(recordName)
%
% [val,Fs,gain,description,ecg_idx,bp_idx] = mvqrs_loadrecord(recordName)
%
% Load physiological data using wfdb toolbox
%
% Required Parameters:
%
% recordName
%       String specifying the name of the record in the WFDB path or
%       in the current directory.
%
%
% Written by Marcus Vollmer, 2015
% Last Modified: February 16, 2015
% Version 0.2
% 
%endOfHelp

try
    load([recordName 'm.mat']); % val contains signals   
    gain = [];
    val = val'; 
catch
    % [~,sig,~] = rdsamp(recordName);
    % gain = ones(1,size(sig,2));
    error(['cannot read ' recordName 'm.mat'])
end 

siginfo = wfdbdesc(recordName); %computationally expensive!  
description = {siginfo.Description};


% get Sampling Frequency
Fs = siginfo.SamplingFrequency;    
    if ~isnumeric(Fs) || isnan(Fs)
        if strcmp(Fs(end-1:end),'Hz')
            tmp = strfind(Fs,' ');
            if isempty(tmp)
                Fs = str2num(Fs(1:end-2));
            else
                Fs = str2num(Fs(1:tmp-1));
            end
        else
            [~,~,Fs] = rdsamp(recordName,1,1,1); 
        end
    end 
Fs = double(Fs);

% Gain
if isempty(gain)
    gain_info = {siginfo.Gain};
    gain = zeros(1,size(gain_info,2));
    for i=1:size(gain_info,2)
        tmp = gain_info{i};
        postmp = strfind(tmp,' ');
        gain(i) = str2num(tmp(1:postmp(1)-1));
    end
    val = val./repmat(gain,size(val,1),1);
end


% Signal description
    load('signallist.mat');
    ecg_idx = zeros(size(signallistecg,1),size(description,2));
    bp_idx = zeros(size(signallistbp,1),size(description,2));

    % remove white space
    description(~cellfun(@isempty,description)) = strtrim(description(~cellfun(@isempty,description)));
   
    % result match description ECG
    for list_ecg=1:size(signallistecg,1)
       ecg_idx(list_ecg,:) = strcmpi(description,strtrim(signallistecg{list_ecg}));
    end
    ecg_idx = find(sum(ecg_idx)>0);

    % result match description BP
    for list_bp=1:size(signallistbp,1)
       bp_idx(list_bp,:) = strcmpi(description,strtrim(signallistbp{list_bp}));
    end
    bp_idx = find(sum(bp_idx)>0);



end
