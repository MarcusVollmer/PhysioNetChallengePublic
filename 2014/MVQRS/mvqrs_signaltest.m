function [ValidSignals,ValidSignalsSetting,Fs,gain,val] = mvqrs_signaltest(varargin)
%
% [ValidSignals,ValidSignalsSetting,Fs,gain,val] = mvqrs_signaltest(varargin)
%
% Load physiological data and determines signals to be used for heart beat
% annotation.
%
% Required Parameters:
%
% recordName
%       String specifying the name of the record in the WFDB path or
%       in the current directory.
%
% Optional Parameters are:
%
% threshold
%       The threshold factor is a value between 0 and 1. The higher the
%       factor the less annotations. (default: 0.5)
% sec
%       An integer variable which specifies the length of the time series
%       which will be used in order to determine valid signals (in
%       seconds). (default: 30)
% downsampling
%       An integer variable which specifies the frequency (in Hz) of the
%       time series which will be triggered by omission. (default: 80)
%
%
%
% Written by Marcus Vollmer, 2015
% Last Modified: January 26, 2015
% Version 0.2
%
% %Example:
% wfdb2mat('mitdb/200')
% [ValidSignals,ValidSignalsSetting,Fs,gain,val] = mvqrs_signaltest('200',30,100);
% 
%endOfHelp


%Set default pararameter values
inputs = {'recordName','threshold','sec','downsampling'};
threshold = 0.5;
sec = 30;
downsampling = 80;
for n=1:nargin
    if(~isempty(varargin{n}))
        eval([inputs{n} '=varargin{n};'])
    end
end
% - add here: checkup for variable formats and existence of record file


% load record
[val,Fs,gain,description,ecg_idx,bp_idx] = mvqrs_loadrecord(recordName);
sig = val;

% read specified seconds of all signals
if ~isempty(sec)
    SL = size(sig,1);  
    if SL/Fs>2*sec
        sig = sig(round(SL/2)-sec*Fs:round(SL/2),:);
    end
end


% Downsampling by omission
factor = floor(Fs/downsampling);
Fs = Fs/factor;
sig = sig(1:factor:end,:); % omit information

   
% initialization part
AnnotationECG = zeros(size(sig));
AnnotationBP = zeros(size(sig));
results = zeros(size(sig,2),4);


% Generate annotation files for each signal and for BP and ECG settings
pct=.25; R=.4;
Beat_min=50; Beat_max=220;

for i=1:size(sig,2) 
    for j=1:2
        switch j
            case 1
                wl_tma = ceil(.2*Fs); wl_we = ceil(.2*Fs);
                [Ann,~] = mvqrs_ann(sig(:,i),Fs,wl_tma,pct,wl_we,Beat_min,Beat_max,threshold,R); % valid annotations
                AnnotationECG(Ann,i) = 1;                
            case 2
                wl_tma = ceil(Fs); wl_we = ceil(.4*Fs);
                [Ann,~] = mvqrs_ann(sig(:,i),Fs,wl_tma,pct,wl_we,Beat_min,Beat_max,threshold,R); % valid annotations
                AnnotationBP(Ann,i) = 1;                
        end    
        
        results(i,j) = size(Ann,1); %number of annotations
        if size(Ann,1)>=floor((size(sig,1)/Fs)*24/60)
            % 24 Hz assumed to be the minimum heart rate
            RRinterval = diff(Ann);
            relRRinterval = RRinterval(2:end)./RRinterval(1:end-1);
            results(i,2+j) = sum(relRRinterval<=1.2 & relRRinterval>=.8)/size(relRRinterval,1);
        end
    end
end


% choose signals and fix setting/ usage of description
ValidSignals = results(:,3:4)>=.8;

% result match description ECG
if ~isempty(ecg_idx)
    ecg_idx = setdiff((ValidSignals(ecg_idx,1)==1)'.*ecg_idx,0);
end

% result match description BP
if ~isempty(bp_idx)
    bp_idx = setdiff((ValidSignals(bp_idx,1)==1)'.*bp_idx,0);
end

% Other signals
other_ind = setdiff(find(sum(ValidSignals,2)>0),[ecg_idx bp_idx]);
if isempty(ecg_idx)    
    tmp = intersect(find(ValidSignals(:,1)),other_ind);
    if ~isempty(tmp)
        ecg_idx = tmp(1);
        other_ind = setdiff(other_ind,tmp(1));
    else
        % take signals as ECG reference with "ECG" description
        ecg_idx = zeros(size(signallistecg,1),size(description,2));
        for list_ecg=1:size(signallistecg,1)
           ecg_idx(list_ecg,:) = strcmp(description,signallistecg{list_ecg});
        end
        ecg_idx = find(sum(ecg_idx)>0);

        if isempty(ecg_idx)
            % take first signal as ECG reference
            ecg_idx = 1;
        end
        other_ind = setdiff(other_ind,ecg_idx);
    end
end
other_bp_ind = intersect(find(ValidSignals(:,2)),other_ind);
other_ecg_ind = setdiff(other_ind,other_bp_ind);


% vector of valid signals
ValidSignals = [ecg_idx(:); other_ecg_ind(:);  bp_idx(:);  other_bp_ind(:); ]';
%description(ValidSignals)

% vector of settings
ValidSignalsSetting = [zeros(1,length(ecg_idx)) ones(1,length(other_ecg_ind)) 2*ones(1,length(bp_idx)) 3*ones(1,length(other_bp_ind))];

Fs = Fs*factor;
end
