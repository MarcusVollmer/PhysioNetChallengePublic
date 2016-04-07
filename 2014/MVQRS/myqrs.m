function myAnn = myqrs(varargin)
%
% myAnn = myqrs(record,...)
%
% Octave-compatible code for QRS detection
%
% Required Parameters:
%
% record
%       Mat file specifying the waveform.
%
% The record can be followed by parameter/value pairs to specify
% additional properties of the QRS detection.
%
% Fs
%       The sampling frequency in Hz. (default: 1000)
% threshold
%       The threshold factor is a value between 0 and 1. The higher the
%       factor the less annotations. (default: 0.5)
% sec
%       An integer variable which specifies the length of the time series
%       which will be used in order to determine valid channels (in
%       seconds). (default: [], no signal selection, all signals will be
%       used)
% downsampling
%       An integer variable which specifies the frequency (in Hz) of the
%       time series which will be triggered by omission. (default: 80)
%
% This function has the output argument myAnn, which is a vector containing
% the annotation for each detected beat in the record.
%
% Dependencies:
%
%       1) The Toolbox is supported only on 64-bit MATLAB 2013a - 2015a
%          and on GNU Octave 3.6.4 and later, on Linux, Mac OS X, and Windows.
%
%       2) On 64-bit MATLAB you will need the Statistics Toolbox (till
%          2014b) or Statistics and Machine Learning Toolbox (since 2015a).
%
%
% Written by Marcus Vollmer, January 15, 2015.
% Adapted code from the 2014 PhysioNet/CinC Challenge function by Ikaro
% Silva, December 10, 2013. 
%
% Last Modified: May 27, 2015
% Version 0.5
%
%
%endOfHelp

%Set default parameter values
    Fs = 1000;   
    RR = .6;% Estimated mean heart rate    
    wl_tma = ceil(.2*Fs);	%ECG references
    %     wl_tma = ceil(Fs);      %BP channel
    threshold = .5;
    downsampling = 80;
    debugmode = 0;  
    pct=.25;
    R=.4;
    Beat_min=50;
    Beat_max=220;
    debug_file = 'record_debug.mat';
 
%Set parameter values by argument
record = varargin{1};
if nargin>1
    inputs = {'Fs','RR','wl_tma','threshold','downsampling','debugmode','pct','R','Beat_min','Beat_max','debug_file'};
    for n=2:2:nargin
        tmp = find(strcmp(varargin{n},inputs), 1);
        if(~isempty(tmp))
            eval([inputs{tmp} ' = varargin{n+1};'])
        else
            error(['''' varargin{n} ''' is not an accepted input argument.'])
        end
    end
end

% Downsampling by filter (moving average)
    factor = max([1 floor(Fs/downsampling)]);                     
    Fs = Fs/factor;
    sig = filter((1/factor)*ones(1,factor),1,record(:));
    sig = sig(1:factor:end,:);
    clear record

% signal with NaNs, linear interpolation
    [pos_r,pos_c] = find(isnan(sig));
    if ~isempty(pos_r) 
        for c=unique(pos_c)'
            r = pos_r(pos_c==c);
            if length(r)>1
                tmp_pos = [r(1) r(diff(r)>1) r(end)];
            else
                tmp_pos = [r(1) r(end)];
            end
            for nan_num=1:length(tmp_pos)-1   
                if tmp_pos(nan_num)>1 && tmp_pos(nan_num+1)<size(sig,1)
                    sig(tmp_pos(nan_num):tmp_pos(nan_num+1),c) = interp1([tmp_pos(nan_num)-1 tmp_pos(nan_num+1)+1],sig([max(1,tmp_pos(nan_num)-1) min(tmp_pos(nan_num+1)+1,size(sig,1))],c),tmp_pos(nan_num):tmp_pos(nan_num+1));
                end
            end
        end
    end

% High-pass filtering and standardization
    [sig_tma,tma] = tma_filter(sig,wl_tma,pct);	% trimmed moving average
    sig_tmzscore  = zscore(sig_tma);            % zscore filter

    
% Estimated QT interval
    QT = .228; %assumed, better: estimate
    QT = QT+0.154*RR; %Alex Sagie et al.
    RT = .85*QT; %R peak to T peak, estimate should be verified
    

% Compute annotations
    wl_we = ceil(RT*Fs); %ECG references

    if debugmode>=1
        [Ann,valid,tmpmin,tmpmax,wmin,wmax,range,thr] = mvqrs_ann(sig_tmzscore,Fs,wl_we,Beat_min,Beat_max,threshold,R);
    else
        [Ann,valid] = mvqrs_ann(sig_tmzscore,Fs,wl_we,Beat_min,Beat_max,threshold,R);  
    end
        
    Ann = Ann(Ann>0 & Ann<=size(sig,1));
    myAnn = unique(sort(Ann));
    myAnn(diff(myAnn)<.05*Fs)=[];

    myAnn = myAnn*factor;  


% Open figures in debugmode
    if debugmode>=1
        Ann = myAnn/factor;
        tolerance = floor(.15*Fs);
        save record_debug.mat record sig tma sig_tmzscore tmpmin tmpmax wmin wmax valid range thr Ann Fs factor tolerance Beat_min Beat_max threshold;
        movefile('record_debug.mat',debug_file);
    end

end