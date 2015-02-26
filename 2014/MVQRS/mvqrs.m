function mvqrs(varargin)
%
% mvqrs(recordName,...)
%
% Octave-compatible code for QRS detection
%
% Required Parameters:
%
% recordName
%       String specifying the name of the record to process.  Do not include
%       the '.dat' or '.hea' suffix in recordName.
%
% The record name can be followed by parameter/value pairs to specify
% additional properties of the QRS detection.
%
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
% debugmode
%       An true-false variable which enables the debug mode of the
%       calculation. Various figures will show false negative and false
%       positive annotations. Therefor a 'atr','ari' or 'ecg' file must
%       exists in the same folder as the record file. (default: 0)
%       By now: only available using Matlab.
%
% This function has no output arguments, but it writes an annotation file
% named "recordName.mvqrs" in the current directory, containing a
% annotation for each detected beat in the record.
%
% Dependencies:
%
%       1) This function requires the WFDB Toolbox 0.9.7 and later for
%          MATLAB and Octave. 
%          For information on how to install the toolbox, please see:
%              http://www.physionet.org/physiotools/matlab/wfdb-app-matlab/
%
%       2) The Toolbox is supported only on 64-bit MATLAB 2013a - 2015a
%          and on GNU Octave 3.6.4 and later, on Linux, Mac OS X, and Windows.
%
%       3) On 64-bit MATLAB you will need the Statistics Toolbox (till
%          2014b) or Statistics and Machine Learning Toolbox (since 2015a).
%
%
% This file was written as part of a sample entry for the 2014 Challenge;  for
% the complete sample entry including all other required components, see
%    http://physionet.org/challenge/2014/octave/
% For a description of the Challenge and how to participate in it, see
%    http://physionet.org/challenge/2014/
%
% Written by Marcus Vollmer, January 15, 2015.
% Adapted code from the 2014 PhysioNet/CinC Challenge function by Ikaro
% Silva, December 10, 2013. 
%
% Last Modified: February 26, 2015
% Version 0.5
%
% %Examples:
% mvqrs('mitdb/200')
%
% mvqrs('mitdb/200',...
%     'threshold',.6,...
%     'sec',20,...
%     'downsampling',150,...
%     'debugmode',2)
%
%endOfHelp

%Get record name 
if(~isempty(varargin))
    eval('recordName=varargin{1};')
    if ~exist([recordName '.dat'],'file')
        error(['Undefined record name or folder ''' recordName '.dat''.' ])
    end    
else
    error('Please specify the record name as the first argument.')
end  

%Set default parameter values
    threshold = .5;
    sec = [];
    downsampling = 80;
    debugmode = 0;  
    pct=.25;
    R=.4;
    Beat_min=50;
    Beat_max=220;
    use_gqrs = 0;
    use_wabp = 0;
    use_gqpost = 0;
    debug_file = [recordName '_debug.mat'];
    output_file = 'mvqrs';
 
%Set parameter values by argument
if nargin>1
    inputs = {'threshold','sec','downsampling','debugmode','pct','R','Beat_min','Beat_max','use_gqrs','use_wabp','use_gqpost','debug_file','output_file'};
    for n=2:2:nargin
        tmp = find(strcmp(varargin{n},inputs), 1);
        if(~isempty(tmp))
            eval([inputs{tmp} ' = varargin{n+1};'])
        else
            error(['''' varargin{n} ''' is not an accepted input argument.'])
        end
    end
end


% debug mode for matlab user
    v = version;
    octave = 0;
    if isempty(strfind(v,'R'))
        octave = 1;
% Using the WFDB Toolbox's wfdbloadlib function, initialize the Toolbox
% configuration to the default values.  (This is required for Octave.)       
        [~,config] = wfdbloadlib;
        if(config.inOctave)
            crash_dumps_octave_core(0);
        end
    end


% Generate the 100m.mat and 100m.hea files from the *.dat and *.hea files
% Determine which channels can be used for annotation
    if exist([recordName 'm.mat'],'file')~=2
        wfdb2mat(recordName);
    end
    if ~isempty(sec)    
        [ValidSignals,ValidSignalsSetting,Fs,gain,val] = ...
            mvqrs_signaltest(recordName,threshold,sec,downsampling); 
    else
        [val,Fs,gain,description,ecg_idx,bp_idx] = ...
            mvqrs_loadrecord(recordName);
        ValidSignals = 1:size(val,2);
        ValidSignalsSetting = 3*ones(1,size(val,2));
        ValidSignalsSetting(ecg_idx) = 0;
        ValidSignalsSetting(bp_idx) = 2;
    end
%     delete([recordName 'm.mat'])
%     delete([recordName 'm.hea'])


% % Downsampling by omission
%     factor = floor(Fs/downsampling);                     
%     Fs = Fs/factor;
%     sig = val(1:factor:end,ValidSignals); % omit information
%     sig = sig./repmat(gain(ValidSignals),size(sig,1),1);
%     clear val
% Downsampling by filter (moving average)
    factor = max([1 floor(Fs/downsampling)]);                     
    Fs = Fs/factor;
    sig = filter((1/factor)*ones(1,factor),1,val(:,ValidSignals));
    sig = sig(1:factor:end,:);
    if ~isempty(sec)
        sig = sig./repmat(gain(ValidSignals),size(sig,1),1);
    end
    clear val    


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
    sig_tmzscore = zeros(size(sig,1),length(ValidSignals)); 
    tma = zeros(size(sig,1),length(ValidSignals));
    for i=1:length(ValidSignals) 
        switch ValidSignalsSetting(i)   % set signal class specific parameter
            case 0; wl_tma = ceil(.2*Fs);	%ECG references
            case 1; wl_tma = ceil(.2*Fs);	%only ECG setting, not ECG reference
            case 2; wl_tma = ceil(Fs);      %BP channel
            case 3; wl_tma = ceil(Fs);      %only BP setting, not BP channel
        end

        [sig_tma,tma(:,i)] = tma_filter(sig(:,i),wl_tma,pct);	% trimmed moving average
        sig_tmzscore(:,i) = zscore(sig_tma);             % zscore filter
        % [sig_tmzscore,~,tmstd] = tmzscore_filter(sig,wl_tma,pct);      % trimmed moving standardization
        % 
        % % set sig_tmzscore to zero if tmstd under signal_precision
        % signal_precision = min(setdiff(unique(diff(sort(sig))),0))
        % sig_tmzscore(tmstd<10*signal_precision) = 0;  
    end
          
    
% Initialization
    Annotation = zeros(size(sig,1),length(ValidSignals));
    valid = ones(size(sig,1),length(ValidSignals)); 
    DELAY = zeros(size(sig,1),length(ValidSignals));
    if debugmode>=1
        tmpmin = zeros(size(sig,1),length(ValidSignals));
        tmpmax = zeros(size(sig,1),length(ValidSignals));
        wmin = zeros(size(sig,1),length(ValidSignals));
        wmax = zeros(size(sig,1),length(ValidSignals));
        range = zeros(size(sig,1),length(ValidSignals));
        thr = zeros(size(sig,1),length(ValidSignals));
    end

    
% Estimate mean heart rate & check for pacemaker
[pacemaker,RR,delay_ECG_val] = mvqrs_pacemakerdetection(recordName,ValidSignalsSetting,Fs,factor,Beat_max,Beat_min);


% Estimated QT interval
    QT = .228; %assumed, better: estimate
    QT = QT+0.154*RR; %Alex Sagie et al.
    RT = .85*QT; %R peak to T peak, estimate should be verified
   
    
% expected algorithmic delay for signals
    delay_BP = 0.2;
    delay_ECG = 0;

% Set pacemaker delay
    if pacemaker==1
        if isnan(delay_ECG_val)
            delay_ECG = -0.15; %assume 150ms ventriculoatrial interval
        else
            delay_ECG = delay_ECG_val;
        end
        if debugmode>=1
            fprintf(['pacemaker detected in ' recordName '. ECG_delay is ' num2str(delay_ECG,'%0.2f') '\n'])
        end        
    end

% Compute annotations for each signal
    for i=1:length(ValidSignals) 
        switch ValidSignalsSetting(i)   % set channel specific parameter
            case 0; delay = delay_ECG; wl_we = ceil(RT*Fs); %ECG references
            case 1; delay = 0; wl_we = ceil(RT*Fs); %only ECG setting, not ECG reference
            case 2; delay = 0; wl_we = [ceil(1.5*RT*Fs) ceil(RT*Fs)]; %BP channel
            case 3; delay = 0; wl_we = [ceil(1.5*RT*Fs) ceil(RT*Fs)]; %only BP setting, not BP channel
        end

        if ValidSignalsSetting(i)==0 && use_gqrs==1
            gqrs(recordName,[],1,i)
            Ann = floor(rdann(recordName,'qrs')/factor);
        elseif ValidSignalsSetting(i)==2 && use_wabp==1
            wabp(recordName,[],[],0,i);
            try
                Ann = floor(rdann(recordName,'wabp')/factor);
            catch
                if debugmode>=1
                    [Ann,valid(:,i),tmpmin(:,i),tmpmax(:,i),wmin(:,i),wmax(:,i),range(:,i),thr(:,i)] = mvqrs_ann(sig_tmzscore(:,i),Fs,wl_we,Beat_min,Beat_max,threshold,R);
                else
                    [Ann,valid(:,i)] = mvqrs_ann(sig_tmzscore(:,i),Fs,wl_we,Beat_min,Beat_max,threshold,R);  
                end                
            end
        else
            if debugmode>=1
                [Ann,valid(:,i),tmpmin(:,i),tmpmax(:,i),wmin(:,i),wmax(:,i),range(:,i),thr(:,i)] = mvqrs_ann(sig_tmzscore(:,i),Fs,wl_we,Beat_min,Beat_max,threshold,R);
            else
                [Ann,valid(:,i)] = mvqrs_ann(sig_tmzscore(:,i),Fs,wl_we,Beat_min,Beat_max,threshold,R);  
            end
        end
        
        Ann = Ann-round(delay*Fs); % delay correction
        Ann = Ann(Ann>0 & Ann<=size(sig,1));
        Annotation(Ann,i)=1;
    end
    

% static delay for annotations
    DELAY(:,ValidSignalsSetting==0) = delay_ECG; %fixed delay for ECG signals
    DELAY(:,ValidSignalsSetting==2) = delay_BP; %fixed delay for BP signals

% apply continuous delay
    [Annotation,DELAY] = mvqrs_delaycontrol(Annotation,Fs,ValidSignalsSetting,delay_BP,DELAY);


%% VALIDATION OF SIGNAL ANNOTATIONS 
% find the best set of signals and generate a multivariate testAnn file and
% candidates

Annotation_full = Annotation;
valid_full = valid;   
myAnn = [];

% From now on the annotations will be checked of all valid channels
% every 30 seconds. Various expected beat positions have to be combined.
    Fm = round(Fs*30);
    
% get signal quality using RR intervals
    [accepted_RR_pct,accepted_RR_count] = mvqrs_RRquality(Annotation_full,ValidSignals,Fm);

% take a subset using the RR quality parameter
    quality_RR = (accepted_RR_pct.^2).*accepted_RR_count;
    quality = quality_RR>=repmat(.8*max(quality_RR),length(ValidSignals),1);
    quality2 = quality;
    quality2(accepted_RR_pct<.9) = 0;

% take signals with more than 90% acceptance rate if there are enough
% signals 
tmp = find(sum(quality)>2 & sum(quality2)>=2);
quality(:,tmp) = quality2(:,tmp);

ts_section = [0 find(sum(abs(diff(quality,[],2)))~=0) ceil(size(Annotation_full,1)/Fm)];
ts_overlap = round(Fs*2);
for i=1:length(ts_section)-1
    subset = find(quality(:,min([size(quality,2) ts_section(i+1)])));

    % get time series section with overlap   
    Annotation = Annotation_full(max([1 (ts_section(i))*Fm+1-ts_overlap]):min([(ts_section(i+1))*Fm+ts_overlap size(Annotation_full,1)]),subset);
    valid = valid_full(max([1 (ts_section(i))*Fm+1-ts_overlap]):min([(ts_section(i+1))*Fm+ts_overlap size(Annotation_full,1)]),subset); 

% A multivariate procedure combines several beat positions and export
% candidates.  
    if length(subset)>1
        [testAnn,candidates] = mvqrs_multivariateAnn(Annotation,valid,ceil(.15*Fs)); 
    else
        testAnn = find(Annotation); candidates=[];
    end
    
% Check candidates within each section
    if ~isempty(candidates)
        if length(candidates)>.2*length(testAnn) 
% More than 50% candidates than merged annotations:   
            if length(candidates)>.5*length(testAnn)
                % case (i) paced signal: take ann of BP signal
                idx = bp_idx(find(ismember(bp_idx,subset)));
                if length(idx)>1
                    [~,tmp] = max(sum(accepted_RR_pct(idx,ts_section(i)+1:min([size(quality_RR,2) ts_section(i+1)])),2));
                    idx = idx(tmp);
                end
            else               
                % case (ii) signal with extrasystoles: take ann of ECG signal
                idx = ecg_idx(find(ismember(ecg_idx,subset)));
                if length(idx)>1
                    [~,tmp] = max(sum(accepted_RR_count(idx,ts_section(i)+1:min([size(quality_RR,2) ts_section(i+1)])),2));
                    idx = idx(tmp);
                end
            end            
            if isempty(idx)
                [~,idx]= max(sum(accepted_RR_pct(:,ts_section(i)+1:min([size(quality_RR,2) ts_section(i+1)])),2)); 
            end            
            testAnn = find(Annotation(:,subset==idx));
        else
            candidates_signal = Annotation(candidates(:,1),:)*[1:size(Annotation,2)]';
            cand_ch = tabulate(candidates_signal);
% check candidates if the number of candidates is less than 30% of accepted
% beats or if there is no dominant channel which generates more than 70% of
% all candidates (false negative generating signal)   
            if sum(cand_ch(:,2)>.7*size(candidates,1))==0 || size(candidates,1)<.3*size(testAnn,1)

% check candidates for remature ventricular contraction (PVC), ventricular
% premature contraction (VPC), ventricular premature beat (VPB), or
% ventricular extrasystole (VES)              
                % needs to be done
                testAnn = mvqrs_checkcandidates(testAnn,Fs,candidates);
            end
        end
    end
    % remove annotation in overlap part
    if i==1
        if i==length(ts_section)-1
            myAnn = testAnn;
        else
            myAnn = [myAnn; testAnn(testAnn<=ts_section(i+1)*Fm)];
        end
    elseif i==length(ts_section)-1
        myAnn = [myAnn; testAnn(testAnn>=ts_overlap)+ts_section(i)*Fm-ts_overlap]; 
    else
        myAnn = [myAnn; testAnn(testAnn>=ts_overlap & testAnn<=diff(ts_section([i i+1]))*Fm+ts_overlap)+ts_section(i)*Fm-ts_overlap]; 
    end
    
end
% % merge close annotations which may occur in overlapping sections 
myAnn = unique(sort(myAnn));
myAnn(diff(myAnn)<.05*Fs)=[];

Annotation = Annotation_full;
valid = valid_full; 
   
myAnn = myAnn*factor;  

% Using the Toolbox's wrann function, write the time-shifted annotations
% into recordName.qrs.

% see Issue ID52:
% wrann(recordName,output_file,myAnn);
if size(myAnn,1)<2500
    wrann(recordName,output_file,myAnn);
else
    wrann(recordName,output_file,myAnn(1:2000)); 
    for ann_num=2:ceil(size(myAnn,1)/2000)
        wrann(recordName,[output_file num2str(ann_num)],myAnn((ann_num-1)*2000+1:min(ann_num*2000,size(myAnn,1))));
        mrgann(recordName,output_file,[output_file num2str(ann_num)],'mvqrs_merge')
        copyfile([recordName '.mvqrs_merge'],[recordName '.' output_file])
        delete([recordName '.' output_file num2str(ann_num)])
    end
    delete([recordName '.mvqrs_merge'])
end
% For more information, run 'help wrann'.

if use_gqpost==1
    gqpost(recordName,output_file);
%     delete([recordName '.' output_file])
%     rename([recordName '.gqp'],[recordName '.' output_file])
end


% Open figures in debugmode
    if debugmode>=1
        Ann = myAnn/factor;
        tolerance = floor(.15*Fs);
        save record_debug.mat recordName sig tma sig_tmzscore tmpmin tmpmax wmin wmax valid range DELAY thr Ann Annotation Fs factor tolerance ValidSignals Beat_min Beat_max threshold quality_RR accepted_RR_pct accepted_RR_count quality ts_section pacemaker;
        movefile('record_debug.mat',debug_file);

        if debugmode==2
            if octave==0
                mvqrs_debugmode_gui(debug_file)
            else
            	mvqrs_debugmode(debug_file)
            end
        end
    end

end