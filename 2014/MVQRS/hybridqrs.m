function [myAnn, allAnn] = hybridqrs(varargin)
%
% myAnn = hybridqrs(record,...)
%
% Octave-compatible code for merging Annotation files
%
% Required Parameters:
%
% record
%       Name of the Annotation files. File extentions must be .qrs_NAME,
%       where NAME are acronyms of considered annotation files.
%
% The record can be followed by parameter/value pairs to specify
% additional properties of the QRS detection.
%
% Fs
%       The sampling frequency in Hz. (default: 1000)
% threshold
%       The threshold factor is a value between 0 and 1. The higher the
%       factor the less annotations. (default: 0.5)
% detectors
%       An array of strings specifing annotation files, which will be
%       merged. (default: empty, which means that all files in terms
%       "record.qrs_* files" will be used)
%
% This function has the output argument myAnn and allAnn.
% myAnn is a vector containing hybrid annotation for each detected beat in
% the record.
% allAnn is a matrix with annotations of several beat detectors.
%
% Dependencies:
%
%       1) The Toolbox is supported only on 64-bit MATLAB 2013a - 2015a
%          and on GNU Octave 3.6.4 and later, on Linux, Mac OS X, and Windows.
%
%       2) You will need the Statistics Toolbox (till 2014b) or Statistics
%          and Machine Learning Toolbox (since 2015a). 
%
%       3) Physionet WFDB Toolbox for loading annotation files.
%
%
% Written by Marcus Vollmer, November 12, 2015.
% Last Modified: December 16, 2015
% Version 0.2
%
%
%endOfHelp

%Set default parameter values
    Fs = 1000;   
    threshold = .5;
    detectors = '';
    output_file = 'qrs_hybrid';
    debugmode = 0;  
    debug_file = 'hybrid_debug.mat';
 
%Set parameter values by argument
record = varargin{1};
if nargin>1
    inputs = {'Fs','threshold','detectors','output_file','debugmode','debug_file'};
    for n=2:2:nargin
        tmp = find(strcmp(varargin{n},inputs), 1);
        if(~isempty(tmp))
            eval([inputs{tmp} ' = varargin{n+1};'])
        else
            error(['''' varargin{n} ''' is not an accepted input argument.'])
        end
    end
end
  
%Load annotations
    if isempty(detectors)
        detectors = dir([record '.qrs_*']);
        detectors = {detectors.name}.';
    else
        detectors = strcat([record '.'],detectors);
    end
    m = 0;
    for i=1:length(detectors)
        tmpAnn = rdann(record,detectors{i}(length(record)+2:end));                
        allAnn.(detectors{i}(length(record)+2:end)) = tmpAnn;
        m = max(m,max(tmpAnn));
    end
    Annotation = zeros(Fs*ceil(m/Fs),length(detectors));
    for i=1:length(detectors)
        Annotation(allAnn.(detectors{i}(length(record)+2:end)),i) = 1;
    end    
    
    
% Evaluate Cross-QRS-Delays using Cross Correlation of filtered annotations
    steps = ceil(2*.15*Fs);
    A_filt = filter(normpdf(-4:(8/steps):4)/normpdf(0),1,Annotation);
    [cr,lgs] = xcorr(A_filt,ceil(Fs),'coeff');
    if debugmode>=2
        figure        
    end
    lags = NaN(length(detectors),length(detectors));
    for row = 1:length(detectors)
        for col = 1:length(detectors)
            nm = length(detectors)*(row-1)+col;

            [~,locs]  = findpeaks(abs(cr(:,nm)));
            loc = locs(cr(locs,nm)>.5);
            lagDiff = lgs(loc);
            loc = loc(abs(lagDiff)==min(abs(lagDiff)));
            lagDiff = lagDiff(abs(lagDiff)==min(abs(lagDiff)));
            if ~isempty(lagDiff)
                lags(row,col) = max(lagDiff);
            end

            if debugmode>=2
                subplot(length(detectors),length(detectors),nm)
                stem(lgs,cr(:,nm),'.')
                title(sprintf('c_{%d%d}',row,col))
                ylim([0 1.3])
                hold on
                plot(lagDiff,cr(loc,nm),'v')
                text(lagDiff,cr(loc,nm)+.05,num2str(lagDiff),...
                    'HorizontalAlignment','center','VerticalAlignment','bottom')
                grid on
                hold off
                if row==1
                    title(detectors{col}(length(record)+2:end),'Interpreter','none')
                end
                if col==1
                    set(get(gca,'YLabel'),'String',detectors{row}(length(record)+2:end),'Interpreter','none')
                end
            end
        end
        if debugmode>=2
            axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized','clipping','off');
            text(0.5, .975,['\bf\fontsize{12}record: ' record],'HorizontalAlignment','center','VerticalAlignment','top')
        end
    end
    discard_detectors = sum(isnan(tril(lags)))>=.5*length(detectors);

    lags_neu = lags;
    lags_neu((abs(lags)>.15*Fs)) = NaN;
    l = length(detectors)-sum(isnan(lags_neu));
    delays = -round(nansum(lags_neu)./l);
    D = delays;
    
    for k=find(isnan(delays) & ~discard_detectors)   
        D(k) = round(nanmean(lags(k,~isnan(delays))+delays(~isnan(delays))));
    end

% Correction of algorithmic delays
    for k=find(~isnan(D) & ~discard_detectors)
        Annotation(:,k) = 0;
        Annotation(max(1,allAnn.(detectors{k}(length(record)+2:end))-D(k)),k) = 1;
    end


%% VALIDATION OF SIGNAL ANNOTATIONS 
% find the best set of heart beat annotations and generate a hybrid
% annotation file and check candidates

Annotation_full = Annotation(:,~discard_detectors);
valid = ones(size(Annotation_full));
ValidSignals = 1:size(Annotation_full,2);
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
    % More than 20% candidates than merged annotations:   
                [~,idx]= max(sum(accepted_RR_pct(:,ts_section(i)+1:min([size(quality_RR,2) ts_section(i+1)])),2));            
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
    
% merge close annotations which may occur in overlapping sections 
    myAnn = unique(sort(myAnn));
    myAnn(diff(myAnn)<.05*Fs)=[];


% Plot annotations
    if debugmode>=2
        figure;
        tmp = repmat(myAnn,1,3)';
        plot(tmp(:),repmat([0;length(detectors)+1;0],length(myAnn),1),'k','linewidth',2)       
        hold on;
        for i=1:length(detectors)
            plot(allAnn.(detectors{i}(length(record)+2:end)),i*ones(size(allAnn.(detectors{i}(length(record)+2:end)))),'o','markersize',10);
        end
        ylim([0,i+1])
        xlim([0 10*Fs])
        set(gca,'YTick',1:i,'YTickLabel',detectors,'TickLabelInterpreter','none')
        set(gca,'XTick',0:Fs:myAnn(end))
    end
       

% Using the Toolbox's wrann function, to write annotations. 
    % see Issue ID52:
    % wrann(record,output_file,myAnn);
    if size(myAnn,1)<2500
        wrann(record,output_file,myAnn);
    else
        wrann(record,output_file,myAnn(1:2000)); 
        for ann_num=2:ceil(size(myAnn,1)/2000)
            wrann(record,[output_file num2str(ann_num)],myAnn((ann_num-1)*2000+1:min(ann_num*2000,size(myAnn,1))));
            mrgann(record,output_file,[output_file num2str(ann_num)],'hybrid_merge')
            copyfile([record '.hybrid_merge'],[record '.' output_file])
            delete([record '.' output_file num2str(ann_num)])
        end
        delete([record '.hybrid_merge'])
    end
    % For more information, run 'help wrann'.


% Open figures in debugmode
    if debugmode>=1
        Annotation = Annotation_full;
%         valid = valid_full; 
        save hybrid_debug_tmp.mat myAnn allAnn Annotation discard_detectors D Annotation_full accepted_RR_pct accepted_RR_count quality_RR quality quality2;
        movefile('hybrid_debug_tmp.mat',debug_file);
    end

end