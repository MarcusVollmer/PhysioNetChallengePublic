function [pacemaker,RR,delay_ECG_val] = mvqrs_pacemakerdetection(recordName,ValidSignalsSetting,Fs,factor,Beat_max,Beat_min)
%
% [pacemaker,RR,delay_ECG_val] = mvqrs_pacemakerdetection(recordName,ValidSignalsSetting,Fs,factor,Beat_max,Beat_min)
%
% Estimates mean heart rate & check for pacemaker.
%
%
%
% Written by Marcus Vollmer, 2015
% Last Modified: February 11, 2015
% Version 0.1
%
%endOfHelp

pacemaker = 0;
pacemaker_ecg_ch_name = '';
delay_ECG_val = NaN;

ecg_ch = find(ValidSignalsSetting==0);
bp_ch = find(ValidSignalsSetting==2);

if length(ecg_ch)+length(bp_ch)>0
    RR = zeros(length(ecg_ch)+length(bp_ch),1);

    for i=1:length(ecg_ch)
        try
            gqrs(recordName,[],1,ecg_ch(i),1,['qrs' num2str(i)])
        catch
        end
        try
            Ann = rdann(recordName,['qrs' num2str(i)]);
            RRinterval = diff(Ann);
            if size(RRinterval,1)>40
                relRR2 = 2*(RRinterval(1:end-2)-RRinterval(3:end))./(RRinterval(1:end-2)+RRinterval(3:end));
                RRinterval_medfilter = medfilt1(RRinterval,20);
                RRinterval_medfilter = RRinterval_medfilter(11:end-9);
                RRinterval_filter = filter(ones(20,1)/20,1,RRinterval);
                RRinterval_filter = RRinterval_filter(20:end);
                
                pct = min([iqr(RRinterval_medfilter)/max([2 ceil(0.01*nanmedian(RRinterval_medfilter))]) ...
                    iqr(RRinterval_filter)/max([2 ceil(0.01*nanmedian(RRinterval_filter))])]) +...
                    iqr(relRR2)/0.025;
                    
                if pct<=2
                    pacemaker = 1;
                    pacemaker_ecg_ch_name = ['qrs' num2str(i)];
                end
            end
        catch
        end
    end

    for i=1:length(bp_ch)
        wabp(recordName,[],[],0,bp_ch(i))
        try
            Ann = rdann(recordName,'wabp');
            RRinterval = diff(Ann);
            if size(RRinterval,1)>40%sum(~isnan(RRinterval))>30 && sum(isnan(RRinterval))/size(RRinterval,1)<.2
                RR(i+length(ecg_ch)) = nanmedian(RRinterval);
                % pacemaker check
                if 60*Fs*factor/RR(i+length(ecg_ch))<120
                    % pacemaker settings with fixed-rate stimulation: Medtronic 85 bpm, Guidant 100 bpm, St. Jude 98 bpm, Biotronik 90 bpm 
                    relRR2 = 2*(RRinterval(1:end-2)-RRinterval(3:end))./(RRinterval(1:end-2)+RRinterval(3:end));
                    RRinterval_medfilter = medfilt1(RRinterval,20);
                    RRinterval_medfilter = RRinterval_medfilter(11:end-9);
                    RRinterval_filter = filter(ones(20,1)/20,1,RRinterval);
                    RRinterval_filter = RRinterval_filter(20:end);

                    pct = min([iqr(RRinterval_medfilter)/max([2 ceil(0.01*nanmedian(RRinterval_medfilter))]) ...
                        iqr(RRinterval_filter)/max([2 ceil(0.01*nanmedian(RRinterval_filter))])]) +...
                        iqr(relRR2)/0.02;
                    if pct<=2
                        pacemaker = 1;
                        if ~isempty(ecg_ch)
                            if isempty(pacemaker_ecg_ch_name)
                                Ann_ECG = rdann(recordName,'qrs1');
                            else
                                Ann_ECG = rdann(recordName,pacemaker_ecg_ch_name);
                            end
                            if size(Ann_ECG,1)<1.5*size(Ann,1)
                                r = [Ann_ECG; Ann];
                                c = [zeros(length(Ann_ECG),1); ones(length(Ann),1)];
                                [~,idx] = sort(r);
                                Pos = [r(idx) c(idx)];                                
                                Posdiff = Pos(2:end,:)-Pos(1:end-1,:);

                                DELAY_tmp = Posdiff(Posdiff(:,2)==1,1)/(Fs*factor);  
                                DELAY_tmp(:,2) = Pos(Posdiff(:,2)==1);            
                                %figure; plot(DELAY_tmp(:,2),DELAY_tmp(:,1))
                                if median(DELAY_tmp(:,1))<.25
                                    delay_ECG_val = -median(DELAY_tmp(:,1))/2;
                                end
                            end
                        end
                    end
                end
            end
        catch
        end
    end  


    tmp = RR(RR>0)/(Fs*factor);
    RR = .5;        
    if ~isempty(tmp)
        tmp = tmp(abs(tmp-median(tmp))/median(tmp)<.2);
        if ~isempty(tmp)  
            RR = median([60/Beat_max mean(tmp) 60/Beat_min]); 
        end
    end
else
    RR = .5;
end
