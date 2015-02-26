function [Annotation,DELAY] = mvqrs_delaycontrol(Annotation,Fs,ValidSignalsSetting,delay_BP,DELAY)
%
% [Annotation,DELAY] = mvqrs_delaycontrol(Annotation,Fs,ValidSignalsSetting,DELAY)
%
% Applies the beat annotation correction using a reference ECG electrode.
% Necessary for blood pressure like signals with expectant pulse transition
% time.
%
% Required Parameters:
%
% Annotation
%       A NxM matrix with {0,1} values with expected beat postions of M
%       different signals.
% Fs
%       The sampling frequency in Hz.
% ValidSignalsSetting
%       A 1xM vector with the settings number of each signal.
% delay_BP
%       A double-precision value which specifies static delay for blood
%       pressure like signals.
%
%
% Written by Marcus Vollmer, 2015
% Last Modified: February 23, 2015
% Version 0.3
%
%endOfHelp



if sum(ValidSignalsSetting>0)>0
    ecg_idx = find(ValidSignalsSetting==0);
    ecgchannel = [];

    if ~isempty(ecg_idx)
        % Determine the best reference ECG electrode
        [accepted_RR_pct,accepted_RR_count] = mvqrs_RRquality(Annotation,ecg_idx,[]);
        quality_RR = accepted_RR_pct.*accepted_RR_count;
        [~,best] = max(quality_RR);

        if accepted_RR_pct(best(1))>.8
            ecgchannel = ecg_idx(best(1));
        end
    end

    if ~isempty(ecgchannel)
        for i=find(ValidSignalsSetting~=0)
            [r,c] = find(Annotation(:,[ecgchannel i])==1);
            [~,idx] = sort(r);
            Pos = [r(idx) c(idx)];
            Posdiff = Pos(2:end,:)-Pos(1:end-1,:);

            % computes the delay from ECG reference to the other signal
            DELAY_tmp = Posdiff(Posdiff(:,2)==1,1)/Fs;  
            DELAY_tmp(:,2) = Pos(Posdiff(:,2)==1);            
            %plot(DELAY_tmp(:,2),DELAY_tmp(:,1))

            % use a +-10 second interval to interpolate the delay
            iv = round(10*Fs);

            AnnPos_All = find(sum(Annotation,2)>0);
            AnnPos = find(Annotation(:,i)==1);
            for j=1:size(AnnPos,1)
                % search for all computed delays inside the interval
                tmp_wd = DELAY_tmp(DELAY_tmp(:,2)<AnnPos(j)+iv & DELAY_tmp(:,2)>AnnPos(j)-iv,1);
                if ~isempty(tmp_wd)
                    % take the median as the continuous delay
                    cont_delay_BP = median(tmp_wd);
                else
                    % increase the interval if tmp_wd is empty
                    tmp_wd = DELAY_tmp(DELAY_tmp(:,2)<AnnPos(j)+2*iv & DELAY_tmp(:,2)>AnnPos(j)-2*iv,1);
                    if ~isempty(tmp_wd)
                        % take the median as the continuous delay if <400ms
                        if median(tmp_wd)<.4
                            cont_delay_BP = median(tmp_wd);
                        else
                            cont_delay_BP = delay_BP;
                        end
                    else
                        % otherwise take the assumed delay of blood
                        % pressure signals (defined above).
                        cont_delay_BP = delay_BP;
                    end
                end
                Annotation(AnnPos(j),i) = 0;
                if AnnPos(j)>round(cont_delay_BP*Fs)
                    Annotation(AnnPos(j)-round(cont_delay_BP*Fs),i) = 1;
                    DELAY(AnnPos(j),i) = cont_delay_BP;
                end
            end

            AnnPos_Other = setdiff(AnnPos_All,AnnPos);
            for j=1:size(AnnPos_Other,1)
                tmp_wd = DELAY_tmp(DELAY_tmp(:,2)<AnnPos_Other(j)+iv & DELAY_tmp(:,2)>AnnPos_Other(j)-iv,1);
                if ~isempty(tmp_wd)
                    cont_delay_BP = median(tmp_wd);
                else
                    tmp_wd = DELAY_tmp(DELAY_tmp(:,2)<AnnPos_Other(j)+2*iv & DELAY_tmp(:,2)>AnnPos_Other(j)-2*iv,1);
                    if ~isempty(tmp_wd)
                        cont_delay_BP = median(tmp_wd);
                    else
                        cont_delay_BP = delay_BP;
                    end
                end

                if AnnPos_Other(j)>round(cont_delay_BP*Fs)
                    DELAY(AnnPos_Other(j),i) = cont_delay_BP;
                end
            end 
        end
    else
        % No trusted ecg signal
        % perform static adjustment
        for i=find(ValidSignalsSetting~=0)
            AnnPos = find(Annotation(:,i)==1);
            DELAY(AnnPos,i) = delay_BP;
            Annotation(AnnPos,i) = 0;
            AnnPosN = AnnPos-round(delay_BP*Fs);
            Annotation(AnnPosN(AnnPosN>0),i) = 1;
        end
        
    end
end    

