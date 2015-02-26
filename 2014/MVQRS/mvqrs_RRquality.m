function [accepted_RR_pct,accepted_RR_count] = mvqrs_RRquality(Annotation,ValidSignals,Fm)
%
% [accepted_RR_pct,accepted_RR_count] = mvqrs_RRquality(Annotation_full,ValidSignals,Fm)
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
% ValidSignals
%       A boolean 1xM vector with the signal to be used.
% Fm
%       The define the section length.
%
%
% Written by Marcus Vollmer, 2015
% Last Modified: February 19, 2015
% Version 0.1
%
%endOfHelp


if isempty(Fm)
    Fm = size(Annotation,1);
end

% From now on the annotations will be checked of all valid channels
% minutely. Various expected beat positions have to be combined.
accepted_RR_count = zeros(length(ValidSignals),max([0 (floor(size(Annotation,1)/Fm)-1)])+1);
accepted_RR_pct = zeros(length(ValidSignals),max([0 (floor(size(Annotation,1)/Fm)-1)])+1);

for ts_section=0:(max([1 floor(size(Annotation,1)/Fm)])-1)
  
    % valid = valid_full(((ts_section*Fm)+1):end,:); 
    % valid = valid_full(((ts_section*Fm)+1):((ts_section+1)*Fm),:);              
    % quality_valid = sum(valid)/size(valid,1);

    % check whether a channel completely failed using RR intervals
    for ch=1:length(ValidSignals)
        % the last full "Fm" section will be analysed till the end of the signals
        if ts_section==max([0 (floor(size(Annotation,1)/Fm)-1)])
            RRinterval = diff(find(Annotation(((ts_section*Fm)+1):end,ch)==1));
        else
            RRinterval = diff(find(Annotation(((ts_section*Fm)+1):((ts_section+1)*Fm),ch)==1));
        end
        if length(RRinterval)>9
              relRR2 = 2*(RRinterval(3:end)-RRinterval(1:end-2))./(RRinterval(3:end)+RRinterval(1:end-2));
              relRR3 = 2*(RRinterval(4:end)-RRinterval(1:end-3))./(RRinterval(4:end)+RRinterval(1:end-3));
              relRR4 = 2*(RRinterval(5:end)-RRinterval(1:end-4))./(RRinterval(5:end)+RRinterval(1:end-4));
              relRR5 = 2*(RRinterval(6:end)-RRinterval(1:end-5))./(RRinterval(6:end)+RRinterval(1:end-5));

            %accepted_RR_count(ch,ts_section+1) = sum(sum(abs([relRR2 relRR3 relRR4 relRR5])<.2,2)>0);
            counts = [sum(abs(relRR2)<.2) sum(abs(relRR3)<.2) sum(abs(relRR4)<.2) sum(abs(relRR5)<.2)];
            [accepted_RR_pct(ch,ts_section+1),tmp] = max(counts./(length(RRinterval)-[2:5]));
            accepted_RR_count(ch,ts_section+1) = counts(tmp)+tmp;
        else
            accepted_RR_count(ch,ts_section+1) = 0;
        end
    end
end


end