function [Ann,candidates] = mvqrs_multivariateAnn(Annotation,valid,interval)
%
% [Ann,candidates] = mvqrs_multivariateAnn(Annotation,valid,interval)
%
% Fusion method of delay corrected beat annotation. Slightly different beat
% positions of different signals will be merged to its mean annotation. 
%
% Required Parameters:
%
% Annotation
%       A NxM matrix with {0,1} values with delay corrected beat postions
%       of M different signals. 
% valid
%       A NxM matrix with {0,1} values of trusted beat decision making of M
%       different signals.
% interval
%       An integer variable which specifies the interval in which the same
%       heart beat appears in each signal.
%
%
% Written by Marcus Vollmer, 2014
% Last Modified: January 16, 2015
% Version 0.1
%
%endOfHelp

% list all beat annotations, sort and compute the number of valid/trusted
% signals
    [r,c] = find(Annotation);
    pos = sortrows([r c]);
    pos = [pos sum(valid(round(pos(:,1)),:),2)];

% compute the difference of beat positions
    d = diff(pos);
    
% difference smaller enough and from different signal origins
    multibeat = find(d(:,1)<=interval & d(:,2)~=0);
    iv = [multibeat multibeat+1];
      
    if ~isempty(iv)
        kk=1;        
        for k=2:size(iv,1)
            if iv(k,1)<=iv(kk,2)
                iv(kk,2) = iv(k,2);
            else
                kk = kk+1;                
                iv(kk,:) = iv(k,:);
            end             
        end
        iv = iv(1:kk,:); 
        
        iv_all = zeros(sum(iv(:,2)-iv(:,1)+1),1);
        tmp = cumsum(iv(:,2)-iv(:,1)+1);
        iv_all(1:tmp(1)) = iv(1,1):iv(1,2);
        for i=2:size(iv,1)
            iv_all(tmp(i-1)+1:tmp(i)) = iv(i,1):iv(i,2);
        end
    else
        iv_all = [];
    end 
    
% beats only in one channel
    pos_one = pos(setdiff(1:size(pos,1),iv_all),1);
    
% check the number of valid signals
    valid_number = sum(valid(pos_one,:),2); 
    Ann = []; 
    
% add beat positions where the number of valid signals is less or equal 50%
    Ann = [Ann; pos_one(valid_number<=.5*size(valid,2))];
    
% add beat positions which appears in more than one signal
    Ann = [Ann; round(mean(pos(iv),2))]; 
    Ann = sort(Ann); 
    
% keep other beat positions as candidates
    candidates = pos_one(valid_number>.5*size(valid,2));


