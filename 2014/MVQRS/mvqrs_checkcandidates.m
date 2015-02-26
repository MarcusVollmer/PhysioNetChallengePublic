function Ann = mvqrs_checkcandidates(testAnn,Fs,candidates)
%
% Ann = mvqrs_checkcandidates(testAnn,Fs,candidates)
%
% Check candidates using RR distances.
%
% Required Parameters:
%
% testAnn
%       A Nx1 vectors with with expected beat positions.
% Fs
%       The sampling frequency in Hz.
% candidates
%       A Mx1 vector with beat positions of candidates.
%
%
% Written by Marcus Vollmer, 2014
% Last Modified: January 30, 2015
% Version 0.1
%
%endOfHelp
add_beats = [];

RR = diff(testAnn);
RR(end+1) = NaN;

for i=1:size(candidates,1)
    candidates(i,2) = sum(testAnn<candidates(i)); 
end

% candidates at start or end will be accepted if distance to testAnn>80%
cand_start = candidates(candidates(:,2)==0,1);
if ~isempty(cand_start)    
    add_beats = cand_start(testAnn(1)-cand_start>.8*RR(1));   
end

cand_end = candidates(candidates(:,2)==size(testAnn,1),1);
if ~isempty(cand_end)    
    add_beats = [add_beats; cand_end(cand_end-testAnn(end)>.8*min(RR(end-2:end-1)))];   
end
candidates = candidates(candidates(:,2)>=1 & candidates(:,2)<size(testAnn,1),:);
    
% prove candidates with RR distances
if ~isempty(candidates)
    % u_candidates defines the neighborhood of candidates
    u_candidates = unique(candidates(:,2));
    u_candidates = [u_candidates-1 u_candidates+1];
    i=1; j=1;
    while j<size(u_candidates,1)
        while j<size(u_candidates,1) && u_candidates(i,2)>u_candidates(j+1,1)
            u_candidates(i,2) = u_candidates(j+1,2);
            j=j+1;
        end
        if j<size(u_candidates,1)
            i=i+1;
            j=j+1;
            u_candidates(i,:) = u_candidates(j,:);
        end
    end
    u_candidates = u_candidates(1:i,:); 
    
    % check every interval with candidates

    tol = round(.075*Fs);
    for i=1:size(u_candidates,1)
        % analyze the RR distances in the neighborhood 
        tmp_pre  = RR(max(1,u_candidates(i,1)-1):u_candidates(i,1));
        tmp      = RR(u_candidates(i,1)+1:u_candidates(i,2)-1);
        tmp_post = RR(u_candidates(i,2));
        if ~isempty(tmp_pre)
            tmp_prepost = [tmp_pre(end); tmp_post];
        else
            tmp_prepost = tmp_post;
        end
        
        % Analyze some gap
        if max(tmp)/nanmin([tmp_pre; tmp_post])>1.6 % candidate confirmed
            candidates_confirmed = candidates(candidates(:,2)>u_candidates(i,1) & candidates(:,2)<u_candidates(i,2),:);

            if isempty(tmp_pre)
                tmp_pre = tmp_post;
            end
            if max(tmp_pre)/min(tmp_pre)<1.2 && max(tmp_prepost)/min(tmp_prepost)<1.2 %harmonic rhythm
                for k=1:size(tmp,1)
                    % determine how many beats are missing
                    num_addbeats_max = tmp(k)/min([tmp_pre; tmp_post]);
                    num_addbeats_min = tmp(k)/max(tmp_prepost);                    
                    if num_addbeats_max>1.7 % threshold for adding beats
                        num_addbeats_iv = max([round(num_addbeats_min-.2) 2]):round(num_addbeats_max+.2); %determine how many beats will be added 
                        A_start = testAnn(u_candidates(i,1)+k);
                        A_end   = testAnn(u_candidates(i,1)+k+1); 
                        step = round((A_end-A_start)./num_addbeats_iv);
                        if size(num_addbeats_iv,2)==1 %case: number of beats definitely
                            num_addbeats = num_addbeats_iv;
                            step = step(1);
                        else %case: number of beats ambiguous
                            decisioncrit = zeros(size(num_addbeats_iv));
                            for l=1:size(num_addbeats_iv,2)
                                [tmp_min, ~] = min(abs(repmat(candidates_confirmed(:,1),1,num_addbeats_iv(l)-1)-repmat(A_start+step(l):step(l):A_start+(num_addbeats_iv(l)-1)*step(l),size(candidates_confirmed,1),1)));
                                decisioncrit(l) = sum(tmp_min-min(tmp_min)<tol)/(num_addbeats_iv(l)-1);
                            end
                            tmp_pos = find(decisioncrit==max(decisioncrit));
                            if size(tmp_pos,2)==1
                                num_addbeats = num_addbeats_iv(tmp_pos);
                                step = step(tmp_pos);
                            else
                                num_addbeats = round(tmp(k)/nanmean([tmp_pre; tmp_post]));
                                step = step(num_addbeats_iv==num_addbeats);
                            end                      
    %                        fprintf('number of beats ambiguous\n')
    %                         keyboard
                        end

                        optionally = A_start+step:step:A_start+(num_addbeats-1)*step;
                        [tmp_min, cand_pos] = min(abs(repmat(candidates_confirmed(:,1),1,size(optionally,2))-repmat(optionally,size(candidates_confirmed,1),1)),[],2);
                        cc_near_optionally = find(tmp_min-min(tmp_min)<=tol);
                        [~, cand_pos2] = min(abs(repmat(candidates_confirmed(:,1),1,size(optionally,2))-repmat(optionally,size(candidates_confirmed,1),1)),[],1); %used if more than one candidate in neighborhood of optional position
                        add_beats = [add_beats; candidates_confirmed(cand_pos2(unique(cand_pos(cc_near_optionally))),1); optionally(setdiff(1:size(optionally,2),cand_pos(cc_near_optionally)))'];  
                    end
                end
            else %Arrhythmia
                % do nothing
                %fprintf('Arrhythmia case\n')
    %             keyboard
                for k=1:size(tmp,1)
                    num_addbeats_max = tmp(k)/nanmin([tmp_pre; tmp_post]);
                    num_addbeats_min = tmp(k)/nanmax([tmp_pre; tmp_post]);
                    if num_addbeats_max>1.7 && num_addbeats_min<50 % threshold for adding beats
                        num_addbeats_iv = max([round(num_addbeats_min-.2) 2]):min([round(num_addbeats_max+.2) 3*max([round(num_addbeats_min-.2) 2])]); %determine how many beats will be added 
                        A_start = testAnn(u_candidates(i,1)+k);
                        A_end   = testAnn(u_candidates(i,1)+k+1);

                %arrhythmic filling option
                    % determine arrhythmic type and compute filling option
                    % temporarily do nothing

                %harmonic filling option
                        step = round((A_end-A_start)./num_addbeats_iv); %harmonic steps                
                        if size(num_addbeats_iv,2)==1 %case: number of beats definitely
                            num_addbeats = num_addbeats_iv;
                            step = step(1);
                            if tmp(k)/nanmean([tmp_pre; tmp_post])>1.5
                                decisioncrit = 1;
                            else
                                decisioncrit = 0;
                            end
                        else %case: number of beats ambiguous
                            decisioncrit = zeros(size(num_addbeats_iv));
                            for l=1:size(num_addbeats_iv,2)
                                [tmp_min, ~] = min(abs(repmat(candidates_confirmed(:,1),1,num_addbeats_iv(l)-1)-repmat(A_start+step(l):step(l):A_start+(num_addbeats_iv(l)-1)*step(l),size(candidates_confirmed,1),1)));
                                decisioncrit(l) = sum(tmp_min-min(tmp_min)<=tol)/(num_addbeats_iv(l)-1);
                            end
                            tmp_pos = find(decisioncrit==max(decisioncrit));
                            if size(tmp_pos,2)==1
                                num_addbeats = num_addbeats_iv(tmp_pos);
                                step = step(tmp_pos);
                            else
                                num_addbeats = num_addbeats_iv(min(tmp_pos));
                                step = step(tmp_pos);
                            end                    
                        end
                        optionally = A_start+step:step:A_start+(num_addbeats-1)*step; %harmonic filling

                %decision arrhythmic or harmonic filling
                        if max(decisioncrit)>=0.5 %.8*decisioncrit_arrhythmic
                            [tmp_min, cand_pos] = min(abs(repmat(candidates_confirmed(:,1),1,size(optionally,2))-repmat(optionally,size(candidates_confirmed,1),1)),[],2);
                            cc_near_optionally = find(tmp_min-min(tmp_min)<=tol);
                            [~, cand_pos2] = min(abs(repmat(candidates_confirmed(:,1),1,size(optionally,2))-repmat(optionally,size(candidates_confirmed,1),1)),[],1); %used if more than one candidate in neighborhood of optional position
                            add_beats = [add_beats; candidates_confirmed(cand_pos2(unique(cand_pos(cc_near_optionally))),1); optionally(setdiff(1:size(optionally,2),cand_pos(cc_near_optionally)))'];  
                            %fprintf('Arrhythmia case - harmonic filling\n')
                        else
                            %fprintf('Arrhythmia case - no filling\n')
                        end

                    end
                end    
    
    
            end

        else
% check for remature ventricular contraction (PVC), ventricular
% premature contraction (VPC), ventricular premature beat (VPB), or
% ventricular extrasystole (VES)

        end
    end  
end

Ann = sortrows([testAnn; add_beats]);

