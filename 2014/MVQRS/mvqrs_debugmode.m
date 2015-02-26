function mvqrs_debugmode(debug_file)
%
% mvqrs_debugmode(debug_file)
%
% Matlab-compatible code for generating Figures in debugmode of mvqrs
%
% Required Parameters:
%
% debug_file
%       .mat file containing all relevant data.
%
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
% Written by Marcus Vollmer, January 20, 2015.
% Last Modified: February 26, 2015
% Version 0.2
%
%endOfHelp

% Load debugfile
load(debug_file) 
load(debug_file,'factor')

% load reference file
if exist([recordName '.atr'],'file')==2
    refAnn_dat = 'atr';
elseif exist([recordName '.ari'],'file')==2
    refAnn_dat = 'ari';
elseif exist([recordName '.ecg'],'file')==2
    refAnn_dat = 'ecg';
end

if exist('refAnn_dat','var')
       
    % get signal names
    varnames = {};
    fid = fopen([recordName '.hea']);
    while (~feof(fid))
        txtline = textscan(fid,'%s',1,'delimiter','\n');
        k = strfind(txtline{:}, ' 0 ');
        if ~isempty(k{:})
            k = k{:}; txtline = txtline{:}; txtline = txtline{:};
            currentvarname = txtline(k(end)+3:end);
            varnames = [varnames,currentvarname];
        end
    end
    fclose(fid);
    
    [refAnn,refAnn_type] = rdann(recordName,refAnn_dat);
    beat_type = {'N','L','R','B','A','a','J','S','V','r','F','e','j','n','E','/','f','Q','?'};
    refAnn = refAnn(ismember(refAnn_type,beat_type));
    refAnn_type = refAnn_type(ismember(refAnn_type,beat_type));
    % N		Normal beat (displayed as "·" by the PhysioBank ATM, LightWAVE, pschart, and psfd)
    % L		Left bundle branch block beat
    % R		Right bundle branch block beat
    % B		Bundle branch block beat (unspecified)
    % A		Atrial premature beat
    % a		Aberrated atrial premature beat
    % J		Nodal (junctional) premature beat
    % S		Supraventricular premature or ectopic beat (atrial or nodal)
    % V		Premature ventricular contraction
    % r		R-on-T premature ventricular contraction
    % F		Fusion of ventricular and normal beat
    % e		Atrial escape beat
    % j		Nodal (junctional) escape beat
    % n		Supraventricular escape beat (atrial or nodal)
    % E		Ventricular escape beat
    % /		Paced beat
    % f		Fusion of paced and normal beat
    % Q		Unclassifiable beat
    % ?		Beat not classified during learning

    refAnnRange = cell2mat(arrayfun(@colon,round(refAnn/factor)-tolerance,round(refAnn/factor)+tolerance, 'Uniform', false));
    AnnRange = cell2mat(arrayfun(@colon,Ann-tolerance,Ann+tolerance, 'Uniform', false));
    fp  = setdiff(Ann,refAnnRange(:))
    fn  = setdiff(round(refAnn/factor),AnnRange(:))


    we = tmpmax-tmpmin;
    shift_wmax = round(Fs*60/Beat_min);
    shift_wmin = round(Fs*60/Beat_max);
    
    figure
    for i=1:length(ValidSignals)
        ax((i-1)*3+1) = subplot(3,length(ValidSignals),i);
        plot(sig(:,i))
        hold on
        scatter(refAnn/factor,sig(round(refAnn/factor),i),'r','filled')
        plot(tma(:,i),'r','linewidth',2)
        try
            title(varnames{ValidSignals(i)})
        catch
            title(num2str(ValidSignals(i)))
        end
        hold off

        ax((i-1)*3+2) = subplot(3,length(ValidSignals),length(ValidSignals)+i);
        plot(sig_tmzscore(:,i))
        hold on
        scatter(refAnn/factor,sig_tmzscore(round(refAnn/factor),i),'r','filled')

        
        ax((i-1)*3+3) = subplot(3,length(ValidSignals),2*length(ValidSignals)+i);
        plot(we(:,i))
        hold on
        scatter(find(Annotation(:,i)),we(find(Annotation(:,i)),i),'MarkerEdgeColor','k','MarkerFaceColor','g','LineWidth',.5)
        h1 = scatter(Ann,we(round(Ann),i),'MarkerEdgeColor','g','LineWidth',1.5);
        h = scatter(refAnn/factor,we(round(refAnn/factor),i),'MarkerEdgeColor','r','LineWidth',1.5);
        set(get(h, 'Children'), 'Markersize', 15)
        set(get(h1, 'Children'), 'Markersize', 10)

        plot(wmin(shift_wmin:end,i),'--k')
        plot(wmax(shift_wmax:end,i),'--k')
        treshold_ts = threshold*(wmax(shift_wmax:end,i)+wmin(shift_wmin:end-shift_wmax+shift_wmin,i));
        plot(treshold_ts,'k')
        plot(valid(:,i),'r')
        plot(range(:,i),'y')
    end

    linkaxes(ax,'x');
    j=1;
    if ~isempty(fp)
        xlim([round(fp(j)/factor)-1000/factor round(fp(j)/factor)+1000/factor])
    elseif ~isempty(fn)
        xlim([round(fn(j)/factor)-1000/factor round(fn(j)/factor)+1000/factor])
    end


end