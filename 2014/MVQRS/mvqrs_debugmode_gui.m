function mvqrs_debugmode_gui(debug_file)

% Initialize and hide the GUI as it is being constructed.
f=figure('Visible','off','Position',[0,0,1200,700],'Units','normalized');%,'PaperSize',[20,13],'ResizeFcn',@setmarkersize

% Construct the components.
% Global variables
global recordName sig tma sig_tmzscore tmpmin tmpmax wmin wmax valid range DELAY thr Ann Annotation Fs factor tolerance ValidSignals Beat_min Beat_max threshold quality_RR accepted_RR_pct accepted_RR_count quality ts_section pacemaker;
global fp fn fp05 fn05 refAnn varnames c1 c2 x xl db_file db_names
c1 = [1 1 1];
c2 = [0 0 0];
x = [];
if isempty(strfind(debug_file,'_debug.mat'))
    debug_file = [debug_file '_debug.mat'];
end
    
db_file = debug_file;
fp = [];
fn = [];

load(db_file)
% sections
% ts_section*30*Fs

recordName
% load reference file
if exist([recordName '.atr'],'file')==2
    refAnn_dat = 'atr';
elseif exist([recordName '.ari'],'file')==2
    refAnn_dat = 'ari';
elseif exist([recordName '.ecg'],'file')==2
    refAnn_dat = 'ecg';
end
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
fp  = setdiff(Ann,refAnnRange(:));
fn  = setdiff(round(refAnn/factor),AnnRange(:));

refAnnRange = cell2mat(arrayfun(@colon,round(refAnn/factor)-round(tolerance/3),round(refAnn/factor)+round(tolerance/3), 'Uniform', false));
AnnRange = cell2mat(arrayfun(@colon,Ann-round(tolerance/3),Ann+round(tolerance/3), 'Uniform', false));
fp05  = setdiff(Ann,refAnnRange(:));
fn05  = setdiff(round(refAnn/factor),AnnRange(:));


pos = strfind(db_file,'/');
if isempty(pos)
    dbf = dir('*_debug.mat');    
else
    dbf = dir([db_file(1:pos(end)) '/*_debug.mat']);
end
db_names = {dbf.name}.';

%Popups   
htext_fp = uicontrol('Style','text','String',['false positive (n=' num2str(length(fp)) ')'],'Units','normalized','Position',[.65 .95 .15 .04],'FontUnits','normalized');   
htext_fn = uicontrol('Style','text','String',['false negative (n=' num2str(length(fn)) ')'],'Units','normalized','Position',[.8 .95 .15 .04],'FontUnits','normalized');   

hpopup_fp = uicontrol('Style','popup','String',fp,'Units','normalized','Position',[.65 .925 .15 .04],'FontUnits','normalized','Callback', @popup_Callback);   
hpopup_fn = uicontrol('Style','popup','String',fn,'Units','normalized','Position',[.8 .925 .15 .04],'FontUnits','normalized','Callback', @popup_Callback);        

hpopup_fp05 = uicontrol('Style','popup','String',fp05,'Units','normalized','Position',[.65 .885 .15 .04],'FontUnits','normalized','Callback', @popup_Callback);   
hpopup_fn05 = uicontrol('Style','popup','String',fn05,'Units','normalized','Position',[.8 .885 .15 .04],'FontUnits','normalized','Callback', @popup_Callback);        

if isempty(fp)
    set(htext_fp,'visible','off')
    set(hpopup_fp,'visible','off')
end
if isempty(fn)
    set(htext_fn,'visible','off')
    set(hpopup_fn,'visible','off')
end
if isempty(fp05)
    set(hpopup_fp05,'visible','off')
end
if isempty(fn05)
    set(hpopup_fn05,'visible','off')
end

htext_recordName = uicontrol('Style','text','String','record name','Units','normalized','Position',[.525 .95 .1 .04],'FontUnits','normalized');
hpopup_file = uicontrol('Style','popup','String',{dbf.name},'Value',find(ismember(db_names,db_file)),'Units','normalized','Position',[.525 .925 .1 .04],'FontUnits','normalized','Callback', @setfile_Callback);   

%Button
hbutton = uicontrol('String','switch colormode','Units','normalized','Position',[.025 .95 .08 .03],'FontUnits','normalized','Callback', @button_Callback);


%Initialisation

% Initialize the GUI.
set(f, 'Color', c2);
set(hpopup_fp,'Backgroundcolor',[.5 1 .5]);
set(hpopup_fn,'Backgroundcolor',[.5 1 .5]);
set(hpopup_file,'Backgroundcolor',[.5 1 .5]);

set(f,'Name','MVQRS debugmode')  % Assign the GUI a name to appear in the window title.
movegui(f,'center')     % Move the GUI to the center of the screen.
set(f,'toolbar','figure');


%Create a plot in the axes.
get_varnames;
show

set(f,'Visible','on');	% Make the GUI visible.



%% POPUPS
function popup_Callback(source,~)
	val = source.Value;
    s = source.String;
    x = str2num(s(val,:)); 
    xlim([round(x/Fs)-5 round(x/Fs)+5])   
end

function setfile_Callback(source,~)
	val = source.Value;
    s = source.String;
    db_file = s{val};
    mvqrs_debugmode_gui(db_file)
end

%% BUTTONS
function button_Callback(~,~)
    ctmp = c1;
    c1 = c2;
    c2 = ctmp;
    set(f, 'Color',c2);
    xl = get(gca,'xlim');
    show
end


% %% MARKERSIZE
% function setmarkersize(source,~)
%     pos = get(source,'Position');
%     h = get(get(source,'children'),'children')
%      pos(3)*.5
%     set(h,'SizeData', pos(3)*.5);
% end

%% GENERAL FUNCTIONS
function show


    we = tmpmax-tmpmin;
    shift_wmax = round(Fs*60/Beat_min);
    shift_wmin = round(Fs*60/Beat_max);
    
    if length(ValidSignals)>1
        gap = .1/(length(ValidSignals)-1);
    else
        gap=0;
    end
    
    
    length(ValidSignals)
    for i=1:length(ValidSignals)
        ax((i-1)*3+1) = subplot('Position',[.075,.075+(i-1)*(.7/length(ValidSignals)+gap),...
            .265,.7/length(ValidSignals)]);
        plot((1:size(sig,1))./Fs,sig(:,i),'color',c1)
        hold on
        scatter(refAnn/(factor*Fs),sig(round(refAnn/factor),i),'r','filled')
        plot((1:size(sig,1))./Fs,tma(:,i),'r','linewidth',2)
        try
            ylabel(gca,varnames{ValidSignals(i)},'Color',c1)
        catch
            ylabel(gca,num2str(ValidSignals(i)),'Color',c1)
        end
        set(gca,'XColor',c1,'YColor',c1,'Color',c2+(.3*c1))
        set(gca,'YTickLabel',[],'XTickLabel',[])
        if i==length(ValidSignals)
            title('downsampled series','Color',c1)
        end
        hold off

        ax((i-1)*3+2) = subplot('Position',[.35,.075+(i-1)*(.7/length(ValidSignals)+gap),...
            .265,.7/length(ValidSignals)]);
        plot((1:size(sig,1))./Fs,sig_tmzscore(:,i),'color',c1)
        hold on
        scatter(refAnn/(factor*Fs),sig_tmzscore(round(refAnn/factor),i),'r','filled')
        plot((1:size(sig,1))./Fs,tmpmin(:,i),'b')
        plot((1:size(sig,1))./Fs,tmpmax(:,i),'b')
        set(gca,'XColor',c1,'YColor',c1,'Color',c2+(.3*c1))
        set(gca,'YTickLabel',[],'XTickLabel',[])
        if i==length(ValidSignals)
            title('tmzscore','Color',c1)
        end
        
        ax((i-1)*3+3) = subplot('Position',[.65,.075+(i-1)*(.7/length(ValidSignals)+gap),...
            .3,.7/length(ValidSignals)]);
        plot((1:size(sig,1))./Fs,we(:,i),'color',[0 0 1])
        hold on
        h2 = scatter(find(Annotation(:,i))./Fs,we(find(Annotation(:,i)),i),...
            'MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',2);
        h1 = scatter(Ann./Fs,we(round(Ann),i),...
            'MarkerEdgeColor','g','LineWidth',2);
        h = scatter(refAnn/(factor*Fs),we(round(refAnn/factor),i),...
            'MarkerEdgeColor','r','LineWidth',2);
        set(h,'SizeData',2^7)
        set(h1,'SizeData',2^6)
        set(h2,'SizeData',2^5)

        plot((1:(size(sig,1)-shift_wmin+1))./Fs,wmin(shift_wmin:end,i),'--k')
        plot((1:(size(sig,1)-shift_wmax+1))./Fs,wmax(shift_wmax:end,i),'--k')
        plot((1:size(sig,1))./Fs,thr(:,i),'k')
        plot((1:size(sig,1))./Fs,valid(:,i),'r')
        plot((1:size(sig,1))./Fs,range(:,i),'r')
        set(gca,'XColor',c1,'YColor',c1,'Color',c2+(.3*c1))
        if i==length(ValidSignals)
            title('filtered','Color',c1)
        end
        if i==1
            xlabel(gca,'sec')
        else
            set(gca,'XTickLabel',[])
        end        
%         set(gca,'YTickLabel',[])
    end

    linkaxes(ax,'x');
    if ~isempty(xl)
        xlim(xl)
    end
end

function get_varnames
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
end


end

