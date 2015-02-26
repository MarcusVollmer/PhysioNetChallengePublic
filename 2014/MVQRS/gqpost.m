function varargout=gqpost(varargin)
% see http://www.physionet.org/physiotools/wag/gqrs-1.htm
% 
%     'usage: gqpost -r RECORD [OPTIONS ...]'
%     'where RECORD is the name of the record to be analyzed, and OPTIONS may'
%     'include any of:'
%     ' -a ANNOTATOR read annotations from the specified ANNOTATOR (default: qrs)'
%     ' -c FILE     initialize parameters from the specified configuration FILE'
%     ' -f TIME     begin processing at specified time'
%     ' -h          print this usage summary'
%     ' -H          read multifrequency signals in high resolution mode'
%     ' -m THRESH   set interpolated event acceptance threshold to THRESH'
%     '              (default: 1)'
%     ' -o ANNOTATOR write annotations to the specified ANNOTATOR (default: gqp)'
%     ' -t TIME     stop processing at specified time'
%     'If too many true beats are rejected, decrease THRESH;  if too many false'
%     'detections are accepted, increase THRESH.'
%     'Note that the output is a complete copy of the input (with rejected events'
%     'flagged as ARFCT).  The -f and -t options only limit the interval during'
%     'which events may be rejected.'
%
% %Example
% N=5000;
% gqrs('mitdb/100',N);
% gqpost('mitdb/100','qrs','gqp');
% ann=rdann('mitdb/100','qrs',[],N);
% ann2=rdann('mitdb/100','qqp',[],N);
% [tm,sig]=rdsamp('mitdb/100',[],N);
% plot(tm,sig(:,1));hold on;grid on
% plot(tm(ann),sig(ann,1),'ro')
% plot(tm(ann),sig(ann2,1),'kd')
%
%endOfHelp

persistent javaWfdbExec
if(isempty(javaWfdbExec))
    javaWfdbExec=getWfdbClass('gqpost');
end

%Set default pararamter values
inputs={'recordName','annotatorLoad','annotatorWrite','timeStart','timeStop','summary','highResolution','file'};
annotatorLoad = 'qrs';
annotatorWrite = 'gqp';
% N=[];
% N0=1;
% signal=[]; %use application default
% threshold=[];%use application default
% outputName=[];
% highResolution=[];
for n=1:nargin
    if(~isempty(varargin{n}))
        eval([inputs{n} '=varargin{n};'])
    end
end

%N0=num2str(N0-1); %-1 is necessary because WFDB is 0 based indexed.
%wfdb_argument={'-r',recordName,'-f',['s' N0]};
wfdb_argument={'-r',recordName,'-a',annotatorLoad,'-o',annotatorWrite};

% if(~isempty(N))
%     wfdb_argument{end+1}='-t';
%     %-1 is necessary because WFDB is 0 based indexed.
%     wfdb_argument{end+1}=['s' num2str(N-1)];
% end
% if(~isempty(signal))
%     wfdb_argument{end+1}='-s';
%     %-1 is necessary because WFDB is 0 based indexed.
%     wfdb_argument{end+1}=num2str(signal-1);
% end
% if(~isempty(threshold))
%     wfdb_argument{end+1}='-m';
%     %-1 is necessary because WFDB is 0 based indexed.
%     wfdb_argument{end+1}=num2str(threshold-1);
% end
% if(~isempty(outputName))
%     wfdb_argument{end+1}='-o';
%     %-1 is necessary because WFDB is 0 based indexed.
%     wfdb_argument{end+1}=outputName;
% end
% if(~isempty(highResolution))
%     if(logical(highResolution))
%         wfdb_argument{end+1}='-H';
%     end
% end

err=javaWfdbExec.execToStringList(wfdb_argument);
if(~isempty(strfind(err.toString,['annopen: can''t'])))
    error(err)
end


