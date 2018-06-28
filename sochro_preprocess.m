function [EDA,ACC,time] = sochro_preprocess(EDA,ACC,time,samplerate,zscore,filters,smoothing,lindetrend,removezero,croptimestart,croptimeend)

%% TRIM TIME AND DATA
if croptimestart~=0 && croptimeend~=0
    time.EDAmins = time.EDAsecs./60;
    time.ACCmins = time.ACCsecs./60;
    %Get data for time windows
    differences1 = time.EDAmins - croptimestart;
    differences2 = time.EDAmins - croptimeend;
    differences3 = time.ACCmins - croptimestart;
    differences4 = time.ACCmins - croptimeend;

    smallest_difference = min(abs(differences1));
    lower_datapoint = find(abs(differences1) == smallest_difference);
    smallest_difference = min(abs(differences2));
    upper_datapoint = find(abs(differences2) == smallest_difference);
    
    time.EDAmins = time.EDAmins(lower_datapoint:upper_datapoint);
    EDA = EDA(lower_datapoint:upper_datapoint,:);
    
    smallest_difference = min(abs(differences3));
    lower_datapoint = find(abs(differences3) == smallest_difference);
    smallest_difference = min(abs(differences4));
    upper_datapoint = find(abs(differences4) == smallest_difference);
    
    time.ACCmins = time.ACCmins(lower_datapoint:upper_datapoint);
    ACC = ACC(lower_datapoint:upper_datapoint,:);
end

if removezero==1
   EDA(EDA==0)=NaN; 
end

if zscore==1
    for ii = 1:length(EDA(1,:))
        EDA(:,ii) = (EDA(:,ii) - nanmean(EDA(:,ii)))./nanstd(EDA(:,ii));
    end
end

if smoothing~=0
   for ii = 1:length(EDA(1,:)) 
      inx = find(~isnan(EDA(:,ii)));
      if isempty(inx)
         disp(['Signal for subject ' num2str(ii) ' has no signal - remove']); 
      else
      EDA(inx,ii)=smooth(EDA(inx,ii),smoothing*samplerate);
      end
   end
end

if sum(filters)~=0
    if filters(2)<0 || filters(3)>=(samplerate/2)
        errordlg(['ERROR - Butterworth filter frequencies must be greater than 0 and less than ' num2str(samplerate/2)]);
        return
    end
    if filters(2)>0
        for ii = 1:length(EDA(1,:))
            inx = EDA(:,ii)==0;
            iinx = isnan(EDA(:,ii));
            EDA(iinx,ii)=0;
            [x,y]=butter(filters(1),filters(2)/(samplerate/2),'high');
            EDA(:,ii)=filter(x,y,EDA(:,ii));
            EDA(inx,ii)=0;
            EDA(iinx,ii)=NaN;
        end
    end
    if filters(3)>0
        for ii = 1:length(EDA(1,:))
            inx = EDA(:,ii)==0;
            iinx = isnan(EDA(:,ii));
            EDA(iinx,ii)=0;
            [x,y]=butter(filters(1),filters(3)/(samplerate/2),'low');
            EDA(:,ii)=filter(x,y,EDA(:,ii));
            EDA(inx,ii)=0;
            EDA(iinx,ii)=NaN;
        end
    end
end

if lindetrend==1
   for ii = 1:length(EDA(1,:)) %for each signal
      inx = find(~isnan(EDA(:,ii)));
      if isempty(inx)
         disp(['Signal for subject ' num2str(ii) ' has no signal - remove']); 
      else
      EDA(inx,ii)=detrend(EDA(inx,ii),'linear');
      end
   end
end




