function [output,errors] = sochro_connectivity(D,statcon,dyncon,mtdcon,surrn,mingroup,minwindow,fisherz,bgfactor,control)

%%%%%%%%%%%%%%%%%%% DATASET CHECK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errors = cell(1);

%% Sort by subject first
if ~isfield(D.parameters,'grouporder') %sort subjects by group specification
    gg = length(D.groups)+1;
    for ii = 1:length(D.group_subj) %for each unassigned subject
        inx = find(strcmp(D.group_subj{ii},{D.data.subj}));
        D.data(inx).group=gg; %assign group to OUTPUT.data
        D.groups{gg}{1}=D.group_subj{ii};
        gg=gg+1;
    end
    for gg = 1:length(D.groups), %for each group
        for ii = 1:length(D.groups{gg}), %for each subject in group gg
            inx = find(strcmp(D.groups{gg}{ii},{D.data.subj}));
            D.data(inx).group=gg; %assign group to OUTPUT.data
        end
    end
    %sort data by group
    for ii = 1:length(D.data), %copy processed EDA to OUTPUT.data
        D.data(ii).EDA = D.EDA(:,ii);
        D.data(ii).ACC = D.ACC(:,ii);
    end
    fields=fieldnames(D.data);
    cellz=struct2cell(D.data);
    sz =size(cellz);
    cellz=reshape(cellz,sz(1),[]);
    cellz=cellz';
    inx = find(strcmp('group',fields));
    cellz=sortrows(cellz,inx);
    cellz=reshape(cellz',sz);
    D.data=cell2struct(cellz,fields,1);
    D.EDA = [D.data.EDA];
    D.ACC = [D.data.ACC];
    D.parameters.grouporder = [D.data.group];
end

%% Checking conditions for minimum length
for kk = 1:length(D.conditions) %go through each condition
    D.conditions(kk).name=strtrim(D.conditions(kk).name); %get rid of spaces in condition names
    dinx=[]; %empty delete indices
    for ww = 1:length(D.conditions(kk).windows)
        if D.conditions(kk).windows{ww}{1}(1)==D.conditions(kk).windows{ww}{1}(2) %if time window is of length zero
            dinx=[dinx ww];
        end
        if (D.conditions(kk).windows{ww}{1}(2)-D.conditions(kk).windows{ww}{1}(1))<(minwindow/60)
            dinx=[dinx ww];
        end
    end
    D.conditions(kk).windows(dinx)=[]; %remove windows
    D.conditions(kk).windowsinx(dinx)=[]; %remove windows indices
end

dinx = [];
for kk = 1:length(D.conditions)
    if isempty(D.conditions(kk).windows) %Check if no windows remain
        dinx=[dinx kk];
    end
end
D.conditions(dinx)=[]; %remove empty conditions
if isempty(D.conditions)
    disp(['No conditions remain for dataset. Connectivity analysis terminated']);
    return
end
%% Load between group factor assignment
if ~isempty(bgfactor)
    [~,~,demo] = xlsread(bgfactor); %load as cell matrix
    inx = find(cellfun(@isnan,demo(:,1)));
    demo(find(cellfun(@isnan,demo(:,1))),:)=[]; %remove subject entries that are NaN
    inx=[];
    for ii = 1:length(demo(1,:))
        if sum(cellfun(@isnan,demo(:,ii)))==length(demo(:,ii))
            inx = [inx ii];
        end
    end
    demo(:,inx)=[]; %remove empty columns
    
    for ii = 1:length(demo(:,1)) %for each ID
        if isnumeric(demo{ii,1}) %if ID is numeric
            demo{ii,1}=num2str(demo{ii,1}); %convert to string for string comparison
        end
    end
    
    for ii = 1:length(D.data) %for each participant
        temp=strsplit(D.data(ii).subj,'_'); %if underscore, remove
%         D.data(ii).subj=temp(1);
        inx = find(strcmpi(temp(1),demo(:,1)));
        if isempty(inx)
            D.data(ii).bgfactor=-1;
        elseif isnan(demo{inx,2})
            D.data(ii).bgfactor=-1;
        elseif demo{inx,2}<1 || demo{inx,2}>2
            D.data(ii).bgfactor=-1;
        else D.data(ii).bgfactor=demo{inx,2};
        end
    end
end

%% Check for dead channels (no change over timewindow)
winlength=D.parameters.EDAsamplerate*10;
winstep=D.parameters.EDAsamplerate*5;
wininx = 1:winstep:length(D.EDA(:,1))-winlength; %indices for entire time series
for ww = 1:length(wininx)
    inx = find(nanstd(D.EDA(wininx(ww):wininx(ww)+winlength,:))==0); %find participants with non-varying EDA
    if ~isempty(inx) %if non-varying participant found in window ww
        dinx = find(~isnan(nanmean(D.EDA(wininx(ww):wininx(ww)+winlength,inx))));
        if ~isempty(dinx)
            D.EDA(wininx(ww):wininx(ww)+winlength,inx(dinx))=NaN;
            for dd = 1:length(dinx)
                %disp(['Datacheck: Dead channel removed in subject ' OUTPUT.data(inx(dinx(dd))).subj ' at ' num2str(OUTPUT.time.EDAmins(wininx(ww))) 'mins'])
            end
        end
    end
end

%% Check for identical channels
corrcheck = corr(D.EDA,'rows','pairwise');
for ii = 1:length(corrcheck)
    corrcheck(ii,ii)=NaN; %remove autocorrelate
end
[dinx,~]=find(corrcheck>=.99);
if ~isempty(dinx)
    disp('Following subjects are identical and will be removed:')
    for dd = 1:length(dinx)
        disp(D.data(dinx(dd)).subj);
        errors{end+1,1}=['Participant ' D.data(dinx(dd)).subj ...
            'was identical to another and removed'];
        %Remove from group assignment
        for gg = 1:length(D.groups) %go through each group to find and remove the subject
            ginx = find(strcmpi(D.data(dinx(dd)).subj,D.groups{gg}));
            if ~isempty(ginx)
                D.groups{gg}(ginx)=[];
            end
        end
        %Remove data
        D.data(dinx(dd))=[];
        D.EDA(:,dinx(dd))=[];
        D.ACC(:,dinx(dd))=[];
    end
end


%%%%%%%%%%%%%%%% SURROGATE DATASET CREATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if surrn>0
    disp('Creating null distribution');
    surr=[];
    surr(1:length(D.EDA),1:surrn,1:length(D.data))=NaN;
    for ii = 1:length(D.data)
        disp(['Participant ' num2str(ii) ' of ' num2str(length(D.data))]);
        X = D.EDA(:,ii); %get subject ii EDA
        ninx = find(isnan(X));
        X(ninx)=[];
        Xsurr=zeros(length(X),surrn);
        Xsurr=sochro_phaseran(X,surrn);
        surr(1:length(X),:,ii)=Xsurr;
    end
end

%%%%%%%%%%%% STATIC CONNECTIVITY ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if statcon==1
    disp('Running static connectivity analysis');
    if isempty(D.conditions)
        errordlg('Conditions must be defined before running connectivity analysis');
        disp(['No conditions specified for ' filename{ff}]);
        return
    end
    
    %% Get condition name
    for kk = 1:length(D.conditions)
       static(kk).condition=D.conditions(kk).name; 
    end
    
    %% Thresholding for graph metrics
    if surrn>1 %if surrogate dataset created
        Xsurr=[];
        Xsurr(1:length(D.data),1:length(D.data),1:surrn)=NaN; %preallocate
        disp('Correlating surrogate datasets for static connectivity analysis');
        for ii=1:surrn %for each surrogate set
            Xsurr(:,:,ii)=corr(squeeze(surr(:,ii,:)),'rows','pairwise');
        end
        Xsurr(isnan(Xsurr))=[]; %remove NaN
        Xsurr(Xsurr==1)=[]; %remove autocorrelates
        if fisherz==1 %Fisher's Z transformation
            Xsurr = .5*(log(1+Xsurr)-log(1-Xsurr));
            Xsurr(Xsurr==Inf)=2.6467; %==1
            Xsurr(Xsurr==-Inf)=-2.6467; %==-1
        end
        for kk = 1:length(D.conditions)
            static(kk).nullmean=nanmean(Xsurr);
            static(kk).nullSD=nanstd(Xsurr);
            static(kk).nullthresh=(1.96*static(1).nullSD)+static(1).nullmean;
            static(kk).nullmethod='surrogatedata';
        end
    else %if no surrogate dataset, compare to mean correlate
        Xsurr=corr(D.EDA,'rows','pairwise');
        Xsurr(Xsurr==1)=[];
        if fisherz==1 %Fisher's Z transformation
            Xsurr = .5*(log(1+Xsurr)-log(1-Xsurr));
            Xsurr(Xsurr==Inf)=2.6467;
            Xsurr(Xsurr==-Inf)=-2.6467;
        end
        for kk = 1:length(D.conditions)
            static(kk).nullmean=nanmean(Xsurr);
            static(kk).nullSD=nanstd(Xsurr);
            static(kk).nullthresh=(1.96*static(kk).nullSD)+static(ykk).nullmean;
            static(kk).nullmethod='Zscore';
        end
    end
    
    %% Connectivity computation
    disp('Running pairwise correlations');
    for kk = 1:length(D.conditions) %for each condition
        for ii = 1:length(D.data) %for each participant
            for jj = 1:length(D.data) %for each matching participant
                x=[]; y=[]; %empty variables for condition(kk) EDA
                if ii==jj %if self-correlation
                    static(kk).R(ii,jj)=NaN;
                    static(kk).Rt(ii,jj)=NaN;
                else
                    for ww = 1:length(D.conditions(kk).windowsinx) %for each condition window
                        win1 = D.conditions(kk).windowsinx{ww}{1}(1):D.conditions(kk).windowsinx{ww}{1}(2); %assign window indices to win1
                        x = [x; D.EDA(win1,ii)]; %unadjusted time series for ii
                        y = [y; D.EDA(win1,jj)]; %unadjusted time series for jj
                    end
                    static(kk).raw(ii,:)=x; %save unadjusted time series for ii
                    static(kk).raw(jj,:)=y; %save unadjusted time series for jj
                    static(kk).R(ii,jj)=corr(x,y,'rows','pairwise'); %correlate
                    static(kk).meanEDA(ii)=nanmean(x); %save mean EDA signal for window
                    if fisherz==1 %Fisher's Z transformation
                        X=.5*(log(1+static(kk).R(ii,jj))-log(1-static(kk).R(ii,jj)));
                        X(X==Inf)=2.6467;
                        X(X==-Inf)=-2.6467;
                        static(kk).R(ii,jj)=X;
                    end
                    if static(kk).R(ii,jj)>static(kk).nullthresh %if a significant correlation (edge)
                        static(kk).Rt(ii,jj)=static(kk).R(ii,jj);
                    else
                        static(kk).Rt(ii,jj)=NaN;
                    end
                end
            end
        end
    end
    
    %Check for conspicuous errors in calculation / missing data from a
    %condition
    if control==1
        dinx = [];
        for kk = 1:length(D.conditions)
            dinx=find(isnan(nanmean(static(kk).R)));
        end
        dinx = unique(dinx);
        if ~isempty(dinx)
            disp('Following subjects have no correlation values in one or more conditions and will be removed:')
            for dd = 1:length(dinx)
                disp(D.data(dinx(dd)).subj{1});
                errors{end+1,1}=['Participant ' D.data(dinx(dd)).subj ...
                    'had no correlation values in static connectivity and was removed'];
                %Remove from group assignment
                for gg = 1:length(D.groups) %go through each group to find and remove the subject
                    ginx = find(strcmpi(D.data(dinx(dd)).subj,D.groups{gg}));
                    if ~isempty(ginx)
                        D.groups{gg}(ginx)=[];
                    end
                end
            end
                %Remove data
                D.data(dinx)=[];
                D.EDA(:,dinx)=[];
                D.ACC(:,dinx)=[];
                for kk = 1:length(D.conditions)
                    static(kk).R(dinx,:)=[];
                    static(kk).R(:,dinx)=[];
                    static(kk).Rt(dinx,:)=[];
                    static(kk).Rt(:,dinx)=[];
                end
        end
    end
    
    %% Minimum group graphs
    groupsize = cellfun(@length, D.groups);
    groupsize(groupsize<=mingroup(1)-1)=[];
    groupsize(groupsize>=mingroup(2)-1)=[];
    if length(groupsize)>1
        for kk = 1:length(D.conditions)
            static(kk).mingroupR = static(kk).R;
            static(kk).mingroupRt = static(kk).Rt;
            dinx=[];
            for gg = 1:length(D.groups) %for each group
                inx=find([D.data.group]==gg); %find ingroup members
                if length(inx)<=mingroup(1)-1 || length(inx)>=mingroup(2)-1
                    dinx = [dinx inx];
                end
            end
            static(kk).mingroupR(dinx,:) = NaN;
            static(kk).mingroupRt(dinx,:) = NaN;
            static(kk).mingroupR(:,dinx) = NaN;
            static(kk).mingroupRt(:,dinx) = NaN;
        end
    else
        disp('Less than two groups meet group size requirements');
        for kk = 1:length(D.conditions)
            static(kk).mingroupR = static(kk).R;
            static(kk).mingroupRt = static(kk).Rt;
            static(kk).mingroupR(:,:)=NaN;
            static(kk).mingroupRt(:,:)=NaN;
        end
    end
    
    %% Calculate subject level metrics
    for kk = 1:length(D.conditions)
        for ii = 1:length(static(kk).R) %for each participant
            inx = find([D.data.group]==D.data(ii).group); %find within group participants
            iinx = setdiff(1:length(static(kk).R),inx); %find between group participants
            inx(inx==ii)=[]; %remove self index from within group
            
            %No Threshold (rows 1-6)
            static(kk).subject(1,ii)=nanmean(static(kk).R(ii,:)); %mean R for participant ii
            static(kk).subjectn(1,ii)=sum(~isnan(static(kk).R(ii,:))); %N for mean R
            static(kk).subject(2,ii)=nanmean(static(kk).R(ii,inx)); %mean within group R for participant ii
            static(kk).subjectn(2,ii)=sum(~isnan(static(kk).R(ii,inx))); %N for within group R
            static(kk).subject(3,ii)=nanmean(static(kk).R(ii,iinx)); %mean between group R for participant ii
            static(kk).subjectn(3,ii)=sum(~isnan(static(kk).R(ii,iinx))); %N for between R
            if length(inx)>=mingroup(1)-1 && length(inx)<=mingroup(2)-1 %If group size is within limits
                static(kk).subject(4,ii)=nanmean(static(kk).mingroupR(ii,:)); %mean R for participant ii
                static(kk).subjectn(4,ii)=sum(~isnan(static(kk).mingroupR(ii,:))); %N for mean R
                static(kk).subject(5,ii)=nanmean(static(kk).mingroupR(ii,inx)); %mean within group R for participant ii
                static(kk).subjectn(5,ii)=sum(~isnan(static(kk).mingroupR(ii,inx))); %N for within group R
                static(kk).subject(6,ii)=nanmean(static(kk).mingroupR(ii,iinx)); %mean between group R for participant ii
                static(kk).subjectn(6,ii)=sum(~isnan(static(kk).mingroupR(ii,iinx))); %N for between R
            else
                static(kk).subject(4:6,ii)=NaN; 
                static(kk).subjectn(4:6,ii)=NaN;
            end
            
            %With threshold (rows 7 - 12)
            static(kk).subject(7,ii)=nanmean(static(kk).Rt(ii,:)); %mean R for participant ii
            static(kk).subjectn(7,ii)=sum(~isnan(static(kk).Rt(ii,:))); %N for mean R
            static(kk).subject(8,ii)=nanmean(static(kk).Rt(ii,inx)); %mean within group R for participant ii
            static(kk).subjectn(8,ii)=sum(~isnan(static(kk).Rt(ii,inx))); %N for within R
            static(kk).subject(9,ii)=nanmean(static(kk).Rt(ii,iinx)); %mean between group R for participant ii
            static(kk).subjectn(9,ii)=sum(~isnan(static(kk).Rt(ii,iinx))); %N for between R
            if length(inx)>=mingroup(1)-1 && length(inx)<=mingroup(2)-1 %If group size is within limits
                static(kk).subject(10,ii)=nanmean(static(kk).mingroupRt(ii,:)); %mean R for participant ii
                static(kk).subjectn(10,ii)=sum(~isnan(static(kk).mingroupRt(ii,:))); %N for mean R
                static(kk).subject(11,ii)=nanmean(static(kk).mingroupRt(ii,inx)); %mean within group R for participant ii
                static(kk).subjectn(11,ii)=sum(~isnan(static(kk).mingroupRt(ii,inx))); %N for within group R
                static(kk).subject(12,ii)=nanmean(static(kk).mingroupRt(ii,iinx)); %mean between group R for participant ii
                static(kk).subjectn(12,ii)=sum(~isnan(static(kk).mingroupRt(ii,iinx))); %N for between R
            else
                static(kk).subject(10:12,ii)=NaN; 
                static(kk).subjectn(10:12,ii)=NaN;
            end
            

        end
    end
    
    %% Calculate graph level metrics
    %Define in-group / out-group edges
    for ii = 1:length(D.data)
        for jj = 1:length(D.data)
            if ii==jj
                ggmatrix(ii,jj)=0;
            elseif D.data(ii).group==D.data(jj).group
                ggmatrix(ii,jj)=1;
            else
                ggmatrix(ii,jj)=2;
            end
        end
    end
    ggmatrix=tril(ggmatrix,-1);
    for kk = 1:length(D.conditions)
        static(kk).graph(1)=nanmean(nanmean(static(kk).R)); %mean R value across all connections
        static(kk).graph(4)=nanmean(nanmean(static(kk).mingroupR)); %mean R value across all connections in mingroup sample
        static(kk).graph(7)=nanmean(nanmean(static(kk).Rt)); %mean R value for significant edges
        static(kk).graph(10)=nanmean(nanmean(static(kk).mingroupRt)); %mean R value for significant edges in mingroup sample
        for gg = 1:length(D.groups) %for each group
            inx=find([D.data.group]==gg); %find ingroup members
            iinx=setdiff(1:length(D.data),inx); %find outgroup members
            static(kk).group(1,gg)=nanmean(nanmean(static(kk).R(inx,inx))); %mean R value within group gg
            static(kk).group(2,gg)=nanmean(nanmean(static(kk).R(inx,iinx))); %mean R value with out-group members for group gg
            static(kk).group(3,gg)=nanmean(nanmean(static(kk).mingroupR(inx,inx))); %mean R value within group gg
            static(kk).group(4,gg)=nanmean(nanmean(static(kk).mingroupR(inx,iinx))); %mean R value with out-group members for group gg
            static(kk).group(5,gg)=nanmean(nanmean(static(kk).Rt(inx,inx))); %mean R value within group gg
            static(kk).group(6,gg)=nanmean(nanmean(static(kk).Rt(inx,iinx))); %mean R value with out-group members for group gg
            static(kk).group(7,gg)=nanmean(nanmean(static(kk).mingroupRt(inx,inx))); %mean R value within group gg
            static(kk).group(8,gg)=nanmean(nanmean(static(kk).mingroupRt(inx,iinx))); %mean R value with out-group members for group gg
        end
        static(kk).graph(2)=nanmean(static(kk).R(ggmatrix==1));
        static(kk).graph(3)=nanmean(static(kk).R(ggmatrix==2));
        static(kk).graph(5)=nanmean(static(kk).mingroupR(ggmatrix==1));
        static(kk).graph(6)=nanmean(static(kk).mingroupR(ggmatrix==2));
        static(kk).graph(8)=nanmean(static(kk).Rt(ggmatrix==1));
        static(kk).graph(9)=nanmean(static(kk).Rt(ggmatrix==2));
        static(kk).graph(11)=nanmean(static(kk).mingroupRt(ggmatrix==1));
        static(kk).graph(12)=nanmean(static(kk).mingroupRt(ggmatrix==2));
        static(kk).graphn(1)=length(find(~isnan(static(kk).R(ggmatrix>0))));
        static(kk).graphn(2)=length(find(~isnan(static(kk).R(ggmatrix==1))));
        static(kk).graphn(3)=length(find(~isnan(static(kk).R(ggmatrix==2))));
        static(kk).graphn(4)=length(find(~isnan(static(kk).mingroupR(ggmatrix>0))));
        static(kk).graphn(5)=length(find(~isnan(static(kk).mingroupR(ggmatrix==1))));
        static(kk).graphn(6)=length(find(~isnan(static(kk).mingroupR(ggmatrix==2))));
        static(kk).graphn(7)=length(find(~isnan(static(kk).Rt(ggmatrix>0))));
        static(kk).graphn(8)=length(find(~isnan(static(kk).Rt(ggmatrix==1))));
        static(kk).graphn(9)=length(find(~isnan(static(kk).Rt(ggmatrix==2))));
        static(kk).graphn(10)=length(find(~isnan(static(kk).mingroupRt(ggmatrix>0))));
        static(kk).graphn(11)=length(find(~isnan(static(kk).mingroupRt(ggmatrix==1))));
        static(kk).graphn(12)=length(find(~isnan(static(kk).mingroupRt(ggmatrix==2))));
    end
    
    %% Calculate between group factor metrics
    if ~isempty(bgfactor)
        %create bgfactor matrix
        for ii = 1:length(D.data)
            for jj = 1:length(D.data)
                if ii==jj,
                    bgmatrix(ii,jj)=0;
                else
                    bgmatrix(ii,jj)=D.data(ii).bgfactor+D.data(jj).bgfactor;
                end
            end
        end
       for kk=1:length(D.conditions)
           static(kk).bgfactor(1,1)=nanmean(static(kk).R(bgmatrix==2));
           static(kk).bgfactor(1,2)=nanmean(static(kk).R(bgmatrix==3));
           static(kk).bgfactor(1,3)=nanmean(static(kk).R(bgmatrix==4));
           static(kk).bgfactor(2,1)=nanmean(static(kk).mingroupR(bgmatrix==2));
           static(kk).bgfactor(2,2)=nanmean(static(kk).mingroupR(bgmatrix==3));
           static(kk).bgfactor(2,3)=nanmean(static(kk).mingroupR(bgmatrix==4));
           static(kk).bgfactor(3,1)=nanmean(static(kk).Rt(bgmatrix==2));
           static(kk).bgfactor(3,2)=nanmean(static(kk).Rt(bgmatrix==3));
           static(kk).bgfactor(3,3)=nanmean(static(kk).Rt(bgmatrix==4));
           static(kk).bgfactor(4,1)=nanmean(static(kk).mingroupRt(bgmatrix==2));
           static(kk).bgfactor(4,2)=nanmean(static(kk).mingroupRt(bgmatrix==3));
           static(kk).bgfactor(4,3)=nanmean(static(kk).mingroupRt(bgmatrix==4));
           
           static(kk).bgfactorn(1,1)=length(find(~isnan(static(kk).R(bgmatrix==2))));
           static(kk).bgfactorn(1,2)=length(find(~isnan(static(kk).R(bgmatrix==3))));
           static(kk).bgfactorn(1,3)=length(find(~isnan(static(kk).R(bgmatrix==4))));
           static(kk).bgfactorn(2,1)=length(find(~isnan(static(kk).mingroupR(bgmatrix==2))));
           static(kk).bgfactorn(2,2)=length(find(~isnan(static(kk).mingroupR(bgmatrix==3))));
           static(kk).bgfactorn(2,3)=length(find(~isnan(static(kk).mingroupR(bgmatrix==4))));
           static(kk).bgfactorn(3,1)=length(find(~isnan(static(kk).Rt(bgmatrix==2))));
           static(kk).bgfactorn(3,2)=length(find(~isnan(static(kk).Rt(bgmatrix==3))));
           static(kk).bgfactorn(3,3)=length(find(~isnan(static(kk).Rt(bgmatrix==4))));
           static(kk).bgfactorn(4,1)=length(find(~isnan(static(kk).mingroupRt(bgmatrix==2))));
           static(kk).bgfactorn(4,2)=length(find(~isnan(static(kk).mingroupRt(bgmatrix==3))));
           static(kk).bgfactorn(4,3)=length(find(~isnan(static(kk).mingroupRt(bgmatrix==4))));
           static(kk).bgfactorn=static(kk).bgfactorn./2; %divide by two for undirected edges
       end 
    end
    
    
    %% Get metadata for each condition
    for kk = 1:length(D.conditions)
        static(kk).metadata.groupassignment = [D.data.group];
        static(kk).metadata.mingroupn = length(groupsize);
        static(kk).metadata.groupsize = cellfun(@length, D.groups);
        static(kk).metadata.winlength = 0; %record timelength of condition
        for ww = 1:length(D.conditions(kk).windows)
            static(kk).metadata.winlength = static(kk).metadata.winlength+(D.conditions(kk).windows{ww}{1}(2)-D.conditions(kk).windows{ww}{1}(1));
        end
        static(kk).metadata.subjects = {D.data.subj};
        static(kk).metadata.N = length(D.data);
        static(kk).metadata.file = D.sdir;
        if ~isempty(bgfactor)
            static(kk).metadata.bgfactor = {D.data.bgfactor};
        end
    end
output.static=static;
end %end of static connectivity analysis

%%%%%%%%%%%% DYNAMIC CONNECTIVITY ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(dyncon)~=0 
    disp('Running dynamic connectivity analysis');
    if isempty(D.conditions)
        errordlg('Conditions must be defined before running connectivity analysis');
        disp(['No conditions specified for ' filename{ff}]);
        return
    end

    winlength=dyncon(1)*D.parameters.EDAsamplerate;
    winlength=ceil(winlength);
    winstep=dyncon(2)*D.parameters.EDAsamplerate;
    winstep=ceil(winstep);

    %% Create indices for dynamic timeseries per condition
    for kk = 1:length(D.conditions) %for each condition kk
        dinx=[];
        for ww = 1:length(D.conditions(kk).windows) %for each window for condition kk
            dynamic(kk).condition = D.conditions(kk).name;
            dynamic(kk).dyninx{ww}=D.conditions(kk).windowsinx{ww}{1}(1):winstep:D.conditions(kk).windowsinx{ww}{1}(2);
            dynamic(kk).dynmins{ww}=D.time.EDAmins(dynamic(kk).dyninx{ww});
            if dynamic(kk).dyninx{ww}(end)-dynamic(kk).dyninx{ww}(1)<(winlength+winstep) ||  (dynamic(kk).dynmins{ww}(end)-dynamic(kk).dynmins{ww}(1))*60<minwindow %if condition window ww is shorter than correlation window
                dinx= [dinx ww];
            end
        end
        dynamic(kk).dyninx(dinx)=[];
        dynamic(kk).dynmins(dinx)=[];
    end
    
    %Remove conditions without valid time windows for analysis
    dinx = [];
    for kk = 1:length(dynamic)
        if isempty(dynamic(kk).dyninx)
            dinx=[dinx kk];
        end
    end
    dynamic(dinx)=[];

    %% Whole session dynamics
    disp('Calculating whole session correlations');
    wininx = 1:winstep:length(D.EDA(:,1))-winlength; %indices for entire time series
    dynamicwhole.dynmins=D.time.EDAmins(wininx); %for forward facing
    for ww = 1:length(wininx)
        dynamicwhole.R(:,:,ww)=corr(D.EDA(wininx(ww):wininx(ww)+winlength,:),'rows','pairwise'); %run correlation on window ww
        if fisherz==1
            X=.5*(log(1+dynamicwhole.R(:,:,ww))-log(1-dynamicwhole.R(:,:,ww)));
            X(X==Inf)=2.6467;
            X(X==-Inf)=-2.6467;
            dynamicwhole.R(:,:,ww)=X;
        end
        for ii = 1:length(D.data) %for each participant
            dynamicwhole.R(ii,ii,ww)=NaN; %remove autocorrelate of participant in matrix
            inx = find([D.data.group]==D.data(ii).group); %find group members
            iinx = setdiff(1:length(dynamicwhole.R(:,:,ww)),inx); %define outgroup members
            inx(inx==ii)=[]; %remove self index from group member indices
            dynamicwhole.subject(ww,ii,1)=nanmean(dynamicwhole.R(:,ii,ww)); %mean R per ww
            dynamicwhole.subject(ww,ii,2)=nanmean(dynamicwhole.R(inx,ii,ww));  %within R per ww
            dynamicwhole.subject(ww,ii,3)=nanmean(dynamicwhole.R(iinx,ii,ww)); %between R per ww
        end
    end
    
    %% Minimum group graphs (whole session dynamics)
    groupsize = cellfun(@length, D.groups);
    groupsize(groupsize<=mingroup(1)-1)=[];
    groupsize(groupsize>=mingroup(2)-1)=[];
    if length(groupsize)>1
            dynamicwhole.mingroupR = dynamicwhole.R;
            dinx=[];
            for gg = 1:length(D.groups) %for each group
                inx=find([D.data.group]==gg); %find ingroup members
                if length(inx)<=mingroup(1)-1 || length(inx)>=mingroup(2)-1
                    dinx = [dinx inx];
                end
            end
            dynamicwhole.mingroupR(dinx,:,:) = NaN; %remove
            dynamicwhole.mingroupR(:,dinx,:) = NaN; %remove
            dynamicwhole.subject(:,:,4) = squeeze(nanmean(dynamicwhole.mingroupR))';
            for ii = 1:length(D.data) %for each participant
                inx = find([D.data.group]==D.data(ii).group); %find group members
                iinx = setdiff(1:length(D.data),inx); %define outgroup members
                inx(inx==ii)=[]; %remove self index from group member indices
                dynamicwhole.subject(:,ii,5) = nanmean(dynamicwhole.mingroupR(inx,ii,:));
                dynamicwhole.subject(:,ii,6) = nanmean(dynamicwhole.mingroupR(iinx,ii,:));
            end
    else
        disp('Less than two groups meet group size requirements');
            dynamicwhole.mingroupR = dynamicwhole.R;
            dynamicwhole.mingroupR(:,:,:)=NaN;
            dynamicwhole.subject(:,:,4:6)=dynamicwhole.subject(:,:,1:3);
            dynamicwhole.subject(:,:,4:6)=NaN;
    end
    disp('Whole session correlations complete!');
 
    %% Find thresholds for significance
    if surrn>1
        disp('Calculating surrogate data correlations');
        Xsurr = [];
        Xsurr(1:length(D.data),1:length(D.data),1:length(wininx),1:surrn)=NaN;
        for pp = 1:surrn %for each surrogate
            for ww = 1:length(wininx) %for each window
                Xsurr(:,:,ww,pp)=corr(squeeze(surr(wininx(ww):wininx(ww)+winlength-1,pp,:)),'rows','pairwise');
            end
        end
        disp('Surrogate data correlations complete!');
        Xsurr(isnan(Xsurr))=[];
        Xsurr(Xsurr==1)=[];
        if fisherz==1
            Xsurr=.5*(log(1+Xsurr)-log(1-Xsurr));
            Xsurr(Xsurr==Inf)=2.6467;
            Xsurr(Xsurr==-Inf)=-2.6467;
        end
        for kk = 1:length(D.conditions)
            dynamic(kk).nullmean=nanmean(Xsurr);
            dynamic(kk).nullSD=nanstd(Xsurr);
            dynamic(kk).nullthresh=(1.96*dynamic(kk).nullSD)+dynamic(kk).nullmean;
            dynamic(kk).nullmethod='surrogatedata';
        end
    else
        Xsurr=dynamicwhole.R;
        Xsurr(isnan(Xsurr))=[];
        for kk = 1:length(D.conditions)
            dynamic(kk).nullmean=nanmean(Xsurr);
            dynamic(kk).nullSD=nanstd(Xsurr);
            dynamic(kk).nullthresh=(1.96*dynamic(kk).nullSD)+dynamic(kk).nullmean;
            dynamic(kk).nullmethod='Zscore';
        end
    end

    %% Condition specific dynamics
    disp('Calculating condition specific correlations');
    for kk = 1:length(dynamic)
        wininx = [dynamic(kk).dyninx{:}]; %get window indices
        dinx = find((wininx+winlength>length(D.EDA)));
        wininx(dinx)=[]; %remove window indices that exceed data indices
        for ww = 1:length(wininx) %for each window
            dynamic(kk).R(:,:,ww)=corr(D.EDA(wininx(ww):wininx(ww)+winlength,:));
            dynamic(kk).raw(:,:,ww)=D.EDA(wininx(ww):wininx(ww)+winlength,:);
            if fisherz==1 %Fisher's Z transformation
                X=.5*(log(1+dynamic(kk).R(:,:,ww))-log(1-dynamic(kk).R(:,:,ww)));
                X(X==Inf)=2.6467;
                X(X==-Inf)=-2.6467;
                dynamic(kk).R(:,:,ww)=X;
            end
            dynamic(kk).Rt(:,:,ww)=dynamic(kk).R(:,:,ww)>dynamic(kk).nullthresh;
        end
        dynamic(kk).Rt=dynamic(kk).Rt.*dynamic(kk).R;
        dynamic(kk).Rt(dynamic(kk).Rt==0)=NaN;
        for ii = 1:length(D.data)
            dynamic(kk).R(ii,ii,:)=NaN;
            dynamic(kk).Rt(ii,ii,:)=NaN;
        end
    end
    
    %Check for conspicuous errors in calculation / missing data from a
    %condition
    if control==1
        dinx=[];
        for kk = 1:length(D.conditions)
            dinx=find(isnan(nanmean(nanmean(dynamic(kk).R,3))));
        end
        dinx = unique(dinx);
        if ~isempty(dinx)
            disp('Following subjects have no correlation values in one or more conditions and will be removed:')
            for dd = 1:length(dinx)
                disp(D.data(dinx(dd)).subj);
                errors{end+1,1}=['Participant ' D.data(dinx(dd)).subj ...
                    'had no correlation values in dynamic connectivity and was removed'];
                %Remove from group assignment
                for gg = 1:length(D.groups) %go through each group to find and remove the subject
                    ginx = find(strcmpi(D.data(dinx(dd)).subj,D.groups{gg}));
                    if ~isempty(ginx)
                        D.groups{gg}(ginx)=[];
                    end
                end
            end
            %Remove data
                D.data(dinx)=[];
                D.EDA(:,dinx)=[];
                D.ACC(:,dinx)=[];
                for kk = 1:length(D.conditions)
                    dynamic(kk).R(dinx,:,:)=[];
                    dynamic(kk).R(:,dinx,:)=[];
                    dynamic(kk).Rt(dinx,:,:)=[];
                    dynamic(kk).Rt(:,dinx,:)=[];
                end
        end
    end
    
    %% Minimum group graphs (condition specific)
    groupsize = cellfun(@length, D.groups);
    groupsize(groupsize<=mingroup(1)-1)=[];
    groupsize(groupsize>=mingroup(2)-1)=[];
    if length(groupsize)>1
        for kk = 1:length(D.conditions)
            dynamic(kk).mingroupR = dynamic(kk).R;
            dynamic(kk).mingroupRt = dynamic(kk).Rt;
            dinx=[];
            for gg = 1:length(D.groups) %for each group
                inx=find([D.data.group]==gg); %find ingroup members
                if length(inx)<=mingroup(1)-1 || length(inx)>=mingroup(2)-1
                    dinx = [dinx inx];
                end
            end
            dynamic(kk).mingroupR(dinx,:,:) = NaN;
            dynamic(kk).mingroupRt(dinx,:,:) = NaN;
            dynamic(kk).mingroupR(:,dinx,:) = NaN;
            dynamic(kk).mingroupRt(:,dinx,:) = NaN;
        end
    else
        disp('Less than two groups meet group size requirements');
        for kk = 1:length(D.conditions)
            dynamic(kk).mingroupR = dynamic(kk).R;
            dynamic(kk).mingroupRt = dynamic(kk).Rt;
            dynamic(kk).mingroupR(:,:,:)=NaN;
            dynamic(kk).mingroupRt(:,:,:)=NaN;
        end
    end
    
    
    
    %% Calculate subject level metrics
    for kk = 1:length(dynamic)
        for ii = 1:length(D.data)
            inx = find([D.data.group]==D.data(ii).group); %find within group participants
            iinx = setdiff(1:length(D.data),inx); %find between group participants
            inx(inx==ii)=[]; %remove self index from within group
            
            %No threshold (rows 1-6)
            dynamic(kk).subject(1,ii)=nanmean(nanmean(dynamic(kk).R(ii,:,:)));
            dynamic(kk).subject(2,ii)=nanmean(nanmean(dynamic(kk).R(ii,inx,:)));
            dynamic(kk).subject(3,ii)=nanmean(nanmean(dynamic(kk).R(ii,iinx,:)));
            dynamic(kk).subjectn(1,ii)=nanmean(nansum(~isnan(dynamic(kk).R(ii,:,:))));
            dynamic(kk).subjectn(2,ii)=nanmean(nansum(~isnan(dynamic(kk).R(ii,inx,:))));
            dynamic(kk).subjectn(3,ii)=nanmean(nansum(~isnan(dynamic(kk).R(ii,iinx,:))));
            if length(inx)>=mingroup(1)-1 && length(inx)<=mingroup(2)-1 %If group size is within limits
                dynamic(kk).subject(4,ii)=nanmean(nanmean(dynamic(kk).mingroupR(ii,:,:)));
                dynamic(kk).subject(5,ii)=nanmean(nanmean(dynamic(kk).mingroupR(ii,inx,:)));
                dynamic(kk).subject(6,ii)=nanmean(nanmean(dynamic(kk).mingroupR(ii,iinx,:)));
                dynamic(kk).subjectn(4,ii)=nanmean(nansum(~isnan(dynamic(kk).mingroupR(ii,:,:))));
                dynamic(kk).subjectn(5,ii)=nanmean(nansum(~isnan(dynamic(kk).mingroupR(ii,inx,:))));
                dynamic(kk).subjectn(6,ii)=nanmean(nansum(~isnan(dynamic(kk).mingroupR(ii,iinx,:))));
            else
                dynamic(kk).subject(4:6,ii)=NaN;
                dynamic(kk).subjectn(4:6,ii)=NaN;
            end
            
            %With threshold (rows 7-12)
            dynamic(kk).subject(7,ii)=nanmean(nanmean(dynamic(kk).Rt(ii,:,:)));
            dynamic(kk).subject(8,ii)=nanmean(nanmean(dynamic(kk).Rt(ii,inx,:)));
            dynamic(kk).subject(9,ii)=nanmean(nanmean(dynamic(kk).Rt(ii,iinx,:)));
            dynamic(kk).subjectn(7,ii)=nanmean(nansum(~isnan(dynamic(kk).Rt(ii,:,:))));
            dynamic(kk).subjectn(8,ii)=nanmean(nansum(~isnan(dynamic(kk).Rt(ii,inx,:))));
            dynamic(kk).subjectn(9,ii)=nanmean(nansum(~isnan(dynamic(kk).Rt(ii,iinx,:))));
            if length(inx)>=mingroup(1)-1 && length(inx)<=mingroup(2)-1 %If group size is within limits
                dynamic(kk).subject(10,ii)=nanmean(nanmean(dynamic(kk).mingroupRt(ii,:,:)));
                dynamic(kk).subject(11,ii)=nanmean(nanmean(dynamic(kk).mingroupRt(ii,inx,:)));
                dynamic(kk).subject(12,ii)=nanmean(nanmean(dynamic(kk).mingroupRt(ii,iinx,:)));
                dynamic(kk).subjectn(10,ii)=nanmean(nansum(~isnan(dynamic(kk).mingroupRt(ii,:,:))));
                dynamic(kk).subjectn(11,ii)=nanmean(nansum(~isnan(dynamic(kk).mingroupRt(ii,inx,:))));
                dynamic(kk).subjectn(12,ii)=nanmean(nansum(~isnan(dynamic(kk).mingroupRt(ii,iinx,:))));
            else
                dynamic(kk).subject(10:12,ii)=NaN;
                dynamic(kk).subjectn(10:12,ii)=NaN;
            end
        end
    end
    
