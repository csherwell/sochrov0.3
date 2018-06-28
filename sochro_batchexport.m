function [GRAPH,SUBJECT,GROUPS] = sochro_batchexport(STATIC,DYNAMIC)
if isempty(STATIC)
    clear STATIC
end
if isempty(DYNAMIC)
    clear DYNAMIC
end

%% COPY DATA TO GRAPH LEVEL MATRIX
OUTPUT(1,1) = {'GraphMetrics'};
condN=2; %start filling in header after session variable
if isstruct(STATIC) %if static connectivity calculated
    for kk = 1:length(STATIC) %for each condition
        tempcond = []; %get condition name
        while isempty(tempcond)
            for ii=1:length(STATIC(kk).condition)
                if ~isempty(STATIC(kk).condition(ii).condition)
                    tempcond = STATIC(kk).condition(ii).condition;
                end
            end
        end
        tempcond= regexprep(tempcond,'\s+',''); %remove spaces from condition name
        OUTPUT(1,condN:condN+27)={ ...
            ['stat_' tempcond '_ISC'] , ['stat_' tempcond '_wISC'] , ['stat_' tempcond '_bISC'] , ... %mean intersubject coherence
            ['stat_' tempcond '_ISCgrp'] , ['stat_' tempcond '_wISCgrp'] , ['stat_' tempcond '_bISCgrp'] , ... %min group mean intersubject coherence
            ['stat_' tempcond '_FC'] , ['stat_' tempcond '_wFC'] , ['stat_' tempcond '_bFC'] ... %graph functional connectivity
            ['stat_' tempcond '_FCgrp'] , ['stat_' tempcond '_wFCgrp'] , ['stat_' tempcond '_bFCgrp'] ... %min group graph functional connectivity
            ['stat_' tempcond '_ND'] , ['stat_' tempcond '_wND'] , ['stat_' tempcond '_bND'] , ... %Node degree
            ['stat_' tempcond '_NDgrp'] , ['stat_' tempcond '_wNDgrp'] , ['stat_' tempcond '_bNDgrp'] , ... %min group Node Degree
            ['stat_' tempcond '_GD'] , ['stat_' tempcond '_wGD'] , ['stat_' tempcond '_bGD'] ... %Graph density
            ['stat_' tempcond '_GDgrp'] , ['stat_' tempcond '_wGDgrp'] , ['stat_' tempcond '_bGDgrp'] ... %min group Graph density
            ['stat_' tempcond '_time'] , ['stat_' tempcond '_nodes'] , ['stat_' tempcond '_Ngroups'] ...
            ['stat_' tempcond '_Nmingroups'] }; 
        if isfield(STATIC(kk).condition(1),'bgfactor')
           OUTPUT(1,condN+28:condN+51)= { ...
           ['stat_' tempcond '_ISC_bg11'] , ['stat_' tempcond '_ISC_bg12'] , ['stat_' tempcond '_ISC_bg22'] , ...
           ['stat_' tempcond '_ISCgrp_bg11'] , ['stat_' tempcond '_ISCgrp_bg12'] , ['stat_' tempcond '_ISCgrp_bg22'] , ...
           ['stat_' tempcond '_FC_bg11'] , ['stat_' tempcond '_FC_bg12'] , ['stat_' tempcond '_FC_bg22'] , ...
           ['stat_' tempcond '_FCgrp_bg11'] , ['stat_' tempcond '_FCgrp_bg12'] , ['stat_' tempcond '_FCgrp_bg22'] , ...
           ['stat_' tempcond '_ND_bg11'] , ['stat_' tempcond '_ND_bg12'] , ['stat_' tempcond '_ND_bg22'] , ...
           ['stat_' tempcond '_NDgrp_bg11'] , ['stat_' tempcond '_NDgrp_bg12'] , ['stat_' tempcond '_NDgrp_bg22'] , ...
           ['stat_' tempcond '_GD_bg11'] , ['stat_' tempcond '_GD_bg12'] , ['stat_' tempcond '_GD_bg22'] , ...
           ['stat_' tempcond '_GDgrp_bg11'] , ['stat_' tempcond '_GDgrp_bg12'] , ['stat_' tempcond '_GDgrp_bg22'] , ...     
           };
        end
        if kk==1 %if first condition
            set = 2; %Start in second row (after header)
            for ii = 1:length(STATIC(kk).condition) %per dataset
                if ~isempty(STATIC(kk).condition(ii).graph) %if data is in this row
                    OUTPUT(set,1)={STATIC(kk).condition(ii).label}; %copy session name
                    for pp = 1:12
                        OUTPUT(set,pp+condN-1)={STATIC(kk).condition(ii).graph(pp)}; % ISC & FC
                    end
                    for pp = 13:18
                       OUTPUT(set,pp+condN-1)={nanmean(STATIC(kk).condition(ii).subjectn(pp-6,:))}; %Mean node degree 
                    end
                    for pp = 19:24
                       OUTPUT(set,pp+condN-1)={STATIC(kk).condition(ii).graphn(pp-12)/STATIC(kk).condition(ii).graphn(pp-18)}; %Graph Density (thresholded edges divided by potential edges)
                    end
                    OUTPUT(set,25+condN-1)={STATIC(kk).condition(ii).metadata.winlength}; %copy condition length
                    OUTPUT(set,26+condN-1)={STATIC(kk).condition(ii).metadata.N}; %Number of students/nodes
                    OUTPUT(set,27+condN-1)={length(STATIC(kk).condition(ii).metadata.groupsize)}; %number of groups
                    OUTPUT(set,28+condN-1)={STATIC(kk).condition(ii).metadata.mingroupn}; %Number of mingroup modules
                    if isfield(STATIC(kk).condition(1),'bgfactor')
                       OUTPUT(set,29+condN-1)={STATIC(kk).condition(ii).bgfactor(1,1)};
                       OUTPUT(set,30+condN-1)={STATIC(kk).condition(ii).bgfactor(1,2)};
                       OUTPUT(set,31+condN-1)={STATIC(kk).condition(ii).bgfactor(1,3)};
                       OUTPUT(set,32+condN-1)={STATIC(kk).condition(ii).bgfactor(2,1)};
                       OUTPUT(set,33+condN-1)={STATIC(kk).condition(ii).bgfactor(2,2)};
                       OUTPUT(set,34+condN-1)={STATIC(kk).condition(ii).bgfactor(2,3)};
                       OUTPUT(set,35+condN-1)={STATIC(kk).condition(ii).bgfactor(3,1)};
                       OUTPUT(set,36+condN-1)={STATIC(kk).condition(ii).bgfactor(3,2)};
                       OUTPUT(set,37+condN-1)={STATIC(kk).condition(ii).bgfactor(3,3)};
                       OUTPUT(set,38+condN-1)={STATIC(kk).condition(ii).bgfactor(4,1)};
                       OUTPUT(set,39+condN-1)={STATIC(kk).condition(ii).bgfactor(4,2)};
                       OUTPUT(set,40+condN-1)={STATIC(kk).condition(ii).bgfactor(4,3)};
                       
                       OUTPUT(set,41+condN-1)={STATIC(kk).condition(ii).bgfactorn(3,1)};
                       OUTPUT(set,42+condN-1)={STATIC(kk).condition(ii).bgfactorn(3,2)};
                       OUTPUT(set,43+condN-1)={STATIC(kk).condition(ii).bgfactorn(3,3)};
                       OUTPUT(set,44+condN-1)={STATIC(kk).condition(ii).bgfactorn(4,1)};
                       OUTPUT(set,45+condN-1)={STATIC(kk).condition(ii).bgfactorn(4,2)};
                       OUTPUT(set,46+condN-1)={STATIC(kk).condition(ii).bgfactorn(4,3)};
                       OUTPUT(set,47+condN-1)={STATIC(kk).condition(ii).bgfactorn(3,1)/STATIC(kk).condition(ii).bgfactorn(1,1)};
                       OUTPUT(set,48+condN-1)={STATIC(kk).condition(ii).bgfactorn(3,2)/STATIC(kk).condition(ii).bgfactorn(1,2)};
                       OUTPUT(set,49+condN-1)={STATIC(kk).condition(ii).bgfactorn(3,3)/STATIC(kk).condition(ii).bgfactorn(1,3)};
                       OUTPUT(set,50+condN-1)={STATIC(kk).condition(ii).bgfactorn(4,1)/STATIC(kk).condition(ii).bgfactorn(2,1)};
                       OUTPUT(set,51+condN-1)={STATIC(kk).condition(ii).bgfactorn(4,2)/STATIC(kk).condition(ii).bgfactorn(2,2)};
                       OUTPUT(set,52+condN-1)={STATIC(kk).condition(ii).bgfactorn(4,3)/STATIC(kk).condition(ii).bgfactorn(2,3)};
                    end
                    set=set+1;
                end
            end
        else
            for ii = 1:length(STATIC(kk).condition) %per dataset
                if ~isempty(STATIC(kk).condition(ii).graph) %if data is this row
                    set=find(strcmp(STATIC(kk).condition(ii).label,OUTPUT(:,1)));
                    if isempty(set)
                        set=length(OUTPUT(:,1))+1;
                        OUTPUT(set,1)= {STATIC(kk).condition(ii).label}; %copy session name
                    end
                    for pp = 1:12
                        OUTPUT(set,pp+condN-1)={STATIC(kk).condition(ii).graph(pp)}; % ISC & FC
                    end
                    for pp = 13:18
                        OUTPUT(set,pp+condN-1)={nanmean(STATIC(kk).condition(ii).subjectn(pp-6,:))}; %Mean node degree 
                    end
                    for pp = 19:24
                       OUTPUT(set,pp+condN-1)={STATIC(kk).condition(ii).graphn(pp-12)/STATIC(kk).condition(ii).graphn(pp-18)}; %Graph Density (thresholded edges divided by potential edges)
                    end
                    OUTPUT(set,25+condN-1)={STATIC(kk).condition(ii).metadata.winlength}; %copy condition length
                    OUTPUT(set,26+condN-1)={STATIC(kk).condition(ii).metadata.N}; %Number of students/nodes
                    OUTPUT(set,27+condN-1)={length(STATIC(kk).condition(ii).metadata.groupsize)}; %number of groups
                    OUTPUT(set,28+condN-1)={STATIC(kk).condition(ii).metadata.mingroupn}; %Number of mingroup modules
                    if isfield(STATIC(kk).condition(1),'bgfactor')
                        OUTPUT(set,29+condN-1)={STATIC(kk).condition(ii).bgfactor(1,1)};
                        OUTPUT(set,30+condN-1)={STATIC(kk).condition(ii).bgfactor(1,2)};
                        OUTPUT(set,31+condN-1)={STATIC(kk).condition(ii).bgfactor(1,3)};
                        OUTPUT(set,32+condN-1)={STATIC(kk).condition(ii).bgfactor(2,1)};
                        OUTPUT(set,33+condN-1)={STATIC(kk).condition(ii).bgfactor(2,2)};
                        OUTPUT(set,34+condN-1)={STATIC(kk).condition(ii).bgfactor(2,3)};
                        OUTPUT(set,35+condN-1)={STATIC(kk).condition(ii).bgfactor(3,1)};
                        OUTPUT(set,36+condN-1)={STATIC(kk).condition(ii).bgfactor(3,2)};
                        OUTPUT(set,37+condN-1)={STATIC(kk).condition(ii).bgfactor(3,3)};
                        OUTPUT(set,38+condN-1)={STATIC(kk).condition(ii).bgfactor(4,1)};
                        OUTPUT(set,39+condN-1)={STATIC(kk).condition(ii).bgfactor(4,2)};
                        OUTPUT(set,40+condN-1)={STATIC(kk).condition(ii).bgfactor(4,3)};
                        
                        OUTPUT(set,41+condN-1)={STATIC(kk).condition(ii).bgfactorn(3,1)};
                        OUTPUT(set,42+condN-1)={STATIC(kk).condition(ii).bgfactorn(3,2)};
                        OUTPUT(set,43+condN-1)={STATIC(kk).condition(ii).bgfactorn(3,3)};
                        OUTPUT(set,44+condN-1)={STATIC(kk).condition(ii).bgfactorn(4,1)};
                        OUTPUT(set,45+condN-1)={STATIC(kk).condition(ii).bgfactorn(4,2)};
                        OUTPUT(set,46+condN-1)={STATIC(kk).condition(ii).bgfactorn(4,3)};
                        OUTPUT(set,47+condN-1)={STATIC(kk).condition(ii).bgfactorn(3,1)/STATIC(kk).condition(ii).bgfactorn(1,1)};
                        OUTPUT(set,48+condN-1)={STATIC(kk).condition(ii).bgfactorn(3,2)/STATIC(kk).condition(ii).bgfactorn(1,2)};
                        OUTPUT(set,49+condN-1)={STATIC(kk).condition(ii).bgfactorn(3,3)/STATIC(kk).condition(ii).bgfactorn(1,3)};
                        OUTPUT(set,50+condN-1)={STATIC(kk).condition(ii).bgfactorn(4,1)/STATIC(kk).condition(ii).bgfactorn(2,1)};
                        OUTPUT(set,51+condN-1)={STATIC(kk).condition(ii).bgfactorn(4,2)/STATIC(kk).condition(ii).bgfactorn(2,2)};
                        OUTPUT(set,52+condN-1)={STATIC(kk).condition(ii).bgfactorn(4,3)/STATIC(kk).condition(ii).bgfactorn(2,3)};
                    end
                end
            end
        end
        if isfield(STATIC(kk).condition(1),'bgfactor')
            condN=condN+52;
        else
            condN=condN+28;
        end
    end %end of per condition
end

for i=1:2
% % % if isstruct(DYNAMIC) %if static connectivity calculated
% % %     for kk = 1:length(DYNAMIC) %for each condition
% % %         tempcond = []; %get condition name
% % %         while isempty(tempcond)
% % %             for ii=1:length(DYNAMIC(kk).condition)
% % %                 if ~isempty(DYNAMIC(kk).condition(ii).condition)
% % %                     tempcond = DYNAMIC(kk).condition(ii).condition;
% % %                 end
% % %             end
% % %         end
% % %         tempcond= regexprep(tempcond,'\s+',''); %remove spaces from condition name
% % %         OUTPUT(1,condN:condN+27)={['dyn_' tempcond '_FC'] , ['dyn_' tempcond '_wFC'] , ['dyn_' tempcond '_bFC'] , ... %raw FC
% % %             ['dyn_' tempcond '_FCgrp'] , ['dyn_' tempcond '_wFCgrp'] , ['dyn_' tempcond '_bFCgrp'] , ... %min group raw FC
% % %             ['dynT_' tempcond '_FC'] , ['dynT_' tempcond '_wFC'] , ['dynT_' tempcond '_bFC'] ... %thresholded FC
% % %             ['dynT_' tempcond '_FCgrp'] , ['dynT_' tempcond '_wFCgrp'] , ['dynT_' tempcond '_bFCgrp'] ... %min group thresholded FC
% % %             ['dyn_' tempcond '_ND'] , ['dyn_' tempcond '_wND'] , ['dyn_' tempcond '_bND'] , ... %Node degree
% % %             ['dyn_' tempcond '_NDgrp'] , ['dyn_' tempcond '_wNDgrp'] , ['dyn_' tempcond '_bNDgrp'] , ... %min group Node Degree
% % %             ['dyn_' tempcond '_DS'] , ['dyn_' tempcond '_wDS'] , ['dyn_' tempcond '_bDS'] ... %Node density
% % %             ['dyn_' tempcond '_DSgrp'] , ['dyn_' tempcond '_wDSgrp'] , ['dyn_' tempcond '_bDSgrp'] ... %min group Node density 
% % %             ['dyn_' tempcond '_time'] , ['dyn_' tempcond '_nodes'] , ['dyn_' tempcond '_Ngroups'] ...
% % %             ['dyn_' tempcond '_Nmingroups'] }; 
% % %         
% % %         if ~isstruct(STATIC) && kk==1 %if the first entries to output
% % %             set = 2; %Start in second row (after header)
% % %             for ii = 1:length(DYNAMIC(kk).condition) %per dataset
% % %                 if ~isempty(DYNAMIC(kk).condition(ii).graph) %if data is this row
% % %                     OUTPUT(set,1)={DYNAMIC(kk).condition(ii).label}; %copy session name
% % %                     for pp = 1:24
% % %                         OUTPUT(set,pp+condN-1)={DYNAMIC(kk).condition(ii).graph(pp)};
% % %                     end
% % %                     winlength = 0;
% % %                     for ww = 1:length(DYNAMIC(kk).condition(ii).dynmins)
% % %                         winlength = winlength + (DYNAMIC(kk).condition(ii).dynmins{ww}(end)-DYNAMIC(kk).condition(ii).dynmins{ww}(1)); 
% % %                     end
% % %                     OUTPUT(set,25+condN-1)={winlength}; %copy condition length
% % %                     OUTPUT(set,26+condN-1)={DYNAMIC(kk).condition(ii).N(1)}; %Number of students/nodes
% % %                     OUTPUT(set,27+condN-1)={length(DYNAMIC(kk).condition(ii).groupsize)}; %number of groups
% % %                     OUTPUT(set,28+condN-1)={DYNAMIC(kk).condition(ii).mingroupn}; %Number of mingroup modules
% % %                     set=set+1;
% % %                 end
% % %             end
% % %         else
% % %             for ii = 1:length(DYNAMIC(kk).condition) %per dataset
% % %                 if ~isempty(DYNAMIC(kk).condition(ii).graph) %if data is this row
% % %                     set=find(strcmp(DYNAMIC(kk).condition(ii).label,OUTPUT(:,1)));
% % %                     if isempty(set)
% % %                         set=length(OUTPUT(:,1))+1;
% % %                         OUTPUT(set,1)= {DYNAMIC(kk).condition(ii).label}; %copy session name
% % %                     end
% % %                     for pp = 1:24
% % %                         OUTPUT(set,pp+condN-1)={DYNAMIC(kk).condition(ii).graph(pp)};
% % %                     end
% % %                     winlength = 0;
% % %                     for ww = 1:length(DYNAMIC(kk).condition(ii).dynmins)
% % %                         winlength = winlength + (DYNAMIC(kk).condition(ii).dynmins{ww}(end)-DYNAMIC(kk).condition(ii).dynmins{ww}(1)); 
% % %                     end
% % %                     OUTPUT(set,25+condN-1)={winlength}; %copy condition length
% % %                     OUTPUT(set,26+condN-1)={DYNAMIC(kk).condition(ii).N(1)}; %Number of students/nodes
% % %                     OUTPUT(set,27+condN-1)={length(DYNAMIC(kk).condition(ii).groupsize)}; %number of groups
% % %                     OUTPUT(set,28+condN-1)={DYNAMIC(kk).condition(ii).mingroupn}; %Number of mingroup modules
% % %                 end
% % %             end
% % %         end
% % %         condN=condN+28;
% % %     end %end of per condition
% % % end
end
GRAPH=OUTPUT;
clear OUTPUT

%% COPY DATA TO GROUP LEVEL MATRIX
OUTPUT(1,1) = {'Session'};
OUTPUT(1,2) = {'Group'};
condN=3; %start filling in header after session variable
if isstruct(STATIC) %if static connectivity calculated
    for kk = 1:length(STATIC) %for each condition
        tempcond = []; %get condition name
        while isempty(tempcond)
            for ii=1:length(STATIC(kk).condition)
                if ~isempty(STATIC(kk).condition(ii).condition)
                    tempcond = STATIC(kk).condition(ii).condition;
                end
            end
        end
        tempcond= regexprep(tempcond,'\s+',''); %remove spaces from condition name
        OUTPUT(1,condN:condN+10)={ ...
            ['stat_' tempcond '_wISC'] , ['stat_' tempcond '_bISC'] , ... %mean intersubject coherence
            ['stat_' tempcond '_wISCgrp'] , ['stat_' tempcond '_bISCgrp'] , ... %min group mean intersubject coherence
            ['stat_' tempcond '_wFC'] , ['stat_' tempcond '_bFC'] ... %graph functional connectivity
            ['stat_' tempcond '_wFCgrp'] , ['stat_' tempcond '_bFCgrp'] ... %min group graph functional connectivity
            ['stat_' tempcond '_wnodes'] , ['stat_' tempcond '_bnodes'] ...
            ['stat_' tempcond '_time']
            };
        if kk==1 %if first condition
            set = 2; %start in second row (after header)
            for ii = 1:length(STATIC(kk).condition) %per dataset
                if ~isempty(STATIC(kk).condition(ii).graph) %if data is in this row
                    for jj = 1:length(STATIC(kk).condition(ii).metadata.groupsize) %for each group
                        OUTPUT(set,1)={STATIC(kk).condition(ii).label}; %copy session name
                        OUTPUT(set,2)={num2str(jj)}; %copy session name
                        for pp = 1:8
                           OUTPUT(set,pp+condN-1)={STATIC(kk).condition(ii).group(pp,jj)}; %ISC & FC 
                        end
                        OUTPUT(set,9+condN-1)={STATIC(kk).condition(ii).metadata.groupsize(jj)};
                        OUTPUT(set,10+condN-1)={STATIC(kk).condition(ii).metadata.N-STATIC(kk).condition(ii).metadata.groupsize(jj)};
                        OUTPUT(set,11+condN-1)={STATIC(kk).condition(ii).metadata.winlength};
                        set=set+1;
                    end
                end
            end
        else
            for ii = 1:length(STATIC(kk).condition) %per dataset
                if ~isempty(STATIC(kk).condition(ii).graph) %if data is in this row
                    for jj = 1:length(STATIC(kk).condition(ii).metadata.groupsize)
                        set=find(strcmpi(STATIC(kk).condition(ii).label,OUTPUT(:,1)));
                        sesh=find(strcmp(num2str(jj),OUTPUT(:,2)));
                        set=intersect(set,sesh);
                        if length(set)>1
                            warndlg(['Multiple matching labels: Group - ' num2str(jj) ' ' STATIC(kk).condition(ii).label], 'Critical error');
                            return
                        elseif isempty(set)
                            set=length(OUTPUT(:,1))+1;
                        end
                        OUTPUT(set,1)={STATIC(kk).condition(ii).label}; %copy session name
                        OUTPUT(set,2)={num2str(jj)}; %copy session name
                        for pp = 1:8
                           OUTPUT(set,pp+condN-1)={STATIC(kk).condition(ii).group(pp,jj)}; %ISC & FC 
                        end
                        OUTPUT(set,9+condN-1)={STATIC(kk).condition(ii).metadata.groupsize(jj)};
                        OUTPUT(set,10+condN-1)={STATIC(kk).condition(ii).metadata.N-STATIC(kk).condition(ii).metadata.groupsize(jj)};
                        OUTPUT(set,11+condN-1)={STATIC(kk).condition(ii).metadata.winlength};
                    end
                end
            end
        end
        condN=condN+11;
    end
end
GROUPS=OUTPUT;
clear OUTPUT

%% COPY DATA TO SUBJECT LEVEL MATRIX
OUTPUT(1,1) = {'Subject'};
OUTPUT(1,2) = {'Session'};
condN=3; %start filling in header after session variable
if isstruct(STATIC) %if static connectivity calculated
    for kk = 1:length(STATIC) %for each condition
        tempcond = []; %get condition name
        while isempty(tempcond)
            for ii=1:length(STATIC(kk).condition)
                if ~isempty(STATIC(kk).condition(ii).condition)
                    tempcond = STATIC(kk).condition(ii).condition;
                end
            end
        end
        tempcond= regexprep(tempcond,'\s+',''); %remove spaces from condition name
        OUTPUT(1,condN:condN+28)={ ...
            ['stat_' tempcond '_ISC'] , ['stat_' tempcond '_wISC'] , ['stat_' tempcond '_bISC'] , ... %mean intersubject coherence
            ['stat_' tempcond '_ISCgrp'] , ['stat_' tempcond '_wISCgrp'] , ['stat_' tempcond '_bISCgrp'] , ... %min group mean intersubject coherence
            ['stat_' tempcond '_FC'] , ['stat_' tempcond '_wFC'] , ['stat_' tempcond '_bFC'] ... %functional connectivity
            ['stat_' tempcond '_FCgrp'] , ['stat_' tempcond '_wFCgrp'] , ['stat_' tempcond '_bFCgrp'] ... %min group functional connectivity
            ['stat_' tempcond '_ND'] , ['stat_' tempcond '_wND'] , ['stat_' tempcond '_bND'] , ... %Node degree
            ['stat_' tempcond '_NDgrp'] , ['stat_' tempcond '_wNDgrp'] , ['stat_' tempcond '_bNDgrp'] , ... %min group Node Degree
            ['stat_' tempcond '_DS'] , ['stat_' tempcond '_wDS'] , ['stat_' tempcond '_bDS'] ... %node density
            ['stat_' tempcond '_DSgrp'] , ['stat_' tempcond '_wDSgrp'] , ['stat_' tempcond '_bDSgrp'] ... %min group node density
            ['stat_' tempcond '_time'] , ['stat_' tempcond '_wN'] , ['stat_' tempcond '_bN'] ...
            ['stat_' tempcond '_wNgrp'] , ['stat_' tempcond '_bNgrp']
            };
        if isfield(STATIC(kk).condition(1),'bgfactor')
            OUTPUT(1,condN+29) = {'bgfactor'};
        end
        if kk==1 %if first condition
            set = 2; %Start in second row (after header)
            for ii = 1:length(STATIC(kk).condition) %per dataset
                if ~isempty(STATIC(kk).condition(ii).subject) %if data is in this row
                    for jj = 1:length(STATIC(kk).condition(ii).subject(1,:)) %per subject
                        OUTPUT(set,1)={STATIC(kk).condition(ii).metadata.subjects{jj}}; %copy subject ID
                        OUTPUT(set,2)={STATIC(kk).condition(ii).label}; %copy session name
                        for pp = 1:12
                            OUTPUT(set,pp+condN-1)={STATIC(kk).condition(ii).subject(pp,jj)}; %ISC & FC
                        end
                        for pp = 13:18
                            OUTPUT(set,pp+condN-1)={STATIC(kk).condition(ii).subjectn(pp-6,jj)}; %Node Degree
                        end
                        for pp = 19:24
                            OUTPUT(set,pp+condN-1)={STATIC(kk).condition(ii).subjectn(pp-12,jj)/STATIC(kk).condition(ii).subjectn(pp-18,jj)}; %Node Density (thresholded edges divided by potential edges)
                        end
                        OUTPUT(set,25+condN-1)={STATIC(kk).condition(ii).metadata.winlength};
                        OUTPUT(set,26+condN-1)={STATIC(kk).condition(ii).subjectn(2,jj)};
                        OUTPUT(set,27+condN-1)={STATIC(kk).condition(ii).subjectn(3,jj)};
                        OUTPUT(set,28+condN-1)={STATIC(kk).condition(ii).subjectn(5,jj)};
                        OUTPUT(set,29+condN-1)={STATIC(kk).condition(ii).subjectn(6,jj)};
                        if isfield(STATIC(kk).condition(1),'bgfactor')
                        OUTPUT(set,30+condN-1)=(STATIC(kk).condition(ii).metadata.bgfactor(jj));
                        end
                        set=set+1; %increase row index
                    end
                end
            end
        else
            for ii = 1:length(STATIC(kk).condition) %per dataset
                if ~isempty(STATIC(kk).condition(ii).graph) %if data is in this row
                    for jj = 1:length(STATIC(kk).condition(ii).subject(1,:)) %per subject
                        set=find(strcmpi(STATIC(kk).condition(ii).metadata.subjects(jj),OUTPUT(:,1)));
                        sesh=find(strcmp(STATIC(kk).condition(ii).label,OUTPUT(:,2)));
                        set=intersect(set,sesh);
                        if length(set)>1
                            warndlg(['Multiple matching labels: ' STATIC(kk).condition(ii).subjects(jj) ' ' STATIC(kk).condition(ii).label], 'Critical error');
                            return
                        elseif isempty(set)
                            set=length(OUTPUT(:,1))+1;
                        end
                        OUTPUT(set,1)={STATIC(kk).condition(ii).metadata.subjects{jj}}; %copy subject ID
                        OUTPUT(set,2)={STATIC(kk).condition(ii).label}; %copy session name
                        for pp = 1:12
                            OUTPUT(set,pp+condN-1)={STATIC(kk).condition(ii).subject(pp,jj)}; %ISC & FC
                        end
                        for pp = 13:18
                            OUTPUT(set,pp+condN-1)={STATIC(kk).condition(ii).subjectn(pp-6,jj)}; %Node Degree
                        end
                        for pp = 19:24
                            OUTPUT(set,pp+condN-1)={STATIC(kk).condition(ii).subjectn(pp-12,jj)/STATIC(kk).condition(ii).subjectn(pp-18,jj)}; %Node Density (thresholded edges divided by potential edges)
                        end
                        OUTPUT(set,25+condN-1)={STATIC(kk).condition(ii).metadata.winlength};
                        OUTPUT(set,26+condN-1)={STATIC(kk).condition(ii).subjectn(2,jj)};
                        OUTPUT(set,27+condN-1)={STATIC(kk).condition(ii).subjectn(3,jj)};
                        OUTPUT(set,28+condN-1)={STATIC(kk).condition(ii).subjectn(5,jj)};
                        OUTPUT(set,29+condN-1)={STATIC(kk).condition(ii).subjectn(6,jj)};
                        if isfield(STATIC(kk).condition(1),'bgfactor')
                        OUTPUT(set,30+condN-1)=STATIC(kk).condition(ii).metadata.bgfactor(jj);
                        end
                    end
                end
            end
        end
        if isfield(STATIC(kk).condition(1),'bgfactor')
            condN=condN+30;
        else
            condN=condN+29;
        end
    end %end of per condition
end

for i =1:2
% if isstruct(DYNAMIC) %if dynamic connectivity calculated
%     for kk = 1:length(DYNAMIC) %for each condition
%         tempcond = []; %get condition name
%         while isempty(tempcond)
%             for ii=1:length(DYNAMIC(kk).condition)
%                 if ~isempty(DYNAMIC(kk).condition(ii).condition)
%                     tempcond = DYNAMIC(kk).condition(ii).condition;
%                 end
%             end
%         end
%         tempcond= regexprep(tempcond,'\s+',''); %remove spaces from condition name
%         OUTPUT(1,condN:condN+24)={['dyn_' tempcond '_FC'] , ['dyn_' tempcond '_wFC'] , ['dyn_' tempcond '_bFC'] , ... %raw FC
%             ['dyn_' tempcond '_FCgrp'] , ['dyn_' tempcond '_wFCgrp'] , ['dyn_' tempcond '_bFCgrp'] , ... %min group raw FC
%             ['dynT_' tempcond '_FC'] , ['dynT_' tempcond '_wFC'] , ['dynT_' tempcond '_bFC'] ... %thresholded FC
%             ['dynT_' tempcond '_FCgrp'] , ['dynT_' tempcond '_wFCgrp'] , ['dynT_' tempcond '_bFCgrp'] ... %min group thresholded FC
%             ['dyn_' tempcond '_ND'] , ['dyn_' tempcond '_wND'] , ['dyn_' tempcond '_bND'] , ... %Node degree
%             ['dyn_' tempcond '_NDgrp'] , ['dyn_' tempcond '_wNDgrp'] , ['dyn_' tempcond '_bNDgrp'] , ... %min group Node Degree
%             ['dyn_' tempcond '_DS'] , ['dyn_' tempcond '_wDS'] , ['dyn_' tempcond '_bDS'] ... %Node density
%             ['dyn_' tempcond '_DSgrp'] , ['dyn_' tempcond '_wDSgrp'] , ['dyn_' tempcond '_bDSgrp'] ...
%             ['dyn_' tempcond '_time'] };
%             %min group Node density
%         if ~isstruct(STATIC) && kk==1 %if first condition
%             set = 2; %Start in second row (after header)
%             for ii = 1:length(DYNAMIC(kk).condition) %per dataset
%                 if ~isempty(DYNAMIC(kk).condition(ii).graph) %if data is in this row
%                     for jj = 1:length(DYNAMIC(kk).condition(ii).subject(1,:)) %per subject
%                         OUTPUT(set,1)={DYNAMIC(kk).condition(ii).subjects{jj}}; %copy subject ID
%                         OUTPUT(set,2)={DYNAMIC(kk).condition(ii).label}; %copy session name
%                         for pp = 1:12
%                             OUTPUT(set,pp+condN-1)={DYNAMIC(kk).condition(ii).subject(pp,jj)};
%                         end
%                         for pp = 1:12
%                             OUTPUT(set,pp+12+condN-1)={DYNAMIC(kk).condition(ii).nodedegree(pp,jj)};
%                         end
%                         winlength = 0;
%                         for ww = 1:length(DYNAMIC(kk).condition(ii).dynmins)
%                             winlength = winlength + (DYNAMIC(kk).condition(ii).dynmins{ww}(end)-DYNAMIC(kk).condition(ii).dynmins{ww}(1));
%                         end
%                         OUTPUT(set,25+condN-1)={winlength}; %copy condition length
%                         set=set+1; %increase row index
%                     end
%                 end
%             end
%         else
%             for ii = 1:length(DYNAMIC(kk).condition) %per dataset
%                 if ~isempty(DYNAMIC(kk).condition(ii).graph) %if data is in this row
%                     for jj = 1:length(DYNAMIC(kk).condition(ii).subject(1,:)) %per subject
%                         set=find(strcmpi(DYNAMIC(kk).condition(ii).subjects(jj),OUTPUT(:,1)));
%                         sesh=find(strcmp(DYNAMIC(kk).condition(ii).label,OUTPUT(:,2)));
%                         set=intersect(set,sesh);
%                         if length(set)>1
%                             warndlg(['Multiple matching labels: ' DYNAMIC(kk).condition(ii).subjects(jj) ' ' DYNAMIC(kk).condition(ii).label], 'Critical error');
%                             return
%                         elseif isempty(set)
%                             set=length(OUTPUT(:,1))+1;
%                         end
%                         OUTPUT(set,1)={DYNAMIC(kk).condition(ii).subjects{jj}}; %copy subject ID
%                         OUTPUT(set,2)={DYNAMIC(kk).condition(ii).label}; %copy session name
%                         for pp = 1:12
%                             OUTPUT(set,pp+condN-1)={DYNAMIC(kk).condition(ii).subject(pp,jj)};
%                         end
%                         for pp = 1:12
%                             OUTPUT(set,pp+12+condN-1)={DYNAMIC(kk).condition(ii).nodedegree(pp,jj)};
%                         end
%                         winlength = 0;
%                         for ww = 1:length(DYNAMIC(kk).condition(ii).dynmins)
%                             winlength = winlength + (DYNAMIC(kk).condition(ii).dynmins{ww}(end)-DYNAMIC(kk).condition(ii).dynmins{ww}(1));
%                         end
%                         OUTPUT(set,25+condN-1)={winlength}; %copy condition length
%                     end
%                 end
%             end
%         end
%         condN=condN+25;
%     end %end of per condition
% end
end
SUBJECT=OUTPUT;

clear OUTPUT


% %% CALCULATE METADATA
% OUTPUT(1,1) = {'Session'};
% condN=2; %start filling in header after session variable
% if isstruct(STATIC) %if static connectivity calculated
%     for kk = 1:length(STATIC) %for each condition
%         tempcond = []; %get condition name
%         while isempty(tempcond)
%             for ii=1:length(STATIC(kk).condition)
%                 if ~isempty(STATIC(kk).condition(ii).condition)
%                     tempcond = STATIC(kk).condition(ii).condition;
%                 end
%             end
%         end
%         tempcond= regexprep(tempcond,'\s+',''); %remove spaces from condition name
%         OUTPUT(1,condN:condN+3)={['stat_' tempcond '_time'] , ['stat_' tempcond '_nodes'] , ['stat_' tempcond '_Ngroups'] , ...
%             ['stat_' tempcond '_Nmingroups'] };
%         
%         if kk==1 %if first condition
%             set = 2; %Start in second row (after header)
%             for ii = 1:length(STATIC(kk).condition) %per dataset
%                 if ~isempty(STATIC(kk).condition(ii).graph) %if data is in this row
%                     OUTPUT(set,1)={STATIC(kk).condition(ii).label}; %copy session name
%                     OUTPUT(set,1+condN-1) = {STATIC(kk).condition(ii).winlength}; %copy condition length
%                     OUTPUT(set,2+condN-1) = {STATIC(kk).condition(ii).N(1)}; %Number of students/nodes
%                     OUTPUT(set,3+condN-1) = {length(STATIC(kk).condition(ii).groupsize)}; %number of groups
%                     OUTPUT(set,4+condN-1) = {STATIC(kk).condition(ii).mingroupn}; %Number of mingroup modules
%                     
%                     set=set+1;
%                 end
%             end
%         else
%             for ii = 1:length(STATIC(kk).condition) %per dataset
%                 if ~isempty(STATIC(kk).condition(ii).graph) %if data is this row
%                     set=find(strcmp(STATIC(kk).condition(ii).label,OUTPUT(:,1)));
%                     if isempty(set)
%                         set=length(OUTPUT(:,1))+1;
%                         OUTPUT(set,1)= {STATIC(kk).condition(ii).label}; %copy session name
%                     end
%                     OUTPUT(set,1+condN-1) = {STATIC(kk).condition(ii).winlength}; %copy condition length
%                     OUTPUT(set,2+condN-1) = {STATIC(kk).condition(ii).N(1)}; %Number of students/nodes
%                     OUTPUT(set,3+condN-1) = {length(STATIC(kk).condition(ii).groupsize)}; %number of groups
%                     OUTPUT(set,4+condN-1) = {STATIC(kk).condition(ii).mingroupn}; %Number of mingroup modules
%                 end
%             end
%         end
%         condN=condN+4;
%     end %end of per condition
% end
% 
% 
% METADATA=OUTPUT;

%% PROMPTS TO DELETE MISSING CASES
ButtonName = questdlg('Would you like to remove cases (sessions) that are missing data from one or more conditions?',...
    'GRAPH LEVEL MEASURES',...
    'Yes', 'No', 'Yes');

if strcmpi(ButtonName,'Yes')
    inx = [];
    for ii = 1:length(GRAPH(:,1))
        if sum(cellfun(@isempty,GRAPH(ii,:)))~=0
            inx = [inx ii];
        end
    end
    GRAPH(inx,:)=[];
end
for ii = 2:length(GRAPH(:,1))
   for jj = 2:length(GRAPH(1,:))
       if ~isempty(GRAPH{ii,jj})
          if ~isnan(GRAPH{ii,jj})
             GRAPH{ii,jj}= (round(GRAPH{ii,jj}*100000))/100000; 
          end
       end
   end
end

ButtonName = questdlg('Would you like to remove cases (subjects) that are missing data from one or more conditions?',...
    'SUBJECT LEVEL MEASURES',...
    'Yes', 'No', 'Yes');

if strcmpi(ButtonName,'Yes')
    inx = [];
    for ii = 1:length(SUBJECT(:,1))
        if sum(cellfun(@isempty,SUBJECT(ii,:)))~=0
            inx = [inx ii];
        end
    end
    SUBJECT(inx,:)=[];
end
for ii = 2:length(SUBJECT(:,1))
   for jj = 3:length(SUBJECT(1,:))
       if ~isempty(SUBJECT{ii,jj})
          if ~isnan(SUBJECT{ii,jj})
             SUBJECT{ii,jj}= (round(SUBJECT{ii,jj}*100000))/100000; 
          end
       end
   end
end

ButtonName = questdlg('Would you like to remove cases (groups) that are missing data from one or more conditions?',...
    'GROUP LEVEL MEASURES',...
    'Yes', 'No', 'Yes');

if strcmpi(ButtonName,'Yes')
    inx = [];
    for ii = 1:length(GROUPS(:,1))
        if sum(cellfun(@isempty,GROUPS(ii,:)))~=0
            inx = [inx ii];
        end
    end
    GROUPS(inx,:)=[];
end
for ii = 2:length(GROUPS(:,1))
   for jj = 3:length(GROUPS(1,:))
       if ~isempty(GROUPS{ii,jj})
          if ~isnan(GROUPS{ii,jj})
             GROUPS{ii,jj}= (round(GROUPS{ii,jj}*100000))/100000; 
          end
       end
   end
end

[filename, pathname] = uiputfile( ...
    {'*.xls;';'*.xlsx'}, ...
    'Save As');
if isequal(filename,0) || isequal(pathname,0)
    disp('User cancelled saving excel data')
else
warning('off','MATLAB:xlswrite:AddSheet')
xlswrite([pathname filename],GRAPH,'GRAPH','A1');
xlswrite([pathname filename],SUBJECT,'SUBJECT','A1');
xlswrite([pathname filename],GROUPS,'GROUPS','A1');
end










