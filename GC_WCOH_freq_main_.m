%% MVGC "GCCA compatibility mode" demo
%
% Demonstrates usage of the MVGC toolbox in "GCCA compatibility mode"; see
% <mvgchelp.html#6 Miscellaneous issues> in the Help documentation. This is
% partly for the benefit of former users of the Granger Causal Connectivity
% Analysis (<http://www.sussex.ac.uk/Users/anils/aks_code.htm GCCA>) Toolbox
% [2], and partly as an implementation of a more "traditional" approach to
% Granger causality computation. The chief difference is that here two separate
% VAR regressions - the _full_ and _reduced_ regressions (see [1]) - are
% explicitly performed (see <GCCA_tsdata_to_pwcgc.html
% |GCCA_tsdata_to_pwcgc|>), in contrast to the MVGC Toolbox preferred
% approach (see <mvgc_demo.html |mvgc_demo|>), which only requires a full
% regression and is consequently more flexible and numerically accurate.
%
% Granger-causal pairwise-conditional analysis is demonstrated on generated
% VAR data for a 5-node network with known causal structure (see
% <var5_test.html |var5_test|>), as in the main MVGC Toolbox demonstration
% script, <mvgc_demo.html |mvgc_demo|>. A drawback of the traditional dual
% regression approach is that in the frequency domain, _conditional_
% spectral causalities cannot be estimated to an acceptable standard; see
% [1] and <GCCA_tsdata_to_smvgc.html |GCCA_tsdata_to_smvgc|> for more
% detail.
%
%% setups
clc;clear all; close all;
% add paths for the data folder and support files for the toolbox
addpath D:\RPI\ResearchWork\Papers_\Data\Data_processed_BU\FLS
addpath D:\RPI\ResearchWork\Papers_\Data\Data_processed_BU\VBLAST
addpath D:\RPI\ResearchWork\Papers_\Effective_Connectivity
addpath D:\RPI\ResearchWork\Papers_\Effective_Connectivity\GCCA
addpath D:\RPI\ResearchWork\Papers_\Effective_Connectivity\GCCA\utilities
addpath D:\RPI\ResearchWork\Papers_\Effective_Connectivity\bsmart
addpath D:\RPI\ResearchWork\Papers_\Arun_data_codes\VBLAST
addpath D:\RPI\ResearchWork\Papers_\Arun_data_codes\FLS
addpath D:\RPI\ResearchWork\Papers_\Effective_Connectivity\toolbox_original\mvgc_v1.0\demo
names = {'LPFC' 'RPFC' 'LPMC' 'RPMC' 'SMA'};
fpath = 'D:\\RPI\\ResearchWork\\Papers_\\Effective_Connectivity\\toolbox_original\\mvgc_v1.0\\Plots\\FLSvsVBLAST';
connection_GC = {'LPFC=>RPFC','LPFC=>LPMC','LPFC=>RPMC','LPFC=>SMA',...
                  'RPFC=>LPFC','RPFC=>LPMC','RPFC=>RPMC','RPFC=>SMA',...
                  'LPMC=>LPFC','LPMC=>RPFC','LPMC=>RPMC','LPMC=>SMA',...
                  'RPMC=>LPFC','RPMC=>RPFC','RPMC=>LPMC','RPMC=>SMA',...
                  'SMA=>LPFC','SMA=>RPFC','SMA=>LPMC','SMA=>RPMC'};
%% #################  WCOH analysis ################
% ##################################################
nos.FLS = 8;
nos.VBLAST = 6;
WCOH = [];
for m = 1:2  % 1 for FLS and 2 for VBLAST
    if m ==1
        N = nos.FLS;
        T = 1;     % days selected for analysis; maximum common in all fls = 10 days
    else
        N = nos.VBLAST;
        T = 1;      % days selected for analysis; maximum common in all fls = 5 days
    end
    for s = 1:N  %no of Subjects
        for t = 1:T    % Days
            if m == 1
                sub_day = sprintf('FLS%d-%d.mat',s,t);
                sub_day_nirs = sprintf('FLS%d-%d.nirs',s,t);
                sub_day_name = sprintf('FLS%d_%d',s,t);
            else
                sub_day = sprintf('VBLAST%d-%d.mat',s,t);
                sub_day_nirs = sprintf('VBLAST%d-%d.nirs',s,t);
                sub_day_name = sprintf('VBLAST%d_%d',s,t);
            end
            if exist(sub_day , 'file') == 2
                DD = load(sub_day);
                da = load(sub_day_nirs,'-mat');
            else
                continue;
            end
            w = find(da.s);    % for trials, stimulus trigger points
            Data = (DD.output.dc.dataTimeSeries)'; % (:,1:33)
            k = 1;
            D = [];
            for i = 1:3:96
                D(k,:) = Data(i,:);
                k      = k+1;
            end
            D = D(:,w(1):w(2));
            LPFC    = D([1,2],:);
            LPFC( any(isnan(LPFC),2),:) = [];
            LPFC    = mean(LPFC,1);
            RPFC    = D([7,8],:);
            RPFC( any(isnan(RPFC),2),:) = [];
            RPFC    = mean(RPFC,1);
            LPMC    = D([10 11],:);                 % with only two channels LPMC  = D([10 11],:); for more chnls D([10 11 12 13],:)
            LPMC( any(isnan(LPMC),2),:) = [];
            LPMC    = mean(LPMC,1);
            RPMC    = D([27 28],:);                        % with only two extreme channels RPMC    = D([27 28],:); or more chnl D([25 26 27 28]
            RPMC( any(isnan(RPMC),2),:) = [];
            RPMC    = mean(RPMC,1);
            SMA     = D([29 30 31],:);
            SMA     = D([32],:);
            Data_Reg = [LPFC;RPFC;LPMC;RPMC;SMA];
            valid_Ch_names = names;
            NaN_rows = find(all(isnan(Data_Reg),2));
            valid_Ch_names(:,NaN_rows) = [];
            Data_Reg( any(isnan(Data_Reg),2),:) = [];
            for i = 1:5
                for j = i+1:5
                    wcoh    = [];
                    F       = [];   % frequency
                    coi     = [];
                    if i ~= j
                        fprintf('i=%d j =%d\n',i,j);
                        [wcoh,~,F,coi] = wcoherence(Data_Reg(i,:),Data_Reg(j,:),25);    %Fs = 25Hz
                        [~, idx1]      = min(abs(F-0.1));
                        [~, idx2]      = min(abs(F-0.07));
                        WCO            = nanmean(nanmean(wcoh(idx1:idx2,:)));
                        connection      = [valid_Ch_names{i},'_', valid_Ch_names{j}];
                        if m ==1
                            WCOH.FLS.(sub_day_name)                        = wcoh;       % for all frequency
                            WCOH_f1mean.FLS.(connection).(sub_day_name)    = WCO;        % for average of the f1 frequency range.
                        else
                            WCOH.VBLaST.(sub_day_name)                     = wcoh;    % for all frequency
                            WCOH_f1mean.VBLaST.(connection).(sub_day_name) = WCO;        % for average of the f1 frequency range.
                        end
                        %                         close all;
                        %                         figure(5)
                        %                         ax = gca;
                        %                         wcoherence(Data_Reg(2,:),Data_Reg(3,:),20);
                        %                         tt = sprintf('%s-%s-%s-WCOH ',sub_trial,valid_Ch_names{i},valid_Ch_names{j});
                        %                         title(tt)
                        %                         ax.XLabel.String= 'Time(sec)';
                        %                         baseFileName = sprintf('%s.png',tt);
                        %                         fullFileName = fullfile(fpath, baseFileName);
                        %                         saveas(figure(5),fullFileName)
                    end
                end
            end
            fprintf(' Analyzing: %s \n',sub_day);
        end
    end
end
%% Physical and virtual score
score = xlsread('FLS_scores.xlsx');
sc_phy = score(10:17,1:10);
sc_vir = score(1:6,1:10);
mean_phy = nanmean(sc_phy,1);
mean_vir = nanmean(sc_vir,1);
figure()
plot(mean_phy,'LineWidth',2);hold on;
plot(mean_vir,'LineWidth',2);hold on;
title('Physical and virtual score: mean of each day');
legend('mean physical\_FLS','mean Virtual\_FLS')
ylabel('FLS score');
xlabel('day');
% figure()
% plot(sc_phy(:,1),'LineWidth',2);hold on;
% plot(sc_vir(:,1),'LineWidth',2); hold on;
% plot(sc_phy(:,2),'LineWidth',2); hold on;
% plot(sc_vir(:,2),'LineWidth',2);hold on;
% plot(sc_phy(:,3),'LineWidth',2);hold on;
% plot(sc_vir(:,3),'LineWidth',2);
% title('Physical and virtual score: first trail different days');
% legend('1D1T physical','1D1T Virtual','2D1T Physical','2D1T  virtual','3D1T Physical','3D1T trial virtual')
% ylabel('FLS score');
% xlabel('subjects');

% Anova on the mean score
[~,NORM_p.sc_phy] = kstest(sc_phy(:,1),'Alpha',0.01);
[~,NORM_p.sc_vir] = kstest(sc_vir(:,1),'Alpha',0.01);
sc_phy_vari = var(sc_phy);
sc_phy_ratio = sc_phy_vari./ mean_phy;
sc_vir_ratio = sc_vir_vari./ mean_vir;
sc_vir_vari = var(sc_vir);
figure()
plot(sc_phy_ratio,'LineWidth',2); hold on;
plot(sc_vir_ratio,'LineWidth',2);
title('Variance to mean ratio');
legend('Physical','Virtual');
xlabel('days');
ylabel('variance to mean');

for i =1:10
    X = [sc_phy(1:6,i),sc_vir(:,i)];
    sc_anova = anova1(X);
end

%% ANOVA FLS vs VBLaST
nos.FLS = 8;  % no of FLS subject 8
nos.VBLAST = 6;  % no of FLS subject 6
FLS_wcoh = [];
VBLaST_wcoh = [];
connection = [];
subject_FLS = {};
subject_VBLAST = {};
for m = 1:2  % 1 for FLS and 2 for VBLAST
    if m ==1
        N = nos.FLS;
        T = 1;     % minimum no of common day in all subjects to include for ANOVA ex: set 1 for first day
    else
        N = nos.VBLAST;
        T = 1;
    end
    r       = 1;    % loop for row or subj_trial
    for s = 1:N  %no of Subjects
        for t = 1:T    %days
            if m == 1
                sub_day = sprintf('FLS%d-%d.mat',s,t);
                sub_day_name = sprintf('FLS%d_%d',s,t);
            else
                sub_day = sprintf('VBLAST%d-%d.mat',s,t);
                sub_day_name = sprintf('VBLAST%d_%d',s,t);
            end
            c       = 1;    % loop for column or the connections
            for i = 1:5
                for j = i+1:5
                    wcoh        = [];
                    F           = [];   % frequency
                    coi         = [];
                    if i ~= j
                        connection{c}         = [valid_Ch_names{i},'_', valid_Ch_names{j}];
                        if m ==1
                            FLS_wcoh(r,c)     = WCOH_f1mean.FLS.(connection{c}).(sub_day_name);
                            subject_FLS{r,1}       = sprintf('FLS%d',s);
                        else
                            VBLaST_wcoh(r,c)  = WCOH_f1mean.VBLaST.(connection{c}).(sub_day_name);
                            subject_VBLAST{r,1}       = sprintf('VBLAST%d',s);
                        end
                        c = c+1;
                    end
                end
            end
            r = r+1;
            fprintf(' Analyzing: %s \n',sub_day);
        end
    end
end
%% Saphiro-Wilk's test of Normality
[H, pValue, W] = swtest(FLS_wcoh(:,1), 0.05);

for i = 1:10
    [~,NORM_p.FLS(i,1)] = kstest(FLS_wcoh(:,i));
end
for i = 1:10
    [~,NORM_p.VBLaST(i,1)] = kstest(VBLaST_wcoh(:,i));
end

%Two way anova
g = [ones(1,size(FLS_wcoh,1))'; zeros(1,size(VBLaST_wcoh,1))'];
p = anovan([FLS_wcoh(:,1);VBLaST_wcoh(:,1)],{g});

%With in FLS;
[p1,tbl1,stats1] = anova1(FLS_wcoh);
title('One way ANOVA between all connection in FLS')
ylabel('FLS-Average WCOH for freq[0.1,0.07]')
xticklabels(connection);
xtickangle(45)

%With in VBLaST;
[p2,tbl2,stats2] = anova1(VBLaST_wcoh);
title('One way ANOVA between all connection in VBLaST')
ylabel('VBLaST-Average WCOH for freq[0.1,0.07]')
xticklabels(connection);
xtickangle(45)
%% Two way anova test
n = size(VBLaST_wcoh,1);
X = [FLS_wcoh(1:n,:); VBLaST_wcoh(1:n,:)];
[~,~,stats] = anova2(X,n);
%% Repeated measure Anova test
% For FLS
tb1 = table(subject_FLS,FLS_wcoh(:,1),FLS_wcoh(:,2),FLS_wcoh(:,3),FLS_wcoh(:,4),FLS_wcoh(:,5),FLS_wcoh(:,6),...
    FLS_wcoh(:,7),FLS_wcoh(:,8),FLS_wcoh(:,9),FLS_wcoh(:,10),...
    'VariableNames',{'subject_FLS','C1','C2','C3','C4','C5','C6','C7','C8','C9','C10'});
Meas = table([1 2 3 4 5 6 7 8 9 10]','VariableNames',{'connection'});
rm = fitrm(tb1,'C1-C10~subject_FLS','WithinDesign',Meas);
ranovatbl = ranova(rm);
%for VBLAST
tb2 = table(subject_VBLAST,VBLaST_wcoh(:,1),VBLaST_wcoh(:,2),VBLaST_wcoh(:,3),VBLaST_wcoh(:,4),VBLaST_wcoh(:,5),...
    VBLaST_wcoh(:,6),VBLaST_wcoh(:,7),VBLaST_wcoh(:,8),VBLaST_wcoh(:,9),VBLaST_wcoh(:,10),...
    'VariableNames',{'subject_FLS','C1','C2','C3','C4','C5','C6','C7','C8','C9','C10'});
Meas = table([1 2 3 4 5 6 7 8 9 10]','VariableNames',{'connection'});
rm = fitrm(tb2,'C1-C10~subject_FLS','WithinDesign',Meas);
ranovatbl = ranova(rm);
% for both
tb = [tb1;tb2];
rm = fitrm(tb,'C1-C10~subject_FLS','WithinDesign',Meas);
ranovatbl = ranova(rm);

%% Anova for same connection between FLS and VBLaST
group = {'FLS';'VBLaST'};
n = size(VBLaST_wcoh,1);
for i = 1:10
    %close all;
    X = [FLS_wcoh(1:n,i), VBLaST_wcoh(1:n,i)];
    tbl = ['tbl',num2str(i)];
    %     figure(10)
    baseFileName = sprintf('FLS vs VBLaST Aonva Conn:-%s.png',connection{i});
    [p1(i),tbl,stats.(connection{i})] = anova1(X);
    title(baseFileName)
    xticklabels(group);
    xtickangle(45)
    fullFileName = fullfile(fpath, baseFileName);
    %saveas(figure,fullFileName)
end
%       ANOVA on trial basis
% for m = 1:1  % 1 for FLS and 2 for VBLAST
%     if m ==1
%         N = nos.FLS;
%         T = 1;     % minimum no of common day in all subjects to include for ANOVA ex: set 1 for first trial
%     else
%         N = nos.VBLAST;
%         T = 1;
%     end
%     r       = 1;    % loop for row or subj_day
%     for s = 1:N  %no of Subjects
%         for t = 1:T    %day
%             if m == 1
%                 sub_day = sprintf('FLS%d-%d.mat',s,t);
%                 sub_day_name = sprintf('FLS%d_%d',s,t);
%             else
%                 sub_day = sprintf('VBLAST%d-%d.mat',s,t);
%                 sub_day_name = sprintf('VBLAST%d_%d',s,t);
%             end
%             c       = 1;    % loop for column or the connections
%             for i = 1:5
%                 for j = i+1:5
%                     wcoh    = [];
%                     F       = [];   % frequency
%                     coi     = [];
%                     if i ~= j
%                         connection{c}      = [valid_Ch_names{i},'_', valid_Ch_names{j}];
%                         if m ==1
%                             FLS_wcoh(r,c)     = WCOH_f1mean.FLS.(connection{c}).(sub_day_name);
%                         else
%                             VBLaST_wcoh(r,c) = WCOH_f1mean.VBLaST.(connection{c}).(sub_trial_name);
%                         end
%                         c = c+1;
%                     end
%                 end
%             end
%             r = r+1;
%             fprintf(' Analyzing: %s \n',sub_day);
%         end
%     end
% end
%% #############################################GC Analysis########################################################%%%
% ############################################
%load('FLS1-1.nirs')
nos.FLS = 8;  % no of FLS subject
nos.VBLAST = 6;  % no of FLS subject
GC_f1 = [];
GC_f2 =[];
B_phy = zeros(nos.FLS,25);               % 25 for all pair of connections
B_vir = zeros(nos.VBLAST,25);
h = 1; hh = 1;  % counters 
for m = 1:2  % 1 for FLS and 2 for VBLAST
    if m ==1
        N = nos.FLS;
    else
        N = nos.VBLAST;
    end
    for s = 1:N  %no of Subjects
        for t = 1:1    %Days
            if m == 1
                sub_day      = sprintf('FLS%d-%d.mat',s,t);
                sub_day_nirs = sprintf('FLS%d-%d.nirs',s,t);
                sub_day_name = sprintf('FLS%d_%d',s,t);
            else
                sub_day      = sprintf('VBLAST%d-%d.mat',s,t);
                sub_day_nirs = sprintf('VBLAST%d-%d.nirs',s,t);
                sub_day_name = sprintf('VBLAST%d_%d',s,t);
            end
            if exist(sub_day , 'file') == 2
                DD = load(sub_day);
                da = load(sub_day_nirs,'-mat');
            else
                continue;
            end
            w = find(da.s);    % for trials, stimulus trigger points
            Data = (DD.output.dc.dataTimeSeries)'; % (:,1:33)
            k = 1;
            D = [];
            for i = 1:3:96
                D(k,:) = Data(i,:);
                k      = k+1;
            end
            D = D(:,w(1):w(2));
            LPFC    = D([1,2],:);
            LPFC( any(isnan(LPFC),2),:) = [];
            LPFC    = mean(LPFC,1);
            RPFC    = D([7,8],:);
            RPFC( any(isnan(RPFC),2),:) = [];
            RPFC    = mean(RPFC,1);
            LPMC    = D([10 11 12 13],:);                 % with only two channels LPMC  = D([10 11],:);
            LPMC( any(isnan(LPMC),2),:) = [];
            LPMC    = mean(LPMC,1);
            RPMC    = D([25 26 27 28],:);                        % with only two extreme channels RPMC    = D([27 28],:);
            RPMC( any(isnan(RPMC),2),:) = [];
            RPMC    = mean(RPMC,1);
            SMA     = D([29 30 31],:);
            SMA     = D([32],:);
            Data_Reg = [LPFC;RPFC;LPMC;RPMC;SMA];
            valid_Ch_names = names;
            NaN_rows = find(all(isnan(Data_Reg),2));
            valid_Ch_names(:,NaN_rows) = [];
            Data_Reg( any(isnan(Data_Reg),2),:) = [];
            %         figure(2)
            %         subplot(2,1,1)
            %         plot(D')
            %         title('Signals')
            %         xlabel('timesteps')
            %         ylabel('Magnitude')
            %         subplot(2,1,2)
            %         plot(Data_Reg')
            %         title('Signals after averaging each region')
            %         xlabel('timesteps')
            %         ylabel('Magnitude')
            
            SMA_freq = fft(SMA);
            %         n = length(SMA);
            %         power_SMA = abs(SMA_freq).^2/n;
            pow_SMA = pwelch(SMA);
            figure(3)
            subplot(2,1,1)
            plot(SMA_freq)
            title('FFT of a signalSMA')
            xlabel('Frequency')
            ylabel('Amplitude')
            subplot(2,1,2)
            plot(pow_SMA)
            title('PSD of the signal SMA')
            xlabel('Frequency')
            ylabel('Power Per frequency')
            baseFileName = sprintf('%s_FFT_PSD.png',sub_day);
            fullFileName = fullfile(fpath, baseFileName);
            %saveas(figure(3),fullFileName)
            
            k = rank(D'); % rank
            Z = null(D'); % null space
            %[m,n,M] = size(D);
            
            %[z, pvalue] = kpsstest(Ynew(:,1)');
            freq1 = [0:0.01:2];                                            %Frequency range for GC calculation
            freq2 = [0:0.01:0.5];                                            %Frequency range for GC calculation
            Y = Data_Reg;%(:,100:1109);
            for f = 1:2
                if f == 1
                    freq = freq1;
                else
                    freq = freq2;
                end
                NL = size(Y,2);
                [GW,COH,pp,waut,cons]= cca_pwcausal(Y,1,NL,20,20,freq, 1);
                if f == 1
                    GC_f1.(sub_day_name) = GW;     % for freq1
                else
                    GC_f2.(sub_day_name) = GW;      % for freq2
                end
                %mean in the neurophysiology frequency band.
                idx1 = find(freq == 0.07);
                idx2 = find(freq == 0.1);
                if f == 2
                    if m == 1
                        GC_fqmean.FLS.(sub_day_name)   = mean(GW(:,:,idx1:idx2),3);
                        sz = numel(GC_fqmean.FLS.(sub_day_name));
                        B_phy(h,:)  = reshape(GC_fqmean.FLS.(sub_day_name),1,sz);
                        h = h+1;
                    else
                        GC_fqmean.VBLaST.(sub_day_name) = mean(GW(:,:,idx1:idx2),3);
                        sz = numel(GC_fqmean.VBLaST.(sub_day_name));
                        B_vir(hh,:)  = reshape(GC_fqmean.VBLaST.(sub_day_name),1,sz);
                        hh = hh+1;
                    end
                end
                figure(5); clf reset;       %Plot GC
                ttl = sprintf('G-Causality for Frequecy range %d',f);
                title(ttl);
                cca_plotcausality_spectral(GW,freq);
                baseFileName = sprintf('%s_GC.png',sub_day);
                fullFileName = fullfile(fpath, baseFileName);
                %saveas(figure(5),fullFileName)
                
                % plot coherence
                figure(6); clf reset;
                ttl = sprintf('Coherence for Frequecy range %d',f);
                title(ttl);
                cca_plotcoherence(COH,freq);
            end
            fprintf('%s \n',sub_day);
        end
    end
end
B_phy(:,[1 7 13 19 25]) = [];       % Remove all the diagonals 
B_vir(:,[1 7 13 19 25]) = [];
%% ################  GC--ANOVA test  ############
% KS normality test
for i = 1:20
    [~,NORM_p.FLS(i,1)] = kstest(B_phy(:,i));
end
for i = 1:20
    [~,NORM_p.VBLaST(i,1)] = kstest(B_vir(:,i));
end
%With in FLS;
[p11,tbl1,stats11] = anova1(B_phy);
title('One way ANOVA between all GC_connection in Physical')
ylabel('GC freq[0.1,0.07]')
xticklabels(connection_GC);
xtickangle(45)

%With in VBLaST;
[p22,tbl2,stats22] = anova1(B_vir);
title('One way ANOVA between all GC_connection in virtual')
ylabel('GC freq[0.1,0.07]')
xticklabels(connection_GC);
xtickangle(45)

% pair wise anova test between Physical and virtual simulator
group = {'FLS';'VBLaST'};
n = size(B_vir,1);
X = [];
for i = 1:20
    %close all;
    X = [B_phy(1:n,i), B_vir(1:n,i)];
    tbl = ['tbl',num2str(i)];
    %     figure(10)
    baseFileName = sprintf('FLS vs VBLaST Aonva Conn:-%s.png',connection_GC{i});
    [p1(i),tbl,stats] = anova1(X);
    title(baseFileName)
    xticklabels(group);
    xtickangle(45)
    fullFileName = fullfile(fpath, baseFileName);
    %saveas(figure,fullFileName)
end

%% Mean and Median value
t_step = size(GC_f2.FLS1_1.Win_1,3); % maximum number of time steps at frq_2
GCmean = [];
GCmedian = [];
tr = 1;         % Days number
for m = 1:2     % for FLS and VBLAST
    if m == 1
        N = nos.FLS;
        n_win = 25;
    else
        N = nos.VBLAST;
        n_win = 10;
    end
    for w = 1:n_win                 % number of windows
        W = sprintf('Win_%d',w);
        for t = 1:t_step
            for i =1:5              % for each pair in the matrix
                for j = 1:5         % for each pair in the matrix
                    temp = [];      % stores vector of all the subject to find mean and median
                    for s = 1:N
                        if m ==1
                            sub_day_name = sprintf('FLS%d_%d',s,tr);
                        else
                            sub_day_name = sprintf('VBLAST%d_%d',s,tr);
                        end
                        temp(s) = GC_f2.(sub_day_name).(W)(i,j,t);
                    end
                    if m ==1
                        GCmean.FLS.(W)(i,j ,t) = mean(temp);
                        GCmedian.FLS.(W)(i,j,t) = nanmedian(temp);
                    else
                        GCmean.VBLAST.(W)(i,j ,t) = mean(temp);
                        GCmedian.VBLAST.(W)(i,j,t) = nanmedian(temp);
                    end
                end
            end
        end
    end
end
%%
freq = freq2;
for m = 1:2     % for FLS and VBLAST
    if m == 1
        n_win = 25;
    else
        n_win = 10;
    end
    for w = 1:n_win                 % number of windows
        W = sprintf('Win_%d',w);
        % plot FLS mean
        if m ==1
            figure(7); clf reset;       %Plot GC
            ttl = sprintf('GC mean for Frequecy range %d',f);
            title(ttl);
            cca_plotcausality_spectral(GCmean.FLS.(W),freq);
            baseFileName = sprintf('%s_FLS_GC_mean.png',W);
            fullFileName = fullfile(fpath, baseFileName);
            saveas(figure(7),fullFileName)
            figure(9); clf reset;       %Plot GC
            ttl = sprintf('G-Causality median for Frequecy range %d',f);
            title(ttl);
            cca_plotcausality_spectral(GCmedian.FLS.(W),freq);
            baseFileName = sprintf('%s_FLS_GC_median.png',W);
            fullFileName = fullfile(fpath, baseFileName);
            saveas(figure(9),fullFileName)
        else
            % plot VBLAST median
            figure(8); clf reset;       %Plot GC
            ttl = sprintf('GC mean for Frequecy range %d',f);
            title(ttl);
            cca_plotcausality_spectral(GCmean.VBLAST.(W),freq);
            baseFileName = sprintf('%s_VBLAST_GC_mean.png',W);
            fullFileName = fullfile(fpath, baseFileName);
            saveas(figure(8),fullFileName)
            figure(10); clf reset;       %Plot GC
            ttl = sprintf('G-Causality median for Frequecy range %d',f);
            title(ttl);
            cca_plotcausality_spectral(GCmedian.VBLAST.(W),freq);
            baseFileName = sprintf('%s_VBLAST_GC_median.png',W);
            fullFileName = fullfile(fpath, baseFileName);
            saveas(figure(10),fullFileName)
        end
        
    end
end
%% Bootstrap resampling
Nr = 1; % no of trial
Nl = size(Y_windowed,2);
nlags = 100;
nBoot = 100;
nBwin = 2340;   % window size; should be selected so as to create at least 10 window
pval =  0.01;
CORRTYPE = 0;   % bootstrap confidenceintervals are given for difference-of-influence terms (=1), or for partial difference-of-influence terms (=2) or not for at all for these terms (=0).
DOIFLAG = 0;
Fs = 25;        % sampling frequency

tic
[ret] = cca_pwcausal_bstrap(Y_windowed,Nr,Nl,nlags,nBoot,nBwin,Fs,freq,pval,CORRTYPE) ;
toc
figure(5); clf reset;       %Plot GC
ttl = sprintf('Bootstrap resampling G-Causality for Frequecy range %d',f);
title(ttl);
cca_plotcausality_spectral(GW, freq, ret.ll, ret.ul);

%% Permutation test
% [ret] = cca_pwcausal_permute(Y_windowed,Nr,Nl,nlags,nBoot,nBwin,Fs,freq,pval,CORRTYPE);
% figure(8); clf reset;       %Plot GC
% ttl = sprintf('Permutation test G-Causality for Frequecy range %d',f);
% title(ttl);
% cca_plotcausality_spectral(GW, freq, ret.st);


%% Wavelet transformation method
level = 6;
wpt = wpdec(Y_windowed(1,:),level,'sym6');
figure;
[SPEC,TIMES,FREQ] = wpspectrum(wpt,Fs,'plot');

%% filterings
%     figure(1)
%     subplot(2,2,1)
%     plot(RSMA)
%     title('Signal without filtering')
%     subplot(2,2,3)
%     plot(EXP_filt(5,:))
%     title('Signal with filtfilt filtering')
%     subplot(2,2,2)
%     plot(RSMA_freq)
%     title('Signal with filtfilt filtering')
% detrending and whitening
%Ynew = detrend(Ynew,1);
%Ynew = diff(Ynew',2,2);


% using filtfilt for filtering
%     d = designfilt('lowpassfir', ...
%         'PassbandFrequency',0.01,'StopbandFrequency',0.7, ...
%         'PassbandRipple',1,'StopbandAttenuation',60, ...
%         'DesignMethod','equiripple');
%     EXP_filt = (filtfilt(d,Data_Reg'))';
%NOTE:Also Use dynamic time wraping to compare two signal.