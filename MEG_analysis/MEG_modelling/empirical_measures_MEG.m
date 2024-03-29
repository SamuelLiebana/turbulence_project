% Calculating the Std of the Reduced Edge Matrix from empirical MEG data
% 
%
%  08.04.2021 Samuel Liebana, Morten Kringelbach and Gustavo Deco

% clean up env for clean run
clearvars; close all; clc
%% Parameters of the data, filtering and distance matrix

NSUB=89; 
NPARCELS=68;
Tmax=44196;
lambda=0.18; % coefficient of Markov Distance Rule
TR = 1/200;

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 8;                  % lowpass frequency of filter (Hz)
fhi = 13;                   % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter

% useful to extract the lower triangular portion of a matrix
Isubdiag = find(tril(ones(NPARCELS),-1));

% load the centres of gravity (saved as 'cog_schaeferNPARCELS.mat')
load(strcat('cog_dkt', num2str(NPARCELS),'.mat'))

% calculate the euclidean distance between parcels and store in rr matrix
rr = zeros(NPARCELS, NPARCELS);
for i=1:NPARCELS
    for j=1:NPARCELS
        rr(i,j)=sqrt((cog(i,1)-cog(j,1))^2+(cog(i,2)-cog(j,2))^2+(cog(i,3)-cog(j,3))^2);
    end
end
%% find the indexes of the edges with different lengths

% get edges
[u,v] = find(triu(ones(NPARCELS),1));    
lengths = zeros(length(u), 1);

for i=1:length(u)
    % lists the lengths of all the edges
    lengths(i)=floor(rr(u(i),v(i)));
end
% sorts the length array and returns the set of indexes 
% so that we can recreate the sort
[adistance,idx2]=sort(lengths);

% finds the indexes (-1) of all the *different* length edges in the network
dffidx100 = find(diff(adistance));
%% DEFINING VARIABLES FOR MEASUREMENTS

% matrix with amplitudes of local kuramoto order paramter 
% for each parcel over time
enstrophy=zeros(NPARCELS, Tmax);
enstrophy_su=zeros(NPARCELS,Tmax);

% temporary variable to store filtered signal
signal_filt=zeros(NPARCELS, Tmax);

% temporary variable to store hilbert transform of 
% filtered signal
Phases=zeros(NPARCELS, Tmax);
Phases_su=zeros(NPARCELS,Tmax);

% reduced edge time series matrix (phase version)
edge_phases=zeros(length(dffidx100)-1, Tmax);
edge_phases_su=zeros(length(dffidx100)-1, Tmax);
% intermediate variable used to store the instantaneous difference in phase 
% between hilbert transforms of node timeseries (i.e. edge time series)
phase_lock_matrix_red=zeros(1, length(dffidx100)-1);

% MEASURES TO BE USED FOR FITTING
% std of amplitude of local kuramoto order parameter
Rspatime=zeros(1, NSUB);
Rspatime_su=zeros(1,NSUB);

% mean and std of the global kuramoto order parameter
mean_global=zeros(1, NSUB);
mean_global_su=zeros(1,NSUB);

std_global=zeros(1, NSUB);
std_global_su=zeros(1,NSUB);

% std of the edge matrix
edge_matrix_std=zeros(1, NSUB);
edge_matrix_std_su=zeros(1, NSUB);

% mean of the correlation matrix for correlations between edge time series
% with different shifts
EdgeSpaTimePredictability=zeros(NSUB,8);
EdgeSpaTimePredictability_su=zeros(NSUB,8);

edgespa=zeros(NSUB,size(edge_phases, 1));
edgetime=zeros(NSUB,Tmax);
acfspa=zeros(NSUB,floor((length(dffidx100)-1)/1.2));
acftime=zeros(NSUB,3001);

edgespa_su=zeros(NSUB,size(edge_phases_su, 1));
edgetime_su=zeros(NSUB,Tmax);
acfspa_su=zeros(NSUB,floor((length(dffidx100)-1)/1.2));
acftime_su=zeros(NSUB,3001);


%% Loading the coupling matrix

load SC_dk68_Petra.mat 
C=C/max(max(C));


C1=C;
for i=1:NPARCELS
    C1(i,i)=1;
end
%% Load the data and compute the measures

% if the data is stored in individual files in a directory
myDir = strcat(pwd, '\data');	% gets directory
myFiles = dir(fullfile(myDir,'*.mat'));	

% otherwise if all subjects are in the same file
%load hcp1003_schaefer100_REST1_LR_all.mat

for sub = 1:NSUB
    sub
	baseFileName = myFiles(sub).name;
	fullFileName = fullfile(myDir, baseFileName);
	
    %clear schaeferts
    load(fullFileName);
    
    ts=X';
    ts = double(ts(:, 1:Tmax));
    
    % load the timeseries for each parcel
    %ts=subject{sub}.schaeferts;
    
    % filter the timeseries and find the hilbert transform
    for seed=1:NPARCELS
        ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
        signal_filt(seed,:)=filtfilt(bfilt,afilt,ts(seed,:));
        Xanalytic = hilbert(detrend(signal_filt(seed,:), 'constant'));
        Phases(seed,:) = angle(Xanalytic);
        Phases_su(seed,:)=Phases(seed,randperm(Tmax));
    end
    
    Phases_su=Phases_su(randperm(NPARCELS),:);
    
    % sum all the complex phases over parcels, weighted by the coupling
    % matrix. Take the abs of the result and store the time series in the
    % enstrophy matrix (original and surrogate)
    for i=1:NPARCELS
        sumphases=sum(repmat(C1(i,:)',1,Tmax).*complex(cos(Phases),sin(Phases)), 'omitnan')/sum(C1(i,:), 'omitnan');
        enstrophy(i,:)=abs(sumphases);
        sumphases=sum(repmat(C1(i,:)',1,Tmax).*complex(cos(Phases_su),sin(Phases_su)), 'omitnan')/sum(C1(i,:), 'omitnan');
        enstrophy_su(i,:)=abs(sumphases);
    end
    
    % calculate the global order parameter as the sum of phases without the C kernel    
    global_order_parameter =abs(nansum(complex(cos(Phases),sin(Phases)),1))/NPARCELS;
    global_order_parameter_su =abs(nansum(complex(cos(Phases_su),sin(Phases_su)),1))/NPARCELS;
    
    % two of the measures we will be using to evaluate fitting
    mean_global(sub) = mean(global_order_parameter, 'omitnan');
    std_global(sub) = std(global_order_parameter, 'omitnan');
    
    mean_global_su(sub) = mean(global_order_parameter_su, 'omitnan');
    std_global_su(sub) = std(global_order_parameter_su, 'omitnan');
    
    % std of the abs(local kuramoto order parameter) - another measure
    Rspatime(sub)=std(enstrophy(:));
    Rspatime_su(sub)=std(enstrophy_su(:));
    
    % here we calculate the edge time series (or edge_phases) as the
    % difference in phase between the different node pairs in network
    for n=1:Tmax
        cphi=cos(Phases(:,n));
        sphi=sin(Phases(:,n));
        restacos=repmat(cphi,1,NPARCELS)-repmat(cphi',NPARCELS,1);
        restasin=repmat(sphi,1,NPARCELS)-repmat(sphi',NPARCELS,1);
        PhLoMa=sqrt((restacos).^2+(restasin).^2);
        phase_lock_matrix=PhLoMa(Isubdiag);
        phase_lock_matrix = phase_lock_matrix(idx2); %%%%%%%%%%%%%%%%%%%%
        for i = 1:length(dffidx100)-1
            phase_lock_matrix_red(i)=mean(phase_lock_matrix(dffidx100(i)+1:dffidx100(i+1)), 'omitnan');
        end
        edge_phases(:,n)=phase_lock_matrix_red;
    end
    
    for n=1:Tmax
        cphi=cos(Phases_su(:,n));
        sphi=sin(Phases_su(:,n));
        restacos=repmat(cphi,1,NPARCELS)-repmat(cphi',NPARCELS,1);
        restasin=repmat(sphi,1,NPARCELS)-repmat(sphi',NPARCELS,1);
        PhLoMa=sqrt((restacos).^2+(restasin).^2);
        phase_lock_matrix=PhLoMa(Isubdiag);
        phase_lock_matrix = phase_lock_matrix(idx2); %%%%%%%%%%%%%%%%%%%%
        for i = 1:length(dffidx100)-1
            phase_lock_matrix_red(i)=mean(phase_lock_matrix(dffidx100(i)+1:dffidx100(i+1)), 'omitnan');
        end
        edge_phases_su(:,n)=phase_lock_matrix_red;
    end
    
    
    % fitting measure i.e. std of reduced edge matrix
    edge_matrix_std(sub)=std(edge_phases(:));
    edge_matrix_std_su(sub)=std(edge_phases_su(:));
    
    % if needed, we also calculate the edge spacetime predictability
    tot=length(dffidx100)-1;
    Isubdiagtot = find(triu(ones(tot),-1));
    gbindist=edge_phases';
    ncomm=zeros(tot,tot);
    for k=0:7
        gbindist_k=gbindist(1:end-k,:);
        gbindist_k2=gbindist(1+k:end,:);
        ncomm=corr(gbindist_k,gbindist_k2,'Rows','pairwise');
%        for i = 1:tot
%            for j = 1:tot
%                ccaux=corrcoef(gbindist(1:end-k,i),gbindist(1+k:end,j),'Rows','pairwise');
%                ncomm(i,j)=ccaux(2);
%                old_ncomm(sub,k+1,i,j)=ccaux(2);
%            end
%        end
       EdgeSpaTimePredictability(sub,k+1)=mean(ncomm(Isubdiagtot), 'omitnan');
    end
    
    gbindist=edge_phases_su';
    ncomm=zeros(tot,tot);
    for k=0:7
        gbindist_k=gbindist(1:end-k,:);
        gbindist_k2=gbindist(1+k:end,:);
        ncomm=corr(gbindist_k,gbindist_k2,'Rows','pairwise');
%        for i = 1:tot
%            for j = 1:tot
%                ccaux=corrcoef(gbindist(1:end-k,i),gbindist(1+k:end,j),'Rows','pairwise');
%                ncomm(i,j)=ccaux(2);
%                old_ncomm(sub,k+1,i,j)=ccaux(2);
%            end
%        end
       EdgeSpaTimePredictability_su(sub,k+1)=mean(ncomm(Isubdiagtot), 'omitnan');
    end
    
    edgespa(sub,:)=nanstd(edge_phases,[],2)';
    edgetime(sub,:)=nanstd(edge_phases,[],1);

    acfspa(sub,:)=autocorr(edgespa(sub,:),floor((length(dffidx100)-1)/1.2)-1);
    acftime(sub,:)=autocorr(edgetime(sub,:),3000);

    edgespa_su(sub,:)=nanstd(edge_phases_su,[],2)';
    edgetime_su(sub,:)=nanstd(edge_phases_su,[],1);

    acfspa_su(sub,:)=autocorr(edgespa_su(sub,:),floor((length(dffidx100)-1)/1.2)-1);
    acftime_su(sub,:)=autocorr(edgetime_su(sub,:),3000);
    
    
end

% for debug only
Rspatime_av = Rspatime;%mean(Rspatime);
mean_global_av = mean_global;%mean(mean_global);
std_global_av = std_global;%mean(std_global);
edge_matrix_std_av = edge_matrix_std;%mean(edge_matrix_std);
EdgeSpaTimePredictability_av = mean(EdgeSpaTimePredictability, 2, 'omitnan'); %mean(mean(EdgeSpaTimePredictability, 'omitnan'), 'omitnan');

Rspatime_av_su = Rspatime_su;%mean(Rspatime);
mean_global_av_su = mean_global_su;%mean(mean_global);
std_global_av_su = std_global_su;%mean(std_global);
edge_matrix_std_av_su = edge_matrix_std_su;%mean(edge_matrix_std);
EdgeSpaTimePredictability_av_su = mean(EdgeSpaTimePredictability_su, 2, 'omitnan'); %mean(mean(EdgeSpaTimePredictability, 'omitnan'), 'omitnan');








% save all the measures in HCP1000_empirical_100_subs.mat
save filtered_matrix_dk68_empirical_89_subs_MEG.mat Rspatime_av mean_global_av std_global_av edge_matrix_std_av EdgeSpaTimePredictability_av;
save filtered_su_matrix_dk68_empirical_89_subs_MEG.mat Rspatime_av_su mean_global_av_su std_global_av_su edge_matrix_std_av_su EdgeSpaTimePredictability_av_su;


Turbu=nanmean(edge_matrix_std_av);
Turbu_su=nanmean(edge_matrix_std_av_su);


p=ranksum(edge_matrix_std_av,edge_matrix_std_av_su)

figure(1)
boxplot([edge_matrix_std' edge_matrix_std_su']);

%% autocorr space
figure(2)
shadedErrorBar(0:floor((length(dffidx100)-1)/1.2),[1, nanmean(acfspa)],[0, nanstd(acfspa)],'lineProps', '-r', 'patchSaturation', 0.2)
%hold on;
%shadedErrorBar(0:floor(size(edge_phases, 1)/1.2),[1, nanmean(acfspa_su)],[0, nanstd(acfspa_su)],'lineProps', '-k','patchSaturation', 0.2)

%% autocorr time 
figure(3)
shadedErrorBar(0:3001,[1, nanmean(acftime)],[0, nanstd(acftime)],'lineProps', '-r','patchSaturation', 0.2)
%hold on;
%shadedErrorBar(0:3001,[1, nanmean(acftime_su)],[0, nanstd(acftime_su)],'lineProps', '-k','patchSaturation', 0.2)


%% check:
for i=2000:5000:42000;%2:600
    i
    figure()
    plot(Phases(:,i),'r*');
    figure()
    plot(Phases_su(:,i),'k*');
end


%% Mean information cascade predictability oevr 8 time shifts

%load turbu_emp.mat

%figure(1);

%plot(EdgeSpaTimePredictability(1,:), 'Linewidth',2); hold on;

%title(['Mean Information cascade predictability']);
%% Information cascade predictability for each time shift

%figure(2);

%old_ncomm = squeeze(mean(old_ncomm, 1, 'omitnan'));

%for k=0:7
%    for i = 1:tot-1
%        ncomm2{i}=old_ncomm(k+1,i,i:end);
%    end;
    
%    tot=length(dffidx100)-1;
%    for i = 1:tot
%        tmp=[];
%        for j = 1:tot-i
%            tmp(j)=ncomm2{j}(i+1);
%        end;
%        avdiagncomm(i)=mean(tmp, 'omitnan');
%    end;
%    plot(avdiagncomm,'Linewidth',2); hold on;
    
%end;
%title(['Information cascade predictability']);
