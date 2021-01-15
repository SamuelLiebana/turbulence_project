% Calculating the Mean of the Global Kuramoto Order Parameter, Std of the
% Global Kuramoto Order Parameter, Std of the Local Kuramoto Order
% Parameter and the Std of the Reduced Edge Matrix
% 
%
%  12.01.2021 Samuel Liebana, Morten Kringelbach and Gustavo Deco

% clean up env for clean run
clearvars; close all; clc
%% Parameters of the data, filtering and distance matrix

NSUB=1003; 
NPARCELS=100;
Tmax=1200;
lambda=0.18; % coefficient of Markov Distance Rule
TR=0.72;  % Repetition Time (seconds)

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                  % lowpass frequency of filter (Hz)
fhi = 0.08;                   % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter

% useful to extract the lower triangular portion of a matrix
Isubdiag = find(tril(ones(NPARCELS),-1));

% load the centres of gravity (saved as 'cog_schaeferNPARCELS.mat')
load(strcat('cog_schaefer', num2str(NPARCELS),'.mat'))

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
% for each parcel over time (original and surrogate)
enstrophy=zeros(NPARCELS,Tmax);
enstrophy_su=zeros(NPARCELS,Tmax);

% temporary variable to store filtered signal
signal_filt=zeros(NPARCELS,Tmax);

% temporary variable to store hilbert transform of 
% filtered signal (original and surrogate)
Phases=zeros(NPARCELS,Tmax);
Phases_su=zeros(NPARCELS,Tmax);

% reduced edge time series matrix (phase version)
edge_phases=zeros(length(dffidx100)-1,Tmax);
% intermediate variable used to store the instantaneous difference in phase 
% between hilbert transforms of node timeseries (i.e. edge time series)
phase_lock_matrix_red=zeros(1,length(dffidx100)-1);

% MEASURES TO BE USED FOR FITTING
% std of amplitude of local kuramoto order parameter (original and surrogate)
Rspatime=zeros(1,NSUB);
Rspatime_su=zeros(1,NSUB);

% mean and std of the global kuramoto order parameter
mean_global=zeros(1,NSUB);
std_global=zeros(1,NSUB);

% std of the edge matrix
edge_matrix_std=zeros(1,NSUB);
% mean of the correlation matrix for correlations between edge time series
% with different shifts
EdgeSpaTimePredictability=zeros(NSUB,8);
%% Loading the coupling matrix

% choose whether to use markov rule or DTI measurements
Markov =1;

if Markov==1
    C=zeros(NPARCELS,NPARCELS);
    for i=1:NPARCELS
        for j=1:NPARCELS
            C(i,j)=exp(-lambda*rr(i,j));
        end
        C(i,i)=1;
    end
else 
    load sc_schaefer.mat
    C=sc_schaefer;
    C=C/max(max(C));
end
%% Load the data and compute the measures

% if the data is stored in individual files in a directory
%myDir = strcat(pwd, '\schaefer1000_100participants');	% gets directory
%myFiles = dir(fullfile(myDir,'*.mat'));	

% otherwise if all subjects are in the same file
load hcp1003_schaefer100_REST1_LR_all.mat

for sub = 1:NSUB
    sub
	%baseFileName = myFiles(sub).name;
	%fullFileName = fullfile(myDir, baseFileName);
	
    %clear schaeferts
    %load(fullFileName);
    
    %ts=schaeferts;
    
    % load the timeseries for each parcel
    ts=subject{sub}.schaeferts;
    
    %removes all NaN rows from the matrix (555 and 908 were all NaNs in schaefer 1000)
    %ts = ts(all(~isnan(ts),2),:);
    
    % sets all NaNs to 0 - another option
    %ts(isnan(ts))=0;
    
    % filter the timeseries and find the hilbert transform
    for seed=1:NPARCELS
        ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
        signal_filt(seed,:)=filtfilt(bfilt,afilt,ts(seed,:));
        Xanalytic = hilbert(detrend(signal_filt(seed,:), 'constant'));
        Phases(seed,:) = angle(Xanalytic);
        Phases_su(seed,:)=Phases(seed,randperm(Tmax));
    end
    
    % sum all the complex phases over parcels, weighted by the coupling
    % matrix. Take the abs of the result and store the time series in the
    % enstrophy matrix (original and surrogate)
    for i=1:NPARCELS
        sumphases=sum(repmat(C(i,:)',1,Tmax).*complex(cos(Phases),sin(Phases)), 'omitnan')/sum(C(i,:), 'omitnan');
        enstrophy(i,:)=abs(sumphases);
        sumphases=sum(repmat(C(i,:)',1,Tmax).*complex(cos(Phases_su),sin(Phases_su)), 'omitnan')/sum(C(i,:), 'omitnan');
        enstrophy_su(i,:)=abs(sumphases);
    end
    
    % calculate the global order parameter as the mean over parcels of the
    % enstrophy matrix
    global_order_parameter = mean(enstrophy, 1, 'omitnan');
    
    % two of the measures we will be using to evaluate fitting
    mean_global(sub) = mean(global_order_parameter, 'omitnan');
    std_global(sub) = std(global_order_parameter, 'omitnan');
    
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
    
    % fitting measure i.e. std of reduced edge matrix
    edge_matrix_std(sub)=std(edge_phases(:));
    
    % if needed, we also calculate the edge spacetime predictability
    tot=length(dffidx100)-1;
    gbindist=edge_phases';
    ncomm=zeros(tot,tot);
    for k=0:7
        for i = 1:tot
            for j = 1:tot
                ccaux=corrcoef(gbindist(1:end-k,i),gbindist(1+k:end,j),'Rows','pairwise');
                ncomm(i,j)=ccaux(2);
                old_ncomm(sub,k+1,i,j)=ccaux(2);
            end
        end
        EdgeSpaTimePredictability(sub,k+1)=mean(mean(ncomm, 'omitnan'), 'omitnan');
    end
end

% for debug only
disp(Rspatime)
disp(mean_global)
disp(std_global)
disp(edge_matrix_std)
disp(EdgeSpaTimePredictability)

% save all the measures in turbu_emp.mat
save turbu_emp.mat Rspatime Rspatime_su mean_global std_global old_ncomm EdgeSpaTimePredictability edge_matrix_std;

% save phases in turbu_emp_Phases.mat
save turbu_emp_Phases.mat Phases;
%% Mean information cascade predictability oevr 8 time shifts

load turbu_emp.mat

figure(1);

plot(EdgeSpaTimePredictability(1,:), 'Linewidth',2); hold on;

title(['Mean Information cascade predictability']);
%% Information cascade predictability for each time shift

figure(2);

old_ncomm = squeeze(mean(old_ncomm, 1, 'omitnan'));

for k=0:7
    for i = 1:tot-1
        ncomm2{i}=old_ncomm(k+1,i,i:end);
    end;
    
    tot=length(dffidx100)-1;
    for i = 1:tot
        tmp=[];
        for j = 1:tot-i
            tmp(j)=ncomm2{j}(i+1);
        end;
        avdiagncomm(i)=mean(tmp, 'omitnan');
    end;
    plot(avdiagncomm,'Linewidth',2); hold on;
    
end;
title(['Information cascade predictability']);
%% statistical significance of turbulence

%p=ranksum(Rspatime,Rspatime_su)

%figure(3)
%boxplot([Rspatime' Rspatime_su'])
%% Phases

figure(4)
plot(Phases(:,500),'r*');
hold on;
plot(Phases_su(:,500),'k*');