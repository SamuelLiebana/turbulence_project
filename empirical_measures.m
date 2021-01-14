% Calculating the Mean of the Global Kuramoto Order Parameter, Std of the
% Global Kuramoto Order Parameter, Std of the Local Kuramoto Order
% Parameter and the Std of the Reduced Edge Matrix
% 
%
%  12.01.2021 Samuel Liebana, Morten Kringelbach and Gustavo Deco

% clean up env for clean run
clearvars; close all; clc

%%

% Parameters of the data
NSUB=1003; 
NPARCELS=100;
Tmax=1200;
lambda=0.18;
TR=0.72;  % Repetition Time (seconds)

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                  % lowpass frequency of filter (Hz)
fhi = 0.08;                   % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter

% useful to extract the lower triangular portion 
% of a matrix
Isubdiag = find(tril(ones(NPARCELS),-1));


% load the centres of gravity
%load cog_schaefer1000.mat;

load(strcat('cog_schaefer', num2str(NPARCELS),'.mat'))

rr = zeros(NPARCELS, NPARCELS);

for i=1:NPARCELS
    for j=1:NPARCELS
        rr(i,j)=sqrt((cog(i,1)-cog(j,1))^2+(cog(i,2)-cog(j,2))^2+(cog(i,3)-cog(j,3))^2);
    end
end

% finds the indexes of the edges with different distances

[u,v] = find(triu(ones(NPARCELS),1));    % get edges
distance = zeros(length(u), 1);

for i=1:length(u)
    % lists the lengths of all the edges
    distance(i)=floor(rr(u(i),v(i)));
end
% sorts the distance array and returns the set of indexes of distance 
% so that we can recreate the sort
[adistance,idx2]=sort(distance);

% finds the indexes (-1) of all the *different* length edges in the network
dffidx100 = find(diff(adistance));

%%

enstrophy=zeros(NPARCELS,Tmax);
enstrophy_su=zeros(NPARCELS,Tmax);
signal_filt=zeros(NPARCELS,Tmax);
Phases=zeros(NPARCELS,Tmax);
Phases_su=zeros(NPARCELS,Tmax);
Rspatime=zeros(1,NSUB);
Rspatime_su=zeros(1,NSUB);
mean_global=zeros(1,NSUB);
std_global=zeros(1,NSUB);
edge_phases100=zeros(length(dffidx100)-1,Tmax);
phase_lock_matrix_red100=zeros(1,length(dffidx100)-1);
Edge_meta100=zeros(1,NSUB);
EdgeSpaTimePredictability100=zeros(NSUB,8);

% load the coupling matrix
%load sc_schaefer.mat
%C=sc_schaefer;
%C=C/max(max(C));

C=zeros(NPARCELS,NPARCELS);
for i=1:NPARCELS
    for j=1:NPARCELS
        C(i,j)=exp(-lambda*rr(i,j));
    end
    C(i,i)=1;
end

%%

%myDir = strcat(pwd, '\schaefer1000_100participants');	% gets directory
%myFiles = dir(fullfile(myDir,'*.mat'));	

load hcp1003_schaefer100_REST1_LR_all.mat

for sub = 1:NSUB
    sub
	%baseFileName = myFiles(sub).name;
	%fullFileName = fullfile(myDir, baseFileName);  % Changed myFolder to myDir
	
    %clear schaeferts
    %load(fullFileName);
    
    %ts=schaeferts;
    
    ts=subject{sub}.schaeferts;
    
    %removes all NaN rows from the matrix
    %ts = ts(all(~isnan(ts),2),:);
    
    % sets all NaNs to 0
    %ts(isnan(ts))=0;

    for seed=1:NPARCELS
        ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
        signal_filt(seed,:)=filtfilt(bfilt,afilt,ts(seed,:));
        Xanalytic = hilbert(detrend(signal_filt(seed,:), 'constant'));
        Phases(seed,:) = angle(Xanalytic);
        Phases_su(seed,:)=Phases(seed,randperm(Tmax));
    end
    
    for i=1:NPARCELS
        sumphases=sum(repmat(C(i,:)',1,Tmax).*complex(cos(Phases),sin(Phases)), 'omitnan')/sum(C(i,:), 'omitnan');
        enstrophy(i,:)=abs(sumphases);
        sumphases=sum(repmat(C(i,:)',1,Tmax).*complex(cos(Phases_su),sin(Phases_su)), 'omitnan')/sum(C(i,:), 'omitnan');
        enstrophy_su(i,:)=abs(sumphases);
    end
    
    global_order_parameter = mean(enstrophy, 1, 'omitnan');
    
    mean_global(sub) = mean(global_order_parameter, 'omitnan');
    std_global(sub) = std(global_order_parameter, 'omitnan');
    
    Rspatime(sub)=std(enstrophy(:));
    Rspatime_su(sub)=std(enstrophy_su(:));
    
    
    
    for n=1:Tmax
        cphi=cos(Phases(:,n));
        sphi=sin(Phases(:,n));
        restacos=repmat(cphi,1,NPARCELS)-repmat(cphi',NPARCELS,1);
        restasin=repmat(sphi,1,NPARCELS)-repmat(sphi',NPARCELS,1);
        PhLoMa=sqrt((restacos).^2+(restasin).^2);
        phase_lock_matrix=PhLoMa(Isubdiag);
        phase_lock_matrix = phase_lock_matrix(idx2); %%%%%%%%%%%%%%%%%%%%
        for i = 1:length(dffidx100)-1
            phase_lock_matrix_red100(i)=mean(phase_lock_matrix(dffidx100(i)+1:dffidx100(i+1)), 'omitnan');
        end
        
        edge_phases100(:,n)=phase_lock_matrix_red100;
    end
    Edge_meta100(sub)=std(edge_phases100(:));
    
    tot=length(dffidx100)-1;
    gbindist=edge_phases100';
    ncomm=zeros(tot,tot);
    for k=0:7
        for i = 1:tot
            for j = 1:tot
                ccaux=corrcoef(gbindist(1:end-k,i),gbindist(1+k:end,j),'Rows','pairwise');
                ncomm(i,j)=ccaux(2);
                old_ncomm(sub,k+1,i,j)=ccaux(2);
            end
        end
        EdgeSpaTimePredictability100(sub,k+1)=mean(mean(ncomm, 'omitnan'), 'omitnan');
    end
    
    
   
end


disp(Rspatime)
disp(mean_global)
disp(std_global)
disp(Edge_meta100)
disp(EdgeSpaTimePredictability100)



save turbu_emp.mat Rspatime Rspatime_su mean_global std_global old_ncomm EdgeSpaTimePredictability100 Edge_meta100;


save turbu_emp_Phases.mat Phases;

%%

load turbu_emp.mat

figure(1);

plot(EdgeSpaTimePredictability100(1,:), 'Linewidth',2); hold on;


title(['Information cascade predictability']);

%%

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

%%

%load turbu_emp_Phases.mat;
%load turbu_emp.mat;


Turbu=mean(Rspatime, 'omitnan')


%p=ranksum(Rspatime,Rspatime_su)

%figure(3)
%boxplot([Rspatime' Rspatime_su'])



%% Phases
figure(4)
plot(Phases(:,500),'r*');
hold on;
plot(Phases_su(:,500),'k*');

%% check:
for i=400:2:600
    i
    plot(Phases(:,i),'r*');
    pause;
end
