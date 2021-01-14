% hopf model simulation for 1003 participants using the Schaefer 1000
% parcellation
%
%  12.01.2021 Samuel Liebana, Morten Kringelbach and Gustavo Deco

% clean up env for clean run
clearvars; close all; clc

%%
load cog_schaefer1000.mat
NPARCELS=1000;


for i=1:NPARCELS
    for j=1:NPARCELS
        rr(i,j)=sqrt((cog(i,1)-cog(j,1))^2+(cog(i,2)-cog(j,2))^2+(cog(i,3)-cog(j,3))^2);
    end;
end;
figure, imagesc(rr);
title('Distance matrix')
axis square;

for k=1:floor(max(max(rr)))
    pairs{k}=find(abs(rr-k)<0.5); 
    elem(k)=size(pairs{k},1);
end;
figure, plot(elem);
title('Distance r of pairs')

% finds the indexes of the 

[u,v] = find(triu(ones(NPARCELS),1));    % get edges
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

NSUB=100; % number of subjects in the empirical data
NR=400;


NRini=20;
NRfin=80;
NSUBSIM=100; % number of subjects to simulate


range=max(max(rr));
delta=range/NR;


% define the coupling matrix (either with Markov rule or empirical)

Markov_SC=1;
if Markov_SC==1
    lambda=0.18;
    for i=1:NPARCELLS
        for j=1:NPARCELLS
            C(i,j)=exp(-lambda*rr(i,j));
        end
        C(i,i)=0;
    end
else
    load sc_schaefer.mat
    C=sc_schaefer;
    C=C/max(max(C));
end

% each parcel is coupled with itself with a magnitude of 1

C1=C;
for i=1:NPARCELLS
    C1(i,i)=1;
end

neighbours=cell(1,NPARCELLS);
for i=1:NPARCELLS
    for j=1:NPARCELLS
        r=rr(i,j);
        index=floor(r/delta)+1;
        if index==NR+1
            index=NR;
        end
        if index>1 && index<=35
            neighbours{i}=[neighbours{i} j];
        end
    end
end

%%
% Parameters of the data
TR=0.72;  % Repetition Time (seconds)

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                    % lowpass frequency of filter (Hz)
fhi = 0.08;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter
Isubdiag = find(tril(ones(NPARCELLS),-1));

%%%%
fce1=zeros(100,NPARCELLS,NPARCELLS);
tss=zeros(NPARCELLS,1200);
PowSpect=zeros(600,NPARCELLS,100);


% skip preprocessing of 1003 participants
load hpcdata1003_f_diff_fce.mat

% Parameters for Hopf model
Tmax=1200;
omega = repmat(2*pi*f_diff',1,2); omega(:,1) = -omega(:,1);
dt=0.1*TR/2;
sig=0.01;
dsig = sqrt(dt)*sig;

%%

linfunc = @(A, x)(A(1)*x+A(2));
options=optimset('MaxFunEvals',10000,'MaxIter',1000,'Display','off');

fcsimul=zeros(NSUBSIM,NPARCELLS,NPARCELLS);
corrfcn=zeros(NSUBSIM,NPARCELLS,NR);
ensspasub=zeros(NSUBSIM,NPARCELLS);
ensspasub1=zeros(NSUBSIM,NPARCELLS);
DTspatime=zeros(NPARCELLS,Tmax);
DTspatime1=zeros(NPARCELLS,Tmax);
Rsub=zeros(1,NSUBSIM);
DTsub=zeros(1,NSUBSIM);

G_range=0.5:0.05:2;
enstrophy_r=zeros(length(G_range),NR);

ii=1;
for G=G_range
    
    G
    wC = G*C;
    sumC = repmat(sum(wC,2),1,2); % for sum Cij*xj
    
    fcsimul=zeros(NSUBSIM,NPARCELLS,NPARCELLS);
    ensspasub=zeros(NSUBSIM, NPARCELLS);
    corrfcn=zeros(NSUBSIM,NPARCELLS,NR);
    %% Hopf Simulation
    for sub=1:NSUBSIM
        a=-0.02*ones(NPARCELLS,2);
        xs=zeros(Tmax,NPARCELLS);
        %number of iterations, 100 willk�hrlich, weil reicht in diesem Fall
        z = 0.1*ones(NPARCELLS,2); % --> x = z(:,1), y = z(:,2)
        nn=0;
        % discard first 2000 time steps
        for t=0:dt:2000
            suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
            zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
            z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(NPARCELLS,2);
        end
        % actual modelling (x=BOLD signal (Interpretation), y some other oscillation)
        for t=0:dt:((Tmax-1)*TR)
            suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
            zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
            z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(NPARCELLS,2);
            if abs(mod(t,TR))<0.01
                nn=nn+1;
                xs(nn,:)=z(:,1)';
            end
        end
        ts=xs';
        signal_filt=zeros(NPARCELLS,Tmax);
        Xanalytic=zeros(NPARCELLS,Tmax);
        Phases=zeros(NPARCELLS,Tmax);
        for seed=1:NPARCELLS
            ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
            signal_filt(seed,1:1200)=filtfilt(bfilt,afilt,ts(seed,:));
            Xanalytic = hilbert(detrend(signal_filt(seed,:), 'constant'));
            Phases(seed,:) = angle(Xanalytic);
        end
        
        %%% Perturbation
        a=(-0.02+0.02*rand(NPARCELLS,1)).*ones(NPARCELLS,2);
        nn=0;
        
        for t=0:dt:((Tmax-1)*TR)
            suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
            zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
            z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(NPARCELLS,2);
            if abs(mod(t,TR))<0.01
                nn=nn+1;
                xs(nn,:)=z(:,1)';
            end
        end
        ts=xs';
        Rspatime1=zeros(NPARCELLS,Tmax);
        Tspatime1=zeros(NPARCELLS,Tmax);
        
        signal_filt2=zeros(NPARCELLS,Tmax);
        Xanalytic2=zeros(NPARCELLS,Tmax);
        Phases2=zeros(NPARCELLS,Tmax);
        for seed=1:NPARCELLS
            ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
            signal_filt2(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
            Xanalytic2 = hilbert(detrend(signal_filt2(seed,:), 'constant'));
            Phases2(seed,:) = angle(Xanalytic2);
        end
        
        parfor i=1:NPARCELLS
            % compute enstrophy
            enstrophy=sum(C1(i,:)'.*complex(cos(Phases2),sin(Phases2)))/sum(C1(i,:));
            Rspatime1(i,:)=abs(enstrophy);
            Tspatime1(i,:)=angle(enstrophy);
        end
        
        DTspatime1=zeros(NPARCELLS,Tmax);
        for i=1:NPARCELLS
            DTspatime1(i,:)=Tspatime1(i,:)-mean(Tspatime1(neighbours{i},:));
        end
        ensspasub1(sub,:)=(mean(Rspatime1,2))';
        
        %%%  end perturbation
        
        fcsimulc=zeros(NPARCELLS,NPARCELLS);
        fcsimulc(:,:)=corrcoef(signal_filt'); %compute current corr of fc
        fcsimul(sub,:,:)=zeros(NPARCELLS,NPARCELLS);
        fcsimul(sub,:,:)=corrcoef(signal_filt'); % save this fc
        
        
        % Calculate enstrophy
        Rspatime=zeros(NPARCELLS,Tmax);
        Tspatime=zeros(NPARCELLS,Tmax);
        corrfcnt=zeros(NPARCELLS,NR);
        parfor i=1:NPARCELLS
            numind=zeros(1,NR);
            corrfcn_1=zeros(1,NR);
            for j=1:NPARCELLS
                r=rr(i,j);
                index=floor(r/delta)+1;
                if index==NR+1
                    index=NR;
                end
                %                mcc=fcsimul(sub,i,j);
                mcc=fcsimulc(i,j); % get the current value
                if ~isnan(mcc)
                    corrfcn_1(index)=corrfcn_1(index)+mcc;
                    numind(index)=numind(index)+1;
                end
            end
            
            corrfcnt(i,:)=corrfcn_1./numind;
            
            %           corrfcn(sub,i,:)=corrfcn_1./numind; % save current corrfcn
            
            %%% enstrophy
            enstrophy=sum(C1(i,:)'.*complex(cos(Phases),sin(Phases)))/sum(C1(i,:));
            Rspatime(i,:)=abs(enstrophy);
            Tspatime(i,:)=angle(enstrophy);
        end
        corrfcn(sub,:,:)=zeros(NPARCELLS,NR);
        corrfcn(sub,:,:)=corrfcnt(:,:);
        
        DTspatime=zeros(NPARCELLS,Tmax);
        for i=1:NPARCELLS
            DTspatime(i,:)=Tspatime(i,:)-mean(Tspatime(neighbours{i},:));
        end
        
        Rsub(sub)=std(Rspatime(:));
        DTsub(sub)=std(DTspatime(:));
        ensspasub(sub,:)=(mean(Rspatime,2))';
        save(['Progress_11_01_2021_G_' num2str(G) '_sub_' num2str(sub)],'G','sub');
    end % subsim

    %% compute measures 

    fcsim=squeeze(mean(fcsimul,1));

    % compute segint, ie mean segregation/integration
    [MM QQ]=community_louvain(fcsim,[],[],'negative_sym');
    modularity=QQ;
    integration=mean(mean(abs(fcsim)-eye(NPARCELLS)));
    segint(ii)=modularity*integration

    % compute mean and std for each simulated matrix
    asegint=zeros(1,NSUBSIM);
    for t=1:NSUBSIM
        currfcsim=squeeze(fcsimul(t,:,:));
        [MM QQ]=community_louvain(currfcsim,[],[],'negative_sym');
        modularity=QQ;
        integration=mean(mean(abs(fcsim)-eye(NPARCELLS)));
        asegint(t)=modularity*integration;
    end;
    asegintmean(ii)=mean(asegint);
    asegintstd(ii)=std(asegint);

    % compute fcfitt
    ccaux=corrcoef(fce(Isubdiag),fcsim(Isubdiag),'rows','pairwise');
    fcfitt(ii)=ccaux(2)
    % compute mean and std for each simulated matrix
    afcfitt=zeros(1,NSUBSIM);
    for t=1:NSUBSIM
        currfcsim=squeeze(fcsimul(t,:,:));
        ccaux=corrcoef(fce(Isubdiag),currfcsim(Isubdiag),'rows','pairwise');
        afcfitt(t)=ccaux(2);
    end;
    afcfittmean(ii)=mean(afcfitt)
    afcfittstd(ii)=std(afcfitt)
    

    % compute fcfitt_hete
    fce2=fce-eye(NPARCELLS);
    GBCemp=mean(fce2,2);
    fcsim2=fcsim-eye(NPARCELLS);
    GBCsim=mean(fcsim2,2);
    for i=1:NPARCELLS
        errgbc1(i)=sqrt((GBCsim(i)-GBCemp(i))^2);
    end
    fcfitt_hete(ii)=mean(errgbc1)
 
    % compute Rmeta, mean and std
    Rmeta(ii)=mean(Rsub)
    Rmetastd(ii)=std(Rsub)

    % compute DTmeta, mean and std
    DTmeta(ii)=mean(DTsub)
    DTmetastd(ii)=std(DTsub)
    
    % compute infocapacity, mean and std
    infocapacity(ii)=mean(std(ensspasub1-ones(NSUBSIM,1)*mean(ensspasub)))
    infocapacitystd(ii)=std(std(ensspasub1-ones(NSUBSIM,1)*mean(ensspasub)))

    % compute susceptibility, mean and std
    susceptibility(ii)=mean(mean(ensspasub1-ones(NSUBSIM,1)*mean(ensspasub)))
    susceptibilitystd(ii)=std(mean(ensspasub1-ones(NSUBSIM,1)*mean(ensspasub)))
    
    corrfcn_he=squeeze(mean(corrfcn));
    grandcorrfcn=squeeze(mean(corrfcn_he));
    %% Comparison
    for i=1:NPARCELLS
        for k=NRini:NR
            err11(k)=sqrt((corrfcn_he(i,k)-empcorrfcn(i,k))^2);
        end
        err1(i)=mean(err11);
    end
    err_hete(ii)=mean(err1)
    
    for k=NRini:NR
        errg1(k)=sqrt((grandcorrfcn(k)-empgrandcorrfcn(k))^2);
    end
    err_grand(ii)=mean(errg1)
    
    %%% Powerlaw a
   
    clear xcoor;
    clear ycoor;
    nn=1;
    for k=NRini:NRfin
        if grandcorrfcn(k)>0
            xcoor(nn)=log(xrange(k));
            ycoor(nn)=log(grandcorrfcn(k)/grandcorrfcn(NRini));
            nn=nn+1;
        end
    end
    A0=[-1 1];
    [Afit Residual]= lsqcurvefit(linfunc,A0,xcoor,ycoor,[-4 -10],[4 10],options);
    GoFpow(ii)=Residual
    SlopeCorrPow(ii)=abs(Afit(1));
    
    %%%
    
    ii=ii+1;
end

%% plot figures

figure(1)
plot(G_range,err_grand);
figure(2)
plot(G_range,err_hete);
figure(3)
plot(G_range,fcfitt);
figure(4)
plot(G_range,fcfitt_hete);
figure(5)
plot(G_range,segint);
figure(6)
plot(G_range,GoFpow);
figure(7)
plot(G_range,SlopeCorrPow);
figure(8)
plot(G_range,infocapacity*5);
figure(9)
plot(G_range,susceptibility);
figure(10)
plot(G_range,Rmeta*10);
figure(11)
plot(G_range,DTmeta);

figure(12)
shadedErrorBar(G_range, asegintmean, asegintstd./sqrt(NSUBSIM))
figure(13)
shadedErrorBar(G_range, afcfittmean, afcfittstd./sqrt(NSUBSIM))

figure(14)
shadedErrorBar(G_range, infocapacity*5, (infocapacitystd./sqrt(NSUBSIM))*5)
figure(15)
shadedErrorBar(G_range, susceptibility, (susceptibilitystd./sqrt(NSUBSIM)))

figure(16)
shadedErrorBar(G_range, Rmeta, Rmeta./sqrt(NSUBSIM))
figure(17)
shadedErrorBar(G_range, DTmeta, DTmeta./sqrt(NSUBSIM))
