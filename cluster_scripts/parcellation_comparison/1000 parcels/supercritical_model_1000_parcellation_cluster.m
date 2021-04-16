function supercritical_model_1000_parcellation_cluster(s)

%load('HCP1000_empirical.mat')
load('hpcdata1003_f_diff_fce_1000.mat')
load('sc_schaefer1000.mat')
load('cog_schaefer1000.mat')

%%

NSUB=100;
NPARCELS = 1000;
Tmax=1200;
TR = 0.72;
lambda=0.18; % coefficient of Markov Distance Rule



% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                  % lowpass frequency of filter (Hz)
fhi = 0.08;                   % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt2,afilt2]=butter(k,Wn);   % construct the filter

% useful to extract the lower triangular portion of a matrix
Isubdiag = find(tril(ones(NPARCELS),-1));

% calculate the euclidean distance between parcels and store in rr matrix
rr = zeros(NPARCELS, NPARCELS);
for i=1:NPARCELS
    for j=1:NPARCELS
        rr(i,j)=sqrt((cog(i,1)-cog(j,1))^2+(cog(i,2)-cog(j,2))^2+(cog(i,3)-cog(j,3))^2);
    end
end
%% Loading the coupling matrix

C=sc_schaefer;
C=C/max(max(C));


C1=C;
for i=1:NPARCELS
    C1(i,i)=1;
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
%% the model

omega1 = repmat(2*pi*f_diff',1,2); omega1(:,1) = -omega1(:,1);
dt=0.1*TR/2;
sig=0.0006;
dsig = sqrt(dt)*sig;


G_range= 0:0.05:0.5;
beta_range = 1:0.2:3; % loop over different beta values 1-3 and optimise

beta_range = beta_range(s)

% matrix with amplitudes of local kuramoto order paramter 
% for each parcel over time
enstrophy=zeros(NPARCELS,Tmax);

% temporary variable to store filtered signal
signal_filt=zeros(NPARCELS,Tmax);
% temporary variable to store hilbert transform of filtered signal
Phases=zeros(NPARCELS,Tmax);

% reduced edge time series matrix (phase version)
edge_phases_sim=zeros(length(dffidx100)-1,Tmax);
% intermediate variable used to store the instantaneous difference in phase 
% between hilbert transforms of node timeseries (i.e. edge time series)
phase_lock_matrix_red=zeros(1,length(dffidx100)-1);

Rspatime_sim=zeros(length(G_range), NSUB);

% mean and std of the global kuramoto order parameter
mean_global_sim=zeros(length(G_range), NSUB);
std_global_sim=zeros(length(G_range), NSUB);

% std of the edge matrix
edge_matrix_std_sim=zeros(length(G_range), NSUB);

% mean of the correlation matrix for correlations between edge time series
% with different shifts
EdgeSpaTimePredictability=zeros(length(G_range), NSUB,8);

%ii = 1;

for beta = beta_range

    omega2=beta*ones(NPARCELS,2);
    omega2(:,1) = -omega2(:,1);
    omega=omega1+omega2;

    kk=1;
    for G = G_range %
        
        
        wC = G*C; % zero diagonal C
        sumC = repmat(sum(wC,2),1,2); % for sum Cij*xj
        % Hopf Simulation
        a=1.3*ones(NPARCELS,2);

        for sub=1:NSUB
            G
            sub
            xs=zeros(Tmax,NPARCELS);
            %number of iterations, 100 willk�hrlich, weil reicht in diesem Fall
            %z = 0.1*ones(NPARCELS,2); % --> x = z(:,1), y = z(:,2)
            phase_in = rand(NPARCELS,1)*2*pi-pi;
            z(:,1)= cos(phase_in);
            z(:,2)=sin(phase_in); % --> x = z(:,1), y = z(:,2)

            nn=0;
            % discard first 2000 time steps
            for t=0:dt:2000
                suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
                zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
                z = z + dt*(a.*z + zz.*omega - (z+omega2.*zz).*(z.*z+zz.*zz) + suma) + dsig*randn(NPARCELS,2);
            end
            % actual modeling (x=BOLD signal (Interpretation), y some other oscillation)
            for t=0:dt:((Tmax-1)*TR)
                suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
                zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
                z = z + dt*(a.*z + zz.*omega - (z+omega2.*zz).*(z.*z+zz.*zz) + suma) + dsig*randn(NPARCELS,2);
                if abs(mod(t,TR))<0.01
                    nn=nn+1;
                    xs(nn,:)=z(:,1)';
                end
            end
            ts=xs';

            for seed=1:NPARCELS
                ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
                signal_filt(seed,:) =filtfilt(bfilt2,afilt2,ts(seed,:));
                Xanalytic = hilbert(detrend(signal_filt(seed,:), 'constant'));
                Phases(seed,:) = angle(Xanalytic);
            end

            for i=1:NPARCELS 
                sumphases=sum(repmat(C1(i,:)',1,Tmax).*complex(cos(Phases),sin(Phases)), 'omitnan')/sum(C1(i,:), 'omitnan');
                enstrophy(i,:)=abs(sumphases);
            end

            Rspatime_sim(kk, sub)=std(enstrophy(:));
            
            % calculate the global order parameter as the sum of phases without the C kernel    
            global_order_parameter =abs(nansum(complex(cos(Phases),sin(Phases)),1))/NPARCELS;

            % two of the measures we will be using to evaluate fitting
            mean_global_sim(kk, sub) = mean(global_order_parameter, 'omitnan');
            std_global_sim(kk, sub) = std(global_order_parameter, 'omitnan');

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
                edge_phases_sim(:,n)=phase_lock_matrix_red;
            end

            % fitting measure i.e. std of reduced edge matrix
            edge_matrix_std_sim(kk, sub)=std(edge_phases_sim(:));
            
             % we also calculate the edge spacetime predictability
            tot=length(dffidx100)-1;
            Isubdiagtot = find(triu(ones(tot),-1));
            gbindist=edge_phases_sim';
            ncomm=zeros(tot,tot);
            for k=0:7
               gbindist_k=gbindist(1:end-k,:);
               gbindist_k2=gbindist(1+k:end,:);
               ncomm=corr(gbindist_k,gbindist_k2,'Rows','pairwise');
        %        for i = 1:tot
        %            for j = 1:tot
        %                
        %                ccaux=corrcoef(gbindist(1:end-k,i),gbindist(1+k:end,j),'Rows','pairwise');
        %                ncomm(i,j)=ccaux(2);
        %                old_ncomm(sub,k+1,i,j)=ccaux(2);
        %            end
        %        end
               %EdgeSpaTimePredictability(sub,k+1)=mean(mean(ncomm, 'omitnan'), 'omitnan');
               EdgeSpaTimePredictability(kk, sub,k+1)=mean(ncomm(Isubdiagtot), 'omitnan');
            end

        end

        kk = kk+1;
    end

    %ii = ii+1;
end


%%

Rspatime_av_sim = Rspatime_sim; % mean(Rspatime_sim, 2);
mean_global_av_sim = mean_global_sim; % mean(mean_global_sim, 2);
std_global_av_sim = std_global_sim; % mean(mean_global_sim, 2);
edge_matrix_std_av_sim = edge_matrix_std_sim; % mean(edge_matrix_std_sim, 2);
EdgeSpaTimePredictability_av_sim = mean(EdgeSpaTimePredictability, 3); % mean(EdgeSpaTimePredictability, [2, 3]);

%Rspatime_diff = abs(Rspatime_av_sim - Rspatime_av);
%mean_global_diff = abs(mean_global_av_sim - mean_global_av);
%std_global_diff = abs(std_global_av_sim - std_global_av);edge_matrix_std_diff = abs(edge_matrix_std_av_sim - edge_matrix_std_av);
%EdgeSpaTimePredictability_diff = abs(EdgeSpaTimePredictability_av_sim - EdgeSpaTimePredictability_av);

save(sprintf('matrix_values_supercrit_1000_beta_%03d.mat',s), 'Rspatime_av_sim', 'mean_global_av_sim', 'std_global_av_sim', 'edge_matrix_std_av_sim', 'EdgeSpaTimePredictability_av_sim');
%save(sprintf('abs_errors_supercrit_1000_beta_%03d.mat',s), 'Rspatime_diff', 'mean_global_diff', 'std_global_diff', 'edge_matrix_std_diff', 'EdgeSpaTimePredictability_diff');

end