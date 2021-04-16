clear all;
load empirical_spacorr_rest.mat;
load cog_schaefer100.mat;

NSUB=1003;
NPARCELLS=100;


%%%
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

load hcp1003_schaefer100_REST1_LR_all.mat

nsub=1;
for sub=1:NSUB
    sub
    
    if sub==101
        
        fce2(1,:,:)=squeeze(mean(fce1));
        Power_Areas2(1,:,:)=squeeze(mean(PowSpect,3));
        nsub=1;

    end
    if sub==201
        
        fce2(2,:,:)=squeeze(mean(fce1));
        Power_Areas2(2,:,:)=squeeze(mean(PowSpect,3));
        nsub=1;
        
    end
    if sub==301
        
        fce2(3,:,:)=squeeze(mean(fce1));
        Power_Areas2(3,:,:)=squeeze(mean(PowSpect,3));
        nsub=1;
        
    end
    if sub==401
        
        fce2(4,:,:)=squeeze(mean(fce1));
        Power_Areas2(4,:,:)=squeeze(mean(PowSpect,3));
        nsub=1;
       
    end
    if sub==501
       
        fce2(5,:,:)=squeeze(mean(fce1));
        Power_Areas2(5,:,:)=squeeze(mean(PowSpect,3));
        nsub=1;
       
    end
    if sub==601
      
        fce2(6,:,:)=squeeze(mean(fce1));
        Power_Areas2(6,:,:)=squeeze(mean(PowSpect,3));
        nsub=1;
        
    end
    if sub==701
        
        fce2(7,:,:)=squeeze(mean(fce1));
        Power_Areas2(7,:,:)=squeeze(mean(PowSpect,3));
        nsub=1;
     
    end
    if sub==801
       
        fce2(8,:,:)=squeeze(mean(fce1));
        Power_Areas2(8,:,:)=squeeze(mean(PowSpect,3));
        nsub=1;
     
    end
    if sub==901
     
        fce2(9,:,:)=squeeze(mean(fce1));
        Power_Areas2(9,:,:)=squeeze(mean(PowSpect,3));
        nsub=1;
        
    end
    
    ts=subject{sub}.schaeferts;
    [Ns, Tmax]=size(ts);
    TT=Tmax;
    Ts = TT*TR;
    freq = (0:TT/2-1)/Ts;
    nfreqs=length(freq);
    
    for seed=1:NPARCELLS
        x=detrend(ts(seed,:)-mean(ts(seed,:)));
        tss(seed,:)=filtfilt(bfilt,afilt,x);
        pw = abs(fft(tss(seed,:)));
        PowSpect(:,seed,nsub) = pw(1:floor(TT/2)).^2/(TT/TR);
    end
    fce1(nsub,:,:)=corrcoef(tss(1:NPARCELLS,:)','rows','pairwise');
    nsub=nsub+1;
end


fce2(10,:,:)=squeeze(mean(fce1));
Power_Areas2(10,:,:)=squeeze(mean(PowSpect,3));

fce=squeeze(mean(fce2));
Power_Areas=squeeze(mean(Power_Areas2));
for seed=1:NPARCELLS
    Power_Areas(:,seed)=gaussfilt(freq,Power_Areas(:,seed)',0.01);
end

[maxpowdata,index]=max(Power_Areas);
f_diff = freq(index);
f_diff(find(f_diff==0))=mean(f_diff(find(f_diff~=0)));

clear fce1;
clear fce2;
clear PowSpect;
clear Power_Areas;
clear Power_Areas2;

save hpcdata1003_f_diff_fce_100.mat f_diff fce;
