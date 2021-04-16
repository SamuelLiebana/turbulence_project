clear all;
load cog_dkt68.mat;

NSUB=89;
NPARCELLS=68;
TR = 1/200;

%%%
Isubdiag = find(tril(ones(NPARCELLS),-1));
Tfinal = 44196;
%%%%

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 8;                    % lowpass frequency of filter (Hz)
fhi = 13;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter


tss=zeros(NPARCELLS,Tfinal);
PowSpect=zeros(floor(Tfinal/2),NPARCELLS, NSUB);

% if the data is stored in individual files in a directory
myDir = strcat(pwd, '\data');	% gets directory
myFiles = dir(fullfile(myDir,'*.mat'));	

for sub=1:NSUB
    sub
	baseFileName = myFiles(sub).name;
	fullFileName = fullfile(myDir, baseFileName);
	
    %clear schaeferts
    load(fullFileName);
    
    ts=X';
    ts = double(ts(:, 1:Tfinal));
    
    [Ns, Tmax]=size(ts);
    TT=Tmax;
    Ts = TT*TR;
    freq = (0:TT/2-1)/Ts;
    nfreqs=length(freq);
    
    for seed=1:NPARCELLS
        x=detrend(ts(seed,:)-mean(ts(seed,:)));
        tss(seed,:)=filtfilt(bfilt,afilt,x);
        pw = abs(fft(tss(seed,:)));
        PowSpect(:,seed,sub) = pw(1:floor(TT/2)).^2/(TT/TR);
    end
    
end


Power_Areas=squeeze(mean(PowSpect,3));

for seed=1:NPARCELLS
    seed
    Power_Areas(:,seed)=gaussfilt(freq,Power_Areas(:,seed)',0.01);
end

[maxpowdata,index]=max(Power_Areas);
f_diff = freq(index);
f_diff(find(f_diff==0))=mean(f_diff(find(f_diff~=0)));

clear PowSpect;
clear Power_Areas;

save MEGdata89_f_diff_fce_68.mat f_diff;
