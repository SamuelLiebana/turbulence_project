#!/bin/bash
#SBATCH --qos=vip
#SBATCH --job-name=SuperRing
#SBATCH --mail-type=END
#SBATCH --mail-user=
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=2
#SBATCH --array=1-630
#SBATCH --output=Super%A_%a.out
#SBATCH --error=Super%A_%a.err

#Load Matlab 2017a module
ml MATLAB

matlab -nojvm -nodisplay<<-EOF

s=str2num(getenv('SLURM_ARRAY_TASK_ID'))

%% Structure ring

lambda=1;
IntegStepSize=0.1;

NPARCELLS=1000;
Isubdiag = find(triu(ones(NPARCELLS),-1));

for i=1:NPARCELLS
    for j=1:NPARCELLS
        d1=abs(j-i);
        d2=NPARCELLS-d1;
        rr(i,j)=min(d1,d2)*IntegStepSize;
    end
end

%% Ordering distances for Edge Matrix
[u,v] = find(triu(ones(NPARCELLS),1));    % get edges
for i=1:length(u)
    distance(i)=floor(2*rr(u(i),v(i)));
end
[adistance,idx2]=sort(distance);
dffidx100 = find(diff(adistance));

for i=1:length(u)
    distance(i)=floor(rr(u(i),v(i)));
end
[adistance,idx2]=sort(distance);
dffidx50 = find(diff(adistance));

%%

for i=1:NPARCELLS
    for j=1:NPARCELLS
        C(i,j)=0.5*exp(-lambda*rr(i,j));
    end
    C(i,i)=0.5;
end


ni=1;
for i=1:20:NPARCELLS
    nj=1;
    for j=1:20:NPARCELLS
        C50(ni,nj)=0.5*exp(-lambda*rr(i,j));
        nj=nj+1;
    end
    C50(ni,ni)=0.5;
    ni=ni+1;
end

ni=1;
for i=1:10:NPARCELLS
    nj=1;
    for j=1:10:NPARCELLS
        C100(ni,nj)=0.5*exp(-lambda*rr(i,j));
        nj=nj+1;
    end
    C100(ni,ni)=0.5;
    ni=ni+1;
end

%% Parameters models

beta=2.3;
alpha=atan(beta);
Tmax=100;
dt=0.01;

K=0.05;
omega=ones(NPARCELLS,1);
coupling=K*sqrt(1+beta^2);

omega=omega/coupling;
coupling=1;

NTmax=Tmax/dt+2;

phi=zeros(NPARCELLS,NTmax);
phi100=zeros(100,NTmax);
phi50=zeros(50,NTmax);
Rspatime1000=zeros(NPARCELLS,NTmax);
Rspatime100=zeros(100,NTmax);
Rspatime50=zeros(50,NTmax);
edge_phases100=zeros(length(dffidx100)-1,NTmax);
edge_phases50=zeros(length(dffidx50)-1,NTmax);
phase_lock_matrix_red100=zeros(1,length(dffidx100)-1);
phase_lock_matrix_red50=zeros(1,length(dffidx50)-1);

Edge_meta100=zeros(1,100);
Edge_meta50=zeros(1,100);
LocalKoP100=zeros(1,100);
LocalKoP50=zeros(1,100);
LocalKoP1000=zeros(1,100);
KoP1000=zeros(1,100);
StdKoP1000=zeros(1,100);
EdgeSpaTimePredictability100=zeros(100,8);
EdgeSpaTimePredictability50=zeros(100,8);

Sigma_range=0.00001:0.00001:0.0063;
sigma=Sigma_range(s);

D_val=sigma*sqrt(1+beta^2)/K;

for sub=1:100
    sub
    dsig=sqrt(2*D_val*dt); 
    phi(:,1)=0.1*ones(NPARCELLS,1);
    for t=0:dt:500
        resta=repmat(phi(:,1),1,NPARCELLS)-repmat(phi(:,1)',NPARCELLS,1)+alpha;
        suma=IntegStepSize*sum(C.*sin(resta),2);
        phi(:,1)=wrapTo2Pi(phi(:,1)+dt*(omega-coupling*suma)+dsig*randn(NPARCELLS,1));
    end
    nn=1;
    for t=0:dt:Tmax
        resta=repmat(phi(:,nn),1,NPARCELLS)-repmat(phi(:,nn)',NPARCELLS,1)+alpha;
        suma=IntegStepSize*sum(C.*sin(resta),2);
        phi(:,nn+1)=wrapTo2Pi(phi(:,nn)+dt*(omega-coupling*suma)+dsig*randn(NPARCELLS,1));
        nn=nn+1;
    end
    
    nn=1;
    for n=1:10:NPARCELLS   
        phi100(nn,:)=wrapTo2Pi(atan2(nanmean(sin(phi(n:n+9,:))),nanmean(cos(phi(n:n+9,:)))));
        nn=nn+1;
    end

    nn=1;
    for n=1:20:NPARCELLS
        phi50(nn,:)=wrapTo2Pi(atan2(nanmean(sin(phi(n:n+19,:))),nanmean(cos(phi(n:n+19,:)))));
        nn=nn+1;
    end
    
    %%%%
    Rglobal=abs(sum(complex(cos(phi),sin(phi)),1))/NPARCELLS;
    KoP1000(sub)=mean(Rglobal);
    StdKoP1000(sub)=std(Rglobal);
    
    for i=1:NPARCELLS
        enstrophy=nansum(repmat(squeeze(C(i,:))',1,NTmax).*complex(cos(phi),sin(phi)))/nansum(C(i,:));
        Rspatime1000(i,:)=abs(enstrophy);
    end
    LocalKoP1000(sub)=std(Rspatime1000(:)); 
    
    for i=1:100
        enstrophy100=nansum(repmat(squeeze(C100(i,:))',1,NTmax).*complex(cos(phi100),sin(phi100)))/nansum(C100(i,:));
        Rspatime100(i,:)=abs(enstrophy100);
    end
    LocalKoP100(sub)=std(Rspatime100(:));
    
    for i=1:50
        enstrophy50=nansum(repmat(squeeze(C50(i,:))',1,NTmax).*complex(cos(phi50),sin(phi50)))/nansum(C50(i,:));
        Rspatime50(i,:)=abs(enstrophy50);
    end
    LocalKoP50(sub)=std(Rspatime50(:));
   
    %%%% Edge-centric 100
    
    for n=1:NTmax
        cphi=cos(phi(:,n));
        sphi=sin(phi(:,n));
        restacos=repmat(cphi,1,NPARCELLS)-repmat(cphi',NPARCELLS,1);
        restasin=repmat(sphi,1,NPARCELLS)-repmat(sphi',NPARCELLS,1);
        PhLoMa=sqrt((restacos).^2+(restasin).^2);
        phase_lock_matrix=PhLoMa(Isubdiag);
        for i = 1:length(dffidx100)-1
            phase_lock_matrix_red100(i)=nanmean(phase_lock_matrix(dffidx100(i)+1:dffidx100(i+1)));
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
            end
        end
        EdgeSpaTimePredictability100(sub,k+1)=nanmean(nanmean(ncomm));
    end
    
        %%%% Edge-centric 50
    
    for n=1:NTmax
        cphi=cos(phi(:,n));
        sphi=sin(phi(:,n));
        restacos=repmat(cphi,1,NPARCELLS)-repmat(cphi',NPARCELLS,1);
        restasin=repmat(sphi,1,NPARCELLS)-repmat(sphi',NPARCELLS,1);
        PhLoMa=sqrt((restacos).^2+(restasin).^2);
        phase_lock_matrix=PhLoMa(Isubdiag);
        for i = 1:length(dffidx50)-1
            phase_lock_matrix_red50(i)=nanmean(phase_lock_matrix(dffidx50(i)+1:dffidx50(i+1)));
        end
        
        edge_phases50(:,n)=phase_lock_matrix_red50;
    end
    Edge_meta50(sub)=std(edge_phases50(:));
    
    tot=length(dffidx50)-1;
    gbindist=edge_phases50';
    ncomm=zeros(tot,tot);
    for k=0:7
        for i = 1:tot
            for j = 1:tot
                ccaux=corrcoef(gbindist(1:end-k,i),gbindist(1+k:end,j),'Rows','pairwise');
                ncomm(i,j)=ccaux(2);
            end
        end
        EdgeSpaTimePredictability50(sub,k+1)=nanmean(nanmean(ncomm));
    end
    
end

save(sprintf('Wsuper_%03d.mat',s),'EdgeSpaTimePredictability100','EdgeSpaTimePredictability50','Edge_meta100','Edge_meta50','LocalKoP100','LocalKoP50','LocalKoP1000','KoP1000','StdKoP1000','D_val');
EOF


