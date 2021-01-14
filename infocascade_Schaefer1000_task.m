% information cascade predictability and over time
%
%  01.12.2020 Morten Kringelbach and Gustavo Deco 
%

% clean up env for clean run
clearvars; close all; clc

%% make edge time series

% load example time series
load lang/100307.mat
ts = double(schaeferts');

%% create distance matrix from centre of gravity for schaefer1000
load cog_schaefer1000.mat
for i=1:1000
    for j=1:1000
        dist(i,j)=sqrt((cog(i,1)-cog(j,1))^2+(cog(i,2)-cog(j,2))^2+(cog(i,3)-cog(j,3))^2);
    end;
end;
figure, imagesc(dist);
title('Distance matrix')
axis square;

for k=1:floor(max(max(dist)))
    pairs{k}=find(abs(dist-k)<0.5); 
    elem(k)=size(pairs{k},1);
end;
figure, plot(elem);
title('Distance r of pairs')

%% create edge time series from regional time series
[T,N] = size(ts);
M = N*(N - 1)/2;

% timeseries may contain nan so define own zscore function
zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));
z=zscor_xnan(ts);
[u,v] = find(triu(ones(N),1));    % get edges
ets = z(:,u).*z(:,v);             % edge ts products

%% visualize edge time series

figure, imagesc(ets',[-4,4]); 
title('Edge time series')

gets=ets';
ws=20;
ws_ets=zeros(size(ts,1)-ws,1);
for k=1:size(ts,1)-ws-1
    aux=0;
    for l=k:k+ws
        aux=aux+nansum(gets(:,l).^2);
    end;
    ws_ets(k)=aux;
end;
figure, plot(ws_ets)

%% visualise ets separating the distances with white line

% compute distances in same order as ets
for i=1:length(u)
    distance(i)=floor(dist(u(i),v(i)));
end;

[a,idx2]=sort(distance);
dffidx2 = find(diff(distance(idx2)));
figure, imagesc(ets(:,idx2)',[-4 4]); hold on;
for i = 1:length(dffidx2)
    plot([0.5,T + 0.5],dffidx2(i)*ones(1,2),'w');
end
title('Edge time series sorted by distance')

%% visualise ts and matrices of mean of ets for each distance
zrets=zscore(ets(:,idx2),0,2);
rets=ets(:,idx2);
dffidx3 = find(diff(a));
bindist=zeros(length(dffidx3),T);
binstddist=zeros(length(dffidx3),T);
binmaxdist=zeros(length(dffidx3),T);
for i = 1:length(dffidx3)-1
   bindist(i,:)=nanmean(abs(rets(:, dffidx3(i)+1:dffidx3(i+1))')); 
   binstddist(i,:)=std(rets(:, dffidx3(i)+1:dffidx3(i+1))'); 
   binmaxdist(i,:)=max(rets(:, dffidx3(i)+1:dffidx3(i+1))'); 
end
figure, imagesc(bindist,[-4 4]); hold on;
title('Edge time series mean abs per distance pairs')
figure, imagesc(binmaxdist,[-4 4]); hold on;
title('Edge time series max per distance pairs')
m2=corrcoef(bindist','Rows','pairwise');
figure, imagesc(m2);
axis square;
title('Matrix absmean per distance pairs')
ms2=corrcoef(binstddist','Rows','pairwise');
figure, imagesc(ms2);
axis square;
title('Matrix stdev per distance pairs')
mm2=corrcoef(binmaxdist','Rows','pairwise');
figure, imagesc(mm2);
axis square;
title('Matrix max per distance pairs')

%% visualise info cascade over mean abs 
for i = 1:length(dffidx3)-1
   avm2{i}=m2(i,i:end);
end;

tot=length(dffidx3);
for i = 1:tot
    tmp=[];
    for j = 1:tot-i
        tmp(j)=avm2{j}(i+1);
    end;
	avdiagm(i)=nanmean(tmp);
end;
figure, plot(avdiagm);
title('Information cascade of mean abs ets')

%% visualise information cascade over time
zbindist=zscore(bindist);
tot=length(dffidx3);
spti=zeros(tot,size(ts,1));
clear nm2 avdiagnm;

for k=1:size(ts,1)-ws-1
%    nm2=bindist(1:tot,k).*bindist(1:tot,k)';
    nm2=corrcoef(bindist(1:tot,k:k+ws)','Rows','pairwise');

    for i = 1:tot-1
       avnm2{i}=nm2(i,i:end);
    end;

    for i = 1:tot
        tmp=[];
        for j = 1:tot-i
            tmp(j)=avnm2{j}(i+1);
        end;
        avdiagnm(i)=nanmean(tmp);
    end;
    
    spti(:,k)=avdiagnm;
end;    
figure, imagesc(spti)
title('Information cascade over time')
load cmap
colormap(cmap);
%% 
ws_spti=zeros(size(ts,1)-ws,1);
for k=1:size(ts,1)-ws-1
    aux=0;
    for l=k:k+ws
        aux=aux+nansum(spti(:,l).^2);
    end;
    ws_spti(k)=aux;
end;
figure, plot(ws_spti)

%% visualise predictability matrices

gbindist=bindist';
ncomm=zeros(8,tot,tot);
figure;
for k=0:7
    for i = 1:tot
        for j = 1:tot
            ccaux=corrcoef(gbindist(1:end-k,i),gbindist(1+k:end,j),'Rows','pairwise');
            ncomm(k+1,i,j)=ccaux(2);
        end;
    end;
    tmp=squeeze(ncomm(k+1,:,:));
    subplot(1,8,k+1); imagesc(tmp);
    axis square;
    title(['Step: ' num2str(k)])
end;
sgtitle(['Predictability matrices'])

%% visualise average info cascade predictability 

figure;
for k=0:7
    for i = 1:tot-1
        ncomm2{i}=ncomm(k+1,i,i:end);
    end;
    
    tot=length(dffidx3);
    for i = 1:tot
        tmp=[];
        for j = 1:tot-i
            tmp(j)=ncomm2{j}(i+1);
        end;
        avdiagncomm(i)=nanmean(tmp);
    end;
    plot(avdiagncomm,'Linewidth',2); hold on;
    
end;
title(['Information cascade predictability']);
