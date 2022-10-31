%% this  scritps creates most of the figures shown in the Paper, only the nucleosome based 
% analysis in Figure 3D, S2C, S3D, S4F and S5A are done with the FigureNucs script
% the preprocessing scripts can be skipped if the provided data files
% "data4Figures500bp.mat" or "data4Figures2000bp.mat" is used 

%% preprocessing only necessary if the outfiles are imported - after import with getNaamaData data contains the genome coverage from each sample (12mio bases * number of samples) 
% and smeta contians the meta data of each sample the meta data is updaded
% and processed and reordered 
clear all
goodExps={'s39.mat';'s43.mat';'s44.mat';'s47.mat';'x35.mat';'x36.mat';'s51.mat';'t51.mat';'c51.mat'}
[data,smeta]=getNaamaData('figure',26);
smeta.exp(cellfun('prodofsize',smeta.folder)==54)=strrep(smeta.exp(cellfun('prodofsize',smeta.folder)==54),'x','s');
smeta.exp(ismember(smeta.exp,'s51.mat')&contains(smeta.tag,'h3x2b'))={'t51.mat'};
smeta.gt(ismember(smeta.gt,'rtt109'))={'rtt'}

keep=contains(smeta.exp,goodExps);
smeta=smeta(keep,:);
data=data(:,keep);

% here we bin the data basedinto 500 basepair bins and collect addtional information on each bin
[binData,binMeta,~]=DNA_norm_bin(data,'bin_size',500,'calc_ribo_frac',true,'indChr',false,'format','abs','clean',false,'n_ori',358,'n_chr',17,'Tn5',false);
clear data

clearvars keep
smeta.tp=str2double(regexp(smeta.con,'(?<=HU)\d+','match','once'));
smeta.tp(contains(smeta.con,'async'))=-10;
smeta.tp(contains(smeta.con,'alpha'))=-5;

% here we generate mean HA and DNA profiles form all time courses - they
% are highlighted bby c51 in the smeta table
c=0;
clear mSmeta meanBin
for ab={'ha','DNA'}
    for gt=unique(smeta.gt(ismember(smeta.ab,ab)&smeta.bad==0&ismember(smeta.gt,[{'mecsml'}    {'rtt'}    {'vps75'}    {'wt'}])))'
        for tp=unique([smeta.tp(ismember(smeta.gt,gt)&ismember(smeta.ab,ab)&smeta.bad==0)])'
            c=c+1;
            selSmp=find(smeta.tp==tp&ismember(smeta.gt,gt)&ismember(smeta.ab,ab)&smeta.bad==0)
            if numel(selSmp)>0
                mSmeta(c,:)=smeta(selSmp(1),:);
                meanBin(:,c)=mean(binData(:,selSmp),2);
            else
                mSmeta(c,:)=smeta(find(smeta.tp==30&ismember(smeta.gt,gt)&ismember(smeta.ab,ab)&smeta.bad==0,1),:);
                meanBin(:,c)=nan(size(binData,1),1);%mean(binData(:,selSmp),2);
                mSmeta.bad(c)=1;
                mSmeta.tp=tp;
            end
        end
    end
end
 mSmeta.exp(:)={'c51.mat'}
binData=[binData,meanBin];
smeta=[smeta;mSmeta];
clearvars -except binData binMeta smeta goodExps


[tcs,~,smeta.tid]=unique(smeta(:,{'gt','exp','ab'}),'stable')
[smeta,idx]=sortrows(smeta,{'tid','tp'});
binData=binData(:,idx);
clearvars idx
GP=load('/home/labs/barkailab/felixj/Documents/MATLAB/scripts/gene/group_imp.mat');
binMeta=array2table(binMeta,'VariableNames',{'pos','tor','oriDis','oriIdx','chr','bad'});
%% start here if you use the provided processed data file : "data4Figures500bp.mat"
% here we create the aFactor normalizd profiles
binResid=nan(size(binData));
for t=find(contains(tcs.ab,{'9ac','k56ac','myc','ha','DNA'})&contains(tcs.exp,[goodExps]))'
    selSmp=find(smeta.tid==t);
    %alphaSmp=find(smeta.tid==t&smeta.tp==-5)
    if ismember(tcs.ab(t),{'ha','DNA'})
        alphaSmp=find(ismember(smeta.ab,tcs.ab(t))&ismember(smeta.exp,setdiff(goodExps,'c51.mat'))&smeta.tp==-5);
    else
        alphaSmp=find(ismember(smeta.ab,tcs.ab(t))&ismember(smeta.gt,tcs.gt(t))&ismember(smeta.exp,goodExps)&smeta.tp==-5);
    end
    
    for s=selSmp'
        [p,r]=robustfit(log2(mean(binData(:,alphaSmp),2)+0.1),log2(0.1+mean(binData(:,s),2)));
        binResid(:,s)=r.resid;
    end
end
clearvars -except binData binMeta binResid goodExps GP smeta tcs
% compare DNA timecourses
%selSmp=find(ismember(smeta.exp,{'s51.mat','t51.mat'})&smeta.bad==0);
selSmp=find(ismember(smeta.exp,setdiff(goodExps,'c51.mat'))&ismember(smeta.ab,{'ha','DNA'})&ismember(smeta.gt,'wt')&smeta.bad==0&ismember(smeta.tp,[15,30,40,50,70,90,110,130,150,180]))
[~,idx]=sortrows(smeta(selSmp,:),{'gtid','exp','tp'})
selSmp=selSmp(idx);

sData=binData(:,selSmp);
sData(binMeta.bad>0,:)=nan;
sData=smoothdata(sData,'rlowess',7);
selBins=all(sData>0.3,2)

hold off
figure
imagesc(corr(log2(clusterData(:,selSmp))-median(log2(clusterData(:,selSmp)),2,'omitnan'),'rows','pairwise'))

% imagesc(corr(log2(sData(selBins,:))-median(log2(sData(selBins,:)),2,'omitnan'),'rows','pairwise'),[-.5 .7])
% imagesc(corr(binResid(:,selSmp)-median(binResid(:,selSmp),2),'rows','pairwise'))

pMat=nan(numel(selSmp),numel(selSmp),2);
for i=1:numel(selSmp)
    for j=1:numel(selSmp)
        pMat(i,j,:)=linortfit2(clusterRes(clusterSize>50,selSmp(i)),clusterRes(clusterSize>50,selSmp(j)));
    end
end
bL=find(diff(smeta.tid(selSmp)))
figure
%imagesc(corr(log2(clusterData(:,selSmp))-median(log2(clusterData(:,selSmp)),2,'omitnan'),'rows','pairwise'))
%cMap=[brewermap(128,'purples');flipud(brewermap(128,'purples'))]
imagesc(log2(pMat(:,:,1)),[-.5 .5])
xT=movmean([0;bL;numel(selSmp)],2,'Endpoints','discard');
hold on
plot(xlim',repmat(bL',2,1)+0.5,'k-','Linewidth',2)
plot(repmat(bL',2,1)+0.5,ylim','k-','Linewidth',2)
xticks(xT+0.5)
xticklabels(smeta.exp(selSmp(xT)))
yticks(movmean([0;bL;numel(selSmp)],2,'Endpoints','discard')+0.5)
yticklabels(smeta.exp(selSmp(xT)))
colormap(cMap)
%ylabel(colorbar(),'correlation along clusters')
%title('correlation comparisson between DNA and HA')
title('slope between the clusters')
ylabel(colorbar(),'slope between the clusters-log2')

[val,idx]=sort(binMeta.tor.*(binMeta.bad==0));
idx=idx(val>0)
figure
temp=movmean(binData,11);
imagesc(corr(temp(idx,smeta.tid==47)'))

%% Genome Raw - THe raw data profile as for example shown in  Fig 1A)
%distOris=[1;diff(oriBin)>40]&[diff(oriBin)>40;1];
%figure;plot(diff(oriBin(distOris)))
surLen=60
oriBin=find(abs(binMeta{:,'oriDis'})<250&[1:size(binMeta,1)]'>surLen&[1:size(binMeta,1)]'<(size(binMeta,1)-surLen) & ~ismember(binMeta.oriIdx,find(contains(GP.oris.Var4,'1216'))))

midBin=oriBin(213)%cumsum(distOris)==64)
xLim=[6339301     6677119];
oriColor=gray(3);
oriLim=[21,25];
xTicks=[6300000     6350000     6400000     6450000     6500000     6550000     6600000     6650000     6700000     6750000];
showRpts=table({'wt','rtt','vps75','mec'}',{{'s47','s51'},{'x35','s51'},{'s44','s51'},{'x35'}}')
close all
for r=1:size(showRpts,1)
    figure('Position',[1 41 1920 962])
    c=0;
    for t=find(contains(tcs.exp,showRpts.Var2{r})&contains(tcs.gt,showRpts.Var1{r}))'%find(contains(tcs.exp,{'39','t51'})&contains(tcs.gt,'rtt'))'
        c=c+1;
        subplot(13,1,c*2+[-1 0])
        hold off
        selSmp=find(smeta.tid==t&contains(smeta.exp,tcs.exp{t})&smeta.tp>=15&smeta.bad==0);
        imagesc(smoothdata(smoothdata(binData(midBin+[-500:500],selSmp),'gaussian',7)','movmean',1),'XData',quantile(binMeta.pos(midBin+[-500,500]),[0 1]))
        title(strjoin(tcs{t,1:3}))
        caxis([1 max(caxis)*0.8])
        if c==1;
            selOri=oriBin(binMeta.tor(oriBin)<23)
            hold on
            scatter(binMeta.pos(selOri),repmat(0.5,numel(selOri),1),[],oriColor(sum(binMeta.tor(selOri)<oriLim,2)+1,:),'filled')
        end
        xticks(xTicks)
        xticklabels(round((xTicks-GP.chrIdx(binMeta.chr(oriBin(213))))/1000))
        colorbar()
        colormap(gca,brewermap(128,'Blues'))
        xlim(xLim)
    end
    t=find(contains(tcs.exp,showRpts.Var2{r})&contains(tcs.gt,showRpts.Var1{r})&contains(tcs.ab,'k9'))'
    if numel(t)==1
        
        c=c+1;
        subplot(13,1,c*2+[-1 0])
        selSmp=find(smeta.tid==t&contains(smeta.exp,tcs.exp{t})&smeta.tp>=15);
        imagesc(smoothdata(smoothdata(binResid(midBin+[-500:500],selSmp),'gaussian',1)','movmean',1),'XData',quantile(binMeta.pos(midBin+[-500,500]),[0 1]),[0 .65])
        title(sprintf('norm %s',strjoin(tcs{t,1:3})))
        colormap(gca,brewermap(128,'Purples'))
        colorbar()
        xlim(xLim)
                xticks(xTicks)
        xticklabels(round((xTicks-GP.chrIdx(binMeta.chr(oriBin(213))))/1000))
    end
    subplot(26,1,26)
    imagesc(smoothdata(smoothdata(binMeta.tor(midBin+[-500:500]),'gaussian',1)','movmean',1),'XData',quantile(binMeta.pos(midBin+[-500,500]),[0 1]))
    colormap(gca,flipud(brewermap(128,'BuPu')))
    colorbar()
    title('RT')
    xlim(xLim)    
                xticks(xTicks)
        xticklabels(round((xTicks-GP.chrIdx(binMeta.chr(oriBin(213))))/1000))
        save_gf(gcf,sprintf('GenRaw_%s',showRpts.Var1{r}),'type',{'fig'},'paper','gilad22','size',[])

end

%% prepare the oricentric analysis as shown in  Fig 1
surLen=60
oriBin=find(abs(binMeta{:,'oriDis'})<250&[1:size(binMeta,1)]'>surLen&[1:size(binMeta,1)]'<(size(binMeta,1)-surLen) & ~ismember(binMeta.oriIdx,find(contains(GP.oris.Var4,'1216'))))
clear tempMat tempRes
for o=1:numel(oriBin)
    tempMat(o,:,:)=binData(round(oriBin(o)+[-surLen:surLen]),:);
    tempRes(o,:,:)=binResid(round(oriBin(o)+[-surLen:surLen]),:);
end
oriCluster=min(sum(binMeta.tor(oriBin)>=linspace(min(binMeta.tor(oriBin)),max(binMeta.tor(oriBin)),8),2),5);
oriCluster(contains(GP.oris.Var4(binMeta.oriIdx(oriBin)),'1216'))=0;
clear oriRes oriMean
for c=1:max(oriCluster)
    selO=oriCluster==c;
    oriMean(c,:,:)=.5*mean(tempMat(oriCluster==c,surLen+1:-1:1,:),1)+.5*mean(tempMat(oriCluster==c,surLen+1:end,:),1);
    oriRes(c,:,:)=.5*mean(tempRes(oriCluster==c,surLen+1:-1:1,:),1)+.5*mean(tempRes(oriCluster==c,surLen+1:end,:),1);
end

parSmp=find(contains(smeta.ab,{'ha','DNA'}))
clearvars forkPointR L5fitR forkPoint L5fit
parfor p=1:numel(parSmp)'%15:226%find(contains(smeta.ab,'ha')&smeta.tp>10)'
    s=parSmp(p)
    for o=1:5%size(oriRes,1)
        temp=L5P(1:surLen+1,smoothdata(oriRes(o,:,s),2,'movmean',5));
        [~,fP]=min(differentiate(temp,1:61));
        forkPointR(o,p)=fP;
        L5fitR{o,p}=temp
    end
end
forkPoint(:,parSmp)=forkPointR;
L5fit(:,parSmp)=L5fitR;
clearvars forkPointR L5fitR


%% Ori figures - e.g. 1C
allExps=unique(tcs(:,{'gt','exp'}))
close all
c=0
conSmp=627;
intTps=unique(smeta.tp(ismember(smeta.ab,'DNA')))';

for r=1:size(showRpts,1)
    figure('Units','normalized','Position',[0.5698 0.2324 0.3344 0.5019], 'color','w');    
    c=0
    for t=find(contains(tcs.exp,showRpts.Var2{r})&contains(tcs.gt,showRpts.Var1{r}))'%find(contains(tcs.exp,{'39','t51'})&contains(tcs.gt,'rtt'))'
        c=c+1
        subplot(2,3,c)
        hold off
        selSmp=find(smeta.tid==t&contains(smeta.exp,tcs.exp{t})&smeta.tp>=15&ismember(smeta.tp,intTps)&smeta.bad==0);
        %imageMat=reshape(permute(smoothdata(oriRes(:,:,selSmp),2,'rloess',1),[3,2,1]),numel(selSmp),[]);
        imageMat=reshape(permute(oriRes(:,:,selSmp),[3,2,1]),numel(selSmp),[]);
        imagesc(imageMat)
        %[~,selHa]=ismember(regexprep(smeta.name(selSmp),{'*','^.*?-'},{'','ha-'}),regexprep(smeta.name,'*',''));
        [~,selHa]=ismember([smeta(selSmp,{'gt','tp'}),smeta(ones(numel(selSmp),1)*conSmp,{'exp','ab'})],smeta(:,{'gt','tp','exp','ab'}))
        hold on
        for o=1:5
            scatter(1+forkPoint(o,selHa(selHa>0))+(o-1).*(surLen+1),find(selHa>0),'r.')
        end
        ylabel('time')
        set(gca,'ytick',1:numel(selSmp),'YTickLabel',smeta.tp(selSmp),'xtick',[])
        xlabel('oriCluster')
        hold on
        plot([surLen+1.5:surLen+1:size(imageMat,2)].*[1;1],repmat(ylim',1,max(oriCluster)-1),'k-')
        colorbar()
        title(strjoin(tcs{t,[1,3,2]}))
        caxis([-0.1 max(caxis).*0.8])
        colormap(gca,brewermap(128,'BuPu'))
        xlim([.5 4*(surLen+1)])
    end
    sgtitle('log2 change')
    save_gf(gcf,sprintf('FigOris_%s',showRpts.Var1{r}),'type',{'fig'},'paper','gilad22','size',[])
end

%% prepare spatial correlation for Figure 1D
conSmp=627; % DNA combined
absCorr=nan(5,size(smeta,1));
diffCorr=nan(5,size(smeta,1));
xCorrAbsDiff=nan(5,size(smeta,1));
for t=find(ismember(tcs.exp,goodExps))'
    selSmp=find(smeta.tid==t&contains(smeta.exp,tcs.exp{t})&smeta.tp>=15&ismember(smeta.tp,intTps));
    [~,selHa]=ismember([smeta(selSmp,{'gt','tp'}),smeta(ones(numel(selSmp),1)*conSmp,{'exp','ab'})],smeta(:,{'gt','tp','exp','ab'}));
    for o=1:5
        for i=find(selHa>0)'
            absCorr(o,selSmp(i))=corr(oriRes(o,:,selSmp(i))',L5fit{o,selHa(i)}([1:1:61]));
            diffCorr(o,selSmp(i))=corr(oriRes(o,:,selSmp(i))',differentiate(L5fit{o,selHa(i)},[1:1:61]));
            xCorrAbsDiff(o,selSmp(i))=corr(L5fit{o,selHa(i)}([1:1:61]),differentiate(L5fit{o,selHa(i)},[1:1:61]));
        end
    end
end


%% plot one genotype - figure 1D
[allAbs,~,tcs.abId]=unique(tcs.ab)
types={'','greys','greens','blues','Browns'}
figure
o=1
intTps=[15,30,40,70,90,110,130,150,180]
for t=find(contains(tcs.gt,'wt')&ismember(tcs.exp,goodExps)&~contains(tcs.exp,'51'))'
    selSmp=find(smeta.tid==t&contains(smeta.exp,tcs.exp{t})&ismember(smeta.tp,intTps));
    tpCols=brewermap(numel(selSmp)+5,types{tcs.abId(t)});
    if true%contains(tcs.exp(t),'43')
        scatter(absCorr(o,selSmp),-diffCorr(o,selSmp),rescale(xCorrAbsDiff(o,selSmp),70,240,'InputMax',0,'InputMin',-1), tpCols(1:end-5,:),'filled','DisplayName',tcs.exp{t})
    else
        scatter(absCorr(o,selSmp),-diffCorr(o,selSmp),rescale(xCorrAbsDiff(o,selSmp),70,240,'InputMax',0,'InputMin',-1), tpCols(1:end-5,:),'DisplayName',tcs.exp{t})
    end
    hold on
end
plot(xlim,xlim,'k--')
axis tight
xlabel('correlation with "DNA"')
ylabel('correlation with d/dx DNA')
scatter([0:0.2:.6],[0 0 0 0 ],rescale([-1:1/3:0],70,240,'InputMax',0,'InputMin',-1),[0 0 0],'Linewidth',2)
save_gf(gcf,'Fig1wtSpread','type',{'fig'},'paper','gilad22','size',[])
%% spread wt vs. mytants
%close all
types={'','greys','greens','blues','Browns'}
o=1
intTps=[15,30,40,70,90,110,130,150,180]
for r=2:3
    figure
    c=0;
    for ab=intersect(unique(smeta.ab(ismember(smeta.gt,showRpts.Var1(r)))),{'myc','k9ac'})'
        c=c+1;
        subplot(2,2,c)
        for t=find(contains(tcs.gt,showRpts.Var1([1,r]))&contains(tcs.ab,ab)&contains(tcs.exp,setdiff(goodExps,'c51.mat')))'%find(contains(tcs.gt,{'wt','rtt'})&contains(tcs.ab,{'myc'})&contains(tcs.exp,setdiff(goodExps,'c51.mat')))'
            selSmp=find(smeta.tid==t&contains(smeta.exp,tcs.exp{t})&ismember(smeta.tp,intTps));
            tpCols=brewermap(numel(selSmp)+5,types{tcs.abId(t)});
            if ~contains(tcs.gt(t),'wt')
                scatter(absCorr(o,selSmp),-diffCorr(o,selSmp),rescale(xCorrAbsDiff(o,selSmp),70,240,'InputMax',0,'InputMin',-1), tpCols(1:end-5,:),'filled','DisplayName',tcs.exp{t})
            else
                scatter(absCorr(o,selSmp),-diffCorr(o,selSmp),rescale(xCorrAbsDiff(o,selSmp),70,240,'InputMax',0,'InputMin',-1), tpCols(1:end-5,:),'DisplayName',tcs.exp{t},'Linewidth',2)
            end
            hold on
        end
        axis tight
        xlabel('correlation with "DNA"')
        ylabel('correlation with d/dx DNA')
        scatter([0:0.2:.6],[0 0 0 0 ],rescale([-1:1/3:0],70,240,'InputMax',0,'InputMin',-1),[0 0 0],'Linewidth',2)
        plot(xlim,xlim,'k--')
        title(ab{1})
    end   
    sgtitle(sprintf('wt vs. %s',showRpts.Var1{r}))
    save_gf(gcf,sprintf('%sVsWtSpread',showRpts.Var1{r}),'type',{'fig'},'paper','gilad22','size',[])
end

%% cluster based on tor for cluster-based analysis, e.g. Fig 2C
clearvars -except smeta tcs binData binMeta binResid goodExps showRpts
dnaCluster=load('DNAcluster.mat');
[~,binIdx]=min(abs(binMeta.pos-dnaCluster.meta(:,1)'),[],2);
binIdx=dnaCluster.finalClu(binIdx);
binIdx(ismember(binIdx,find(accumarray(binIdx(binIdx>0),1)<20)))=NaN;

clusterData=nan(max(binIdx),size(binData,2));
clusterRes=nan(max(binIdx),size(binResid,2));

for b=1:max(binIdx)    
    clusterData(b,:)=median(binData(binIdx==b&binMeta.bad==0,:),'omitnan');
    clusterRes(b,:)=median(binResid(binIdx==b&binMeta.bad==0,:),'omitnan');
end
clusterTor=accumarray(binIdx(binIdx>0),binMeta.tor(binIdx>0),[],@(x)median(x,'omitnan'),nan);
clusterSize=accumarray(binIdx(binIdx>0),1)
binMeta.binIdx=binIdx;
clear dnaCluster binIdx

%% absolute levels
nClu=max(binMeta.binIdx)
dTot=zeros(1,size(clusterData,2))
binForTot={[1;3],[nClu+[-1:0]]}
close all
c=0;
for g=find(ismember(tcs.ab,{'ha','DNA'})&accumarray(smeta.tid,smeta.tp,[],@max)>5)'
    idx=find(smeta.tid==g&smeta.tp>-5 & ~smeta.bad);      
    [~,binForTot{2}]=mink(sum(clusterRes(:,idx),2),3);
    intMeans=[0 mean(log2(clusterData(binForTot{1},idx)),1);0 mean(log2(clusterData(binForTot{2},idx)),1)];
    dMeans=diff(intMeans,1,2)';
    [~,turningPoint]=min(movsum(dMeans(:),size(dMeans,1),'Endpoints','discard'));
     dTot(1,idx)=cumsum([max(-dMeans([numel(idx)+1:numel(idx)+turningPoint-1,turningPoint:numel(idx)]),0)]);       
end
smeta.dTot=dTot';
absClu=clusterData.*2.^smeta.dTot';

%absCluF=NaN(nClu,size(smeta,1));
absCludT=NaN(nClu,size(smeta,1));
fit=0
for g=find(ismember(tcs.ab,{'ha','DNA'})&accumarray(smeta.tid,smeta.tp,[],@max)>5)'
    if fit==0
        selSmp=find(smeta.tid==g&smeta.tp>-5&~smeta.bad);
        temp=diff([ones(nClu,1) absClu(:,selSmp)],1,2)./diff([0;smeta.tp(selSmp)])'*100;
        absCludT(:,selSmp)=temp/2+temp(:,[2:end,end])/2;
        %absCluF(:,selSmp)=absClu(:,selSmp);
    else
        selSmp=find(smeta.tid==g&smeta.tp>-5);
        [temp,tempdT,tempFit]=fitAbsHA([0,smeta.tp(selSmp(smeta.bad(selSmp)==0))'],[ones(nClu,1) absClu(:,selSmp(smeta.bad(selSmp)==0))],'figure',false);
        for b=find(~cellfun('isempty',tempFit))
            absCludT(b,selSmp)=differentiate(tempFit{b},smeta.tp(selSmp));
            absCluF(b,selSmp)=tempFit{b}(smeta.tp(selSmp));
        end        
    end
end
for s=find(smeta.bad==1&contains(smeta.ab,'ha'))'
    cmpTc=[find(smeta.tid==smeta.tid(s)&smeta.tp<smeta.tp(s),1,'last'),find(smeta.tid==smeta.tid(s)&smeta.tp>smeta.tp(s),1)]
    weight=abs(smeta.tp(cmpTc)-smeta.tp(s))
    weight=(sum(weight)-weight)./sum(weight)
    absClu(:,s)=sum(absClu(:,cmpTc).*weight',2);
    clusterData(:,s)=sum(clusterData(:,cmpTc).*weight',2);
end
clearvars b binFor binIdx c dMeans dnaCluster dTot sTot fit temp turningPoint idx g intMeans selSmp weight cmptTc

%% Thsi creates the  replacement vs. replication rate or replicated fraction plot, e.g. Fig 2C
bgCols=(brewermap(150,'greens'))%brewermap(100,'oranges');
dnaTp=unique(smeta.tp(ismember(smeta.ab,'DNA')&smeta.tp>0))
close all
plotTp=[15,30,40,50,70,90,110,130,150]%[70,90,110,130,150]%
plotCols=5;
cmpSmp={'DNA','c51.mat'}
for r=1:size(showRpts,1)
    selSmp=find(ismember(smeta.gt,showRpts.Var1{r})&ismember(smeta.tp,dnaTp)&contains(smeta.exp,showRpts.Var2{r})&ismember(smeta.ab,'ha'));
    if numel(selSmp)>0
        figure('Position',[1 41 1920/1.53 962], 'color','w');;
    end
    c=0;
    for s=selSmp'
        sMyc=find(contains(smeta.name,regexprep(smeta.name{s},{'^ha','*'},{'myc',''}))&strcmp(smeta.exp,smeta.exp{s}));
        sMe=find(smeta.tp==smeta.tp(s)&ismember(smeta.gt,smeta.gt(s))&ismember(smeta.ab,cmpSmp(1))&ismember(smeta.exp,cmpSmp(2)));
        if numel(sMyc)==1 & numel(sMe)==1           
            cVal=clusterTor;%absClu(:,sMe);
            xVal=absCludT(:,sMe);
            yVal=exp(log2(clusterData(:,sMyc))-log2(clusterData(:,s)));
            [~,topBin]=maxk(xVal.*(absClu(:,sMe)>1.1),1);
            topBin=min(topBin)
            selbin=find((absClu(:,sMe)<=absClu(topBin,sMe))&~isnan(xVal)&~isnan(yVal));            
            selbin=selbin((yVal(selbin)>=min(maxk(yVal(selbin),70)))&(xVal(selbin)>=min(maxk(xVal(selbin),70))))
            toFit(s,:)=robustfit(xVal(selbin),yVal(selbin));        
            if ismember(smeta.tp(s),plotTp)
                c=c+1;
                subplot(4,plotCols,c+2*plotCols)
                hold off
                scatter(xVal(:),yVal,[],cVal,'filled','Linewidth',1,'MarkerEdgeColor',[0 0 0])
                hold on
                scatter(xVal(1),yVal(1),[],[1,0,0],'Linewidth',2,'MarkerEdgeColor',[1 0 0])
                axis tight
                plot(xlim,xlim.*toFit(s,2)+toFit(s,1),'k-')
                ylabel('myc/HA')
                xlabel(sprintf('rep Rte %s',cmpSmp{1}))
                title(strjoin(smeta{s,{'gt','con','exp'}}))
                bad(s)=false;
                text(xlim*[.95;.05],ylim*[0;0.95],sprintf('%.2fx%+.2x',toFit(s,2)))
                set(gca,'Color',bgCols(round(rescale(max(xVal),1,50,'InputMin',0,'InputMax',2)),:))
                colormap(gca,flipud(brewermap(128,'buPu')))
                %xlim([-.15 1.6])
            end
            xVal=absClu(:,sMe)-1;
            yVal=exp(log2(clusterData(:,sMyc))-log2(clusterData(:,s)));
            [~,topBin]=max(xVal);
            selbin=absClu(:,sMe)<absClu(topBin,sMe);            
            toFitAbs(s,:)=robustfit(xVal(selbin),yVal(selbin));           
            if ismember(smeta.tp(s),plotTp)
                subplot(4,plotCols,c+0*plotCols)
                scatter(xVal(:),yVal,[],cVal,'filled','Linewidth',1,'MarkerEdgeColor',[0 0 0])
                hold on
                scatter(xVal(1),yVal(1),[],[1,0,0],'Linewidth',2,'MarkerEdgeColor',[1 0 0])
                hold on
                axis tight
                ylabel('myc/HA')
                xlabel(sprintf('replicated fraction %s',cmpSmp{1}))
               title(strjoin(smeta{s,{'gt','con','exp'}}))
                bad(s)=false;
                text(xlim*[.95;.05],ylim*[0;0.95],sprintf('%.2fx%+.2x',toFitAbs(s,2)))
                set(gca,'Color',bgCols(round(rescale(max(xVal),1,50,'InputMin',0,'InputMax',2)),:))
                colormap(gca,flipud(brewermap(128,'buPu')))
                xlim([-.1 .8])
            end            
        else
            bad(s)=true;
        end
    end
    save_gf(gcf,sprintf('%s_turnoverSup',showRpts.Var1{r}),'type',{'fig'},'paper','gilad22','size',[])
end


%% slope calculation - for Figure 2D
toFit=nan(size(smeta,1),2);
smeta.tp(smeta.tp==190)=180;
cmpSmp={'DNA','c51.mat'}

for s=find(ismember(smeta.exp,setdiff(goodExps,'c51.mat'))&ismember(smeta.ab,'ha')&smeta.bad>-1&smeta.tp>15)'
    sMyc=find(contains(smeta.name,regexprep(smeta.name{s},{'^ha','*'},{'myc',''}))&strcmp(smeta.exp,smeta.exp{s}));
    sMe=find(smeta.tp==smeta.tp(s)&ismember(smeta.gt,smeta.gt(s))&ismember(smeta.ab,cmpSmp(1))&ismember(smeta.exp,cmpSmp(2)));
    if numel(sMyc)==1 & numel(sMe)==1
        cVal=clusterTor;%absClu(:,sMe);
        xVal=absCludT(:,sMe);
        yVal=exp(log2(clusterData(:,sMyc))-log2(clusterData(:,s)));
        [~,topBin]=maxk(xVal.*(absClu(:,sMe)>1.1),1);
        topBin=min(topBin)
        selbin=find((absClu(:,sMe)<=absClu(topBin,sMe))&~isnan(xVal)&~isnan(yVal));
        selbin=selbin((yVal(selbin)>=min(maxk(yVal(selbin),70)))&(xVal(selbin)>=min(maxk(xVal(selbin),70))))
        toFit(s,:)=robustfit(xVal(selbin),yVal(selbin));
    end
end

%% this create the slope plot - e.g. Figure 2D

%%meanplots
[~,smeta.gtid]=ismember(smeta.gt,showRpts.Var1);
gtColScheme={'Blues','Purples','Greens'}
gtCols=cell2mat(arrayfun(@(r)brewermap(1,gtColScheme{r}),[1:3]','UniformOutput',false))%lines(size(showRpts,1));


repFrac=sum(clusterSize.*absClu,'omitnan')./sum(clusterSize);
plotTp=unique(smeta.tp(ismember(smeta.ab,'DNA')&smeta.tp>10))
figure
for r=1:size(showRpts,1)
    selSmp=find(contains(smeta.gt,showRpts.Var1{r})&ismember(smeta.tp,plotTp)&smeta.bad>-1&contains(smeta.ab,'ha')&contains(smeta.exp,setdiff(goodExps,'c51.mat')))
    [uTps,~,tpId]=unique(smeta.tp(selSmp));
    meanVal=accumarray(tpId,toFit(selSmp,2),[],@mean,NaN);
    stdError=accumarray(tpId,toFit(selSmp,2),[],@(x)std(x)/sqrt(numel(x)-1),NaN)
    uTpsP=uTps(~isnan(stdError));
    meanVal=meanVal(~isnan(stdError));
    stdError=stdError(~isnan(stdError));
    plot(uTpsP,meanVal,'Linewidth',2,'Color',brighten(gtCols(smeta.gtid(selSmp(1)),:),-.4),'Marker','.')
    hold on
    fill([uTpsP;flip(uTpsP)],[meanVal-stdError;flip(meanVal+stdError)],[gtCols(smeta.gtid(selSmp(1)),:)],'FaceAlpha',.5,'LineStyle','none')
    axis tight
    ylabel('slop vs. Rate')
    xlabel('time in HU')
    ylim([0 max(ylim)])
end
save_gf(gcf,sprintf('Fig2SlopesumWtRttVps'),'type',{'fig'},'paper','gilad22','size',[])


%% this plot compares different time courses e.g. Fig 3C
intTps=[30,40,70,90,110,130,150]
gtColScheme={'Blues','Purples','Greens'}
figure
for r=2:3
    subplot(2,2,1)
            selWts=find(contains(smeta.gt,'wt')&ismember(smeta.tp,intTps)&smeta.bad==0&contains(smeta.ab,'ha')&contains(smeta.exp,setdiff(goodExps,'c51.mat')))
            selMts=find(contains(smeta.gt,showRpts.Var1{r})&ismember(smeta.tp,intTps)&smeta.bad==0&contains(smeta.ab,'ha')&contains(smeta.exp,setdiff(goodExps,'c51.mat')))
            [~,tpIdWt]=ismember(smeta.tp(selWts),intTps)
            [~,tpIdMt]=ismember(smeta.tp(selMts),intTps)
            %selDna=find(contains(smeta.gt,showRpts.Var1([1,r]))&ismember(smeta.tp,intTps)&smeta.bad==0&contains(smeta.ab,'DNA')&contains(smeta.exp,'c51.mat'))
            selDna=find(contains(smeta.gt,showRpts.Var1([1,r]))&ismember(smeta.tp,intTps)&smeta.bad==0&contains(smeta.ab,'ha')&~contains(smeta.exp,'c51.mat'))
            [~,tpIdDna]=ismember(smeta.tp(selDna),intTps)
            
            plotVal=@(x)x;
            
            meanValWt=accumarray(tpIdWt,plotVal(toFit(selWts,2)),[],@mean,NaN);
            stdErrorWt=accumarray(tpIdWt,plotVal(toFit(selWts,2)),[],@(x)std(x)/sqrt(numel(x)-1),NaN)
            meanValMt=accumarray(tpIdMt,plotVal(toFit(selMts,2)),[],@mean,NaN);
            stdErrorMt=accumarray(tpIdMt,plotVal(toFit(selMts,2)),[],@(x)std(x)/sqrt(numel(x)-1),NaN)
            %scCols=brewermap(numel(meanValWt)+3,gtColScheme{r});
            meanFrac=accumarray(tpIdDna,absClu(1,selDna),[],@(x)mean(x,'omitnan'),NaN)
            if r==2
                cScheme=brewermap(75,gtColScheme{r});
                cScheme=cScheme([1:3:39,39:75],:)
            else
                cScheme=brewermap(50,gtColScheme{r});
            end
            scCols=cScheme(round(rescale(meanFrac-1,0.5,50.4,'InputMin',0.5,'InputMax',1)),:)
            
            errorbar(meanValWt,meanValMt,stdErrorMt,stdErrorMt,stdErrorWt,stdErrorWt,'LineStyle','none','Color',scCols(end-3,:),'Linewidth',1,'CapSize',0)
            hold on
            scatter(meanValWt,meanValMt,(7+[1:numel(intTps)]).^2,scCols,'filled','MarkerEdgeColor',scCols(end-3,:),'DisplayName',showRpts.Var1{r})
            text(meanValWt,meanValMt,num2str(intTps'))
            xlim([0 4])
            ylim([0 4])
            plot(xlim,xlim,'k--')
            xlabel('(Wt rate)')
            ylabel('(mt rate)')
            subplot(2,2,r)            
            
            scatter(1:2:7,zeros(1,4),(7+[1:2:numel(intTps)]).^2,'k')
            text(1:2:7,zeros(1,4),num2str(intTps(1:2:7)'))
            caxis([1.5 2]-1)
            colormap(gca,cScheme)
            colorbar()
end
save_gf(gcf,sprintf('Fig3SlopeXsumRttVps'),'type',{'fig'},'paper','gilad22','size',[])


%% myc along time for S2
cMap=[flipud(brewermap(64,'purples'));(brewermap(64,'blues'))]
figure;
c=0;
cmpSmp=627
for d=[1:2]
    for t=find(contains(tcs.ab,'myc')&contains(tcs.gt,'wt')&contains(tcs.ab,'myc')&ismember(tcs.exp,goodExps))
        c=c+1;
        subplot(2,10,(c-1)*5+[1:4])
        selSmp=find(ismember(smeta.tid,t)&smeta.tp>=-5);
        [~,idx]=sortrows(table(smeta.tp(selSmp),str2double(regexp(smeta.exp(selSmp),'\d+','match','once'))));
        selSmp=selSmp(idx);
        [~,haIdx]=ismember([smeta(selSmp,{'gt','tp'}),repmat(smeta(cmpSmp,{'ab','exp'}),numel(selSmp),1)],smeta(:,{'gt','tp','ab','exp'}))
        [~,mIdx(haIdx>0)]=max(movmean(absCludT(1:2:end,haIdx(haIdx>0)),3),[],1)
        hold off
        if d==1
            imagesc((movmean(clusterRes(1:2:end,selSmp),1)),[-1.2 1.2])
             title('FC vs. G1 - log2')
        else
            imagesc((movmean(log2(clusterData(1:2:end,selSmp)),1)),[-1.2 1.2])
            title('log2 norm Data')
        end
        hold on
        scatter(1:numel(mIdx),mIdx,[],[.5 .5 .5],'o','filled','MarkerFaceAlpha',0.5)
        colormap(gca,cMap(21:108,:))
        ylabel('bin order by ToR')
        xlabel('Timepoint')
        xticks(movmean([0;find(smeta.tp(selSmp(1:end-1))~=smeta.tp(selSmp(2:end)));numel(selSmp)],2,'Endpoints','discard')+0.5)
        xticklabels(unique(smeta.tp(selSmp)))
        title(strjoin(tcs{t(1),1:2}))
        ylabel(colorbar(),'(log2(myc))')
    end
end
save_gf(gcf,sprintf('Fig1AMycrepl'),'type',{'fig'},'paper','gilad22','size',[])


%% turnover correlation
% calculate turnover all wt samples
selMyc=find(ismember(smeta.ab,'myc')&ismember(smeta.gt,'wt')&ismember(smeta.tp,smeta.tp(ismember(smeta.exp,'t51.mat'))))
haCmp=180
[~,selHa]=ismember([repmat(smeta(haCmp,'ab'),numel(selMyc),1),smeta(selMyc,{'gt','exp','tp'})],smeta(:,{'ab','gt','exp','tp'}))
turnOver=log2(clusterData(:,selMyc))-log2(clusterData(:,selHa));
dnaCmp=627;
[~,selDna]=ismember([repmat(smeta(dnaCmp,{'ab','exp'}),numel(selMyc),1),smeta(selMyc,{'gt','tp'})],smeta(:,{'ab','exp','gt','tp'}))
crAbs=diag(corr(turnOver,absClu(:,selDna),'rows','pairwise'));
crRate=diag(corr(turnOver,absCludT(:,selDna),'rows','pairwise'));
figure
for t=unique(smeta.tid(selMyc))'
    selSmp=smeta.tid(selMyc)==t
    plot(smeta.tp(selMyc(selSmp)),crAbs(selSmp),'--','COlor',.7*[1 1 1])
    hold on
    plot(smeta.tp(selMyc(selSmp)),crRate(selSmp),'-','COlor',.0*[1 1 1])
end
[utps,~,uid]=unique(smeta.tp(selMyc))
mCrAbs=accumarray(uid,crAbs,[],@mean)
mCrRate=accumarray(uid,crRate,[],@mean)
plot(utps,mCrAbs,'--','COlor',.7*[1 1 1],'Linewidth',2)
plot(utps,mCrRate,'-','COlor',.0*[1 1 1],'Linewidth',2)
xlim([15 180])
ylim([0,1])
save_gf(gcf,sprintf('FigS2CrFig'),'type',{'fig'},'paper','gilad22','size',[])

    
%% calcylate temporal correlation as done in Fig. 1E and similiar use the 2000bp bined data for this part

clear all
[data,smeta]=getNaamaData('figure',26);
smeta.exp(cellfun('prodofsize',smeta.folder)==54)=strrep(smeta.exp(cellfun('prodofsize',smeta.folder)==54),'x','s');
smeta.exp(ismember(smeta.exp,'s51.mat')&contains(smeta.tag,'h3x2b'))={'t51.mat'}
smeta.gt(ismember(smeta.gt,'rtt109'))={'rtt'}
smeta.tp=str2double(regexp(smeta.con,'(?<=HU)\d+','match','once'));
smeta.tp(contains(smeta.con,'async'))=-10;
smeta.tp(contains(smeta.con,'alpha'))=-5;
goodExps={'s39.mat';'s43.mat';'s44.mat';'s47.mat';'x35.mat';'x36.mat';'s51.mat';'t51.mat';;'c51.mat'}
keepSmp=ismember(smeta.exp,goodExps)
smeta=smeta(keepSmp,:);
data=data(:,keepSmp);

binSize=2000;
[binData,binMeta,~]=DNA_norm_bin(data,'bin_size',binSize,'calc_ribo_frac',true,'indChr',false,'format','abs','clean',false,'n_ori',358,'n_chr',17,'Tn5',true);
clear data

c=0
for ab={'ha','DNA'}
    for g={'wt','rtt','vps75','mecsml'}
        for tp=unique(smeta.tp(ismember(smeta.ab,ab)&ismember(smeta.gt,g)))'
            c=c+1;
            selSmp=find(ismember(smeta.gt,g)&smeta.tp==tp&ismember(smeta.ab,ab)&smeta.bad==0);
            meanBin(:,c)=mean(binData(:,selSmp),2);
            mSmeta(c,:)=smeta(selSmp(1),:);
            mSmeta.exp(c)={'c51.mat'}
        end
    end
end
smeta=[smeta;mSmeta];
binData=[binData,meanBin];
clearvars -except smeta binData binMeta goodExps

[tcs,~,smeta.tid]=unique(smeta(:,{'gt','exp','ab'}),'stable')
[smeta,idx]=sortrows(smeta,{'tid','tp'});
binData=binData(:,idx);
clearvars idx
GP=load('/home/labs/barkailab/felixj/Documents/MATLAB/scripts/gene/group_imp.mat');
binMeta=array2table(binMeta,'VariableNames',{'pos','tor','oriDis','oriIdx','chr','bad'});

%% start here when loading the prepared data ()
binResid=nan(size(binData));
for t=1:size(tcs,1)%find(contains(tcs.ab,{'9ac','k56ac','myc','ha'})&contains(tcs.exp,{'31','36'}))'
    selSmp=find(smeta.tid==t);
    if ismember(tcs.ab(t),{'ha','DNA'})
        alphaSmp=find(ismember(smeta.ab,tcs.ab(t))&ismember(smeta.exp,setdiff(goodExps,'c51.mat'))&smeta.tp==-5);
    else
        alphaSmp=find(ismember(smeta(:,{'ab','gt'}),tcs(t,{'ab','gt'}))&smeta.tp==-5&ismember(smeta.exp,setdiff(goodExps,'c51.mat')));
    end
    if numel(alphaSmp)>0
        for s=selSmp'
            [p,r]=robustfit(log2(mean(binData(binMeta.bad==0,alphaSmp),2)+0.1),log2(0.1+mean(binData(binMeta.bad==0,s),2)));
            binResid(binMeta.bad==0,s)=r.resid;
        end
    else
        tcs(t,:)
    end
end
clearvars t ans p r s selSmp t
%% prepare fits
clear valHRes diffHRes
maxTp=180;
fitTps=0:5:maxTp;
parfor t=1:size(tcs,1)
    selSmp=find(smeta.tid==t&smeta.tp>-10&smeta.bad==0&smeta.tp<maxTp);
    tempMat=smoothdata(binResid(:,selSmp),'rlowess',6);
    bFit=cell(12157,1);
    %valHRes{t}=nan(12157,37);
    diffHRes{t}=nan(size(binMeta,1),numel(fitTps));
    for b=1:size(tempMat,1)  
        try
        bFit{b}=fit(smeta.tp(selSmp),tempMat(b,:)','smoothingspline');
                %valHRes{t}(b,:)=bFit{b}(0:5:180);
        diffHRes{t}(b,:)=differentiate(bFit{b},(fitTps));
        end
    end
end
%%  compare fits
abOrder={'ha','k9ac','k56ac','myc'}
cMax1=0.75;
cMax2=1;
getPos=@(x,y)x(y);
startTP=7;
endTP=30;%180/5+1;
tauRange=15;
spanMin=endTP-startTP-tauRange+1;
corrLine=nan(size(tcs,1),tauRange*2+1);
close all
figure
cc=0;
cmpPars={'DNA','c51.mat'}
clear saveTable
for t1=find(ismember(tcs.ab,cmpPars(1))&ismember(tcs.exp,cmpPars(2)))'
    for e=unique(tcs.exp(ismember(tcs.gt,tcs.gt(t1))&ismember(tcs.exp,setdiff(goodExps,'c51.mat'))&ismember(tcs.ab,abOrder)))'
        intTc=find(contains(tcs.gt,tcs.gt(t1))&ismember(tcs.exp,e)&ismember(tcs.ab,abOrder));
        [~,idx]=ismember(tcs.ab(intTc),abOrder);
        [~,idx]=sortrows([tcs(intTc,'exp'),table(idx)])
        intTc=intTc(idx);
        cc=cc+1;
        imageMat=nan(2*tauRange+1,numel(intTc)-1);
        c=0;
        normVec=[];
        crType=[];
        for t2=setdiff(intTc',t1,'stable')
            c=c+1;
            corrMat=max(corr(diffHRes{t1}(binMeta.bad==0,7:20),diffHRes{t2}(binMeta.bad==0,7:end),'rows','pairwise'),[0]);
            imageMat(:,c)=arrayfun(@(x)mean(diag(corrMat,x)),[-tauRange:tauRange]);
            %corrMat=max(corr(diffHRes{t1}(binMeta.bad==0,startTP:endTP),diffHRes{t2}(binMeta.bad==0,startTP:endTP)),0);
            %imageMat(:,c)=arrayfun(@(k)mean((diag(corrMat,k))),-tauRange:tauRange,'UniformOutput',true);
            corrLine(t2,:)=imageMat(:,c);
            if contains(tcs.ab(t2),'k9')
                normVec(1,c)=cMax1;
            else
                normVec(1,c)=cMax2;
            end
            [~,crType(1,c)]=ismember(tcs.ab(t2),abOrder);
        end
        [maxVal,maxTp]=max(imageMat);
        maxTp=(maxTp-tauRange-1)*5;
        subplot(4,5,cc)
        xlabel('tau')
        ylabel('cr(tau)')
        imagesc(imageMat./normVec,'YData',5*[-tauRange,tauRange],[0 1])
        hold on;
        plot(xlim,[0 0],'k-','linewidth',2)
        scatter(1:numel(maxTp),maxTp,[],[0 0 0],'filled')
        ylabel('\tau')
        xlabel([])
        title(strrep(sprintf('%s vs %s',strjoin(tcs{t1,:}),e{1}),'.mat',''))
        %colorbar()
        set(gca,'xtick',1:numel(maxTp),'xticklabel',strcat(tcs.ab(setdiff(intTc',t1,'stable')),'-',extractBefore(tcs.exp(setdiff(intTc',t1,'stable')),'.mat')))
        set(gca,'Position',get(gca,'Position').*[1 1 1/5*numel(normVec) 1])
        colormap(gca,brewermap(128,'blues'))
        ylim(min(67.5,tauRange*5).*[-1 1])
        %     pause(.5)
        saveTable{cc}=[tcs(setdiff(intTc',t1,'stable'),:),table(maxVal',maxTp','VariableNames',{'maxCorr','maxTP'})];
    end
end
saveTable=cat(1,saveTable{:});
save(sprintf('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Gilad-CC/yoav-%dbp.mat',binSize),'saveTable')
save_gf(gcf,'YoavAll','type',{'fig'},'paper','gilad22','size',[])
saveTable.Properties.VariableNames=strrep(saveTable.Properties.VariableNames,'max','max2k')

bp500=load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Gilad-CC/yoav-500bp.mat')
[~,oldId]=ismember(saveTable(:,1:3),bp.saveTable(:,1:3))
saveTable=[saveTable,bp500.saveTable(oldId,4:5)]
figure
scatter(saveTable.max2kTP,saveTable.maxTP,100,saveTable.max2kCorr-saveTable.maxCorr,'filled')
hold on
colormap(gca,brewermap(128,'purples'))
xlabel('time delay 2kb-bins')
ylabel('time delay 500bp-bins')
ylabel(colorbar(),'\Delta Cr')
save_gf(gcf,'YoavBinSize','type',{'fig'},'paper','gilad22','size',[])


%% calcuate ToR and compare

tempData=binData;%smoothdata(binData,[],)
tempData(binMeta.bad==1,:)=nan;

dnaCluster=load('DNAcluster.mat');
[~,binIdx]=min(abs(binMeta.pos-dnaCluster.meta(:,1)'),[],2);
binIdx=dnaCluster.finalClu(binIdx);
binIdx(ismember(binIdx,find(accumarray(binIdx(binIdx>0),1)<20)))=NaN;

clusterTemp=nan(max(binIdx),size(tempData,2));

for b=1:max(binIdx)    
    clusterTemp(b,:)=median(binData(binIdx==b&binMeta.bad==0,:),'omitnan');
end
clusterTor=accumarray(binIdx(binIdx>0),binMeta.tor(binIdx>0),[],@(x)median(x,'omitnan'),nan);
clusterSize=accumarray(binIdx(binIdx>0),1)
binMeta.binIdx=binIdx;
clear dnaCluster binIdx

%% absolute lebels
nClu=max(binMeta.binIdx)
dTot=zeros(1,size(clusterTemp,2))
binForTot={[1;3],find(clusterSize>20,2,'last')}
close all
c=0;
for g=find(ismember(tcs.ab,{'ha','DNA'})&accumarray(smeta.tid,smeta.tp,[],@max)>5)'
    idx=find(smeta.tid==g&smeta.tp>-5 & ~smeta.bad);      
    [~,binForTot{2}]=mink(sum(clusterTemp(:,idx),2),3);
    intMeans=[0 mean(log2(clusterTemp(binForTot{1},idx)),1);0 mean(log2(clusterTemp(binForTot{2},idx)),1)];
    dMeans=diff(intMeans,1,2)';
    [~,turningPoint]=min(movsum(dMeans(:),size(dMeans,1),'Endpoints','discard'));
     dTot(1,idx)=cumsum([max(-dMeans([numel(idx)+1:numel(idx)+turningPoint-1,turningPoint:numel(idx)]),0)]);       
end
smeta.dTot=dTot';

smData=smoothdata(tempData,1,'movmean',5);
selSmp=smeta.tid==45
selTp=max(smeta.tp(selSmp),0)
totBin=log2(smData(:,selSmp)')+smeta.dTot(selSmp);

goodBin=find(binMeta.bad==0 & binMeta.tor>0 & all(maxk(totBin,2)>0.2,1)'& (min(totBin,[],1)<0.09)')
repFitR=cell(numel(goodBin),1);
repTimeR=nan(size(goodBin));

%%
for p=1:numel(repFitR)%15:226%find(contains(smeta.ab,'ha')&smeta.tp>10)'    
    b=goodBin(p);
    try
    temp=L5P(selTp,max(movmean(totBin(:,b),3),0));
    repFitR{p}=temp;
    repTimeR(p)=max([find(temp(0:200)>0.2,1),nan]);
    end
end

repTime=nan(size(binMeta,1),1);
repTime(goodBin)=repTimeR;
oriBin=find(abs(binMeta{:,'oriDis'})<1000)%&[1:size(binMeta,1)]'>surLen&[1:size(binMeta,1)]'<(size(binMeta,1)-surLen) & ~ismember(binMeta.oriIdx,find(contains(GP.oris.Var4,'1216'))))

figure;

cMap=flipud(brewermap(128,'purples'))
subplot(2,2,1)
hold off
clearvars nRep
nRep(1,:)=histcounts(binMeta.tor(~isnan(repTime)),[16.5:41.5]);
nRep(2,:)=histcounts(binMeta.tor(isnan(repTime)),[16.5:41.5]);
nPlot=cumsum(cumsum(nRep,2),1)
plot([17:41]',nPlot(2,:),'Linewidth',2,'DisplayName','all annotated','Color',cMap(90,:))
hold on
plot([17:41]',nPlot(1,:),'Linewidth',2,'DisplayName','replicated','Color',cMap(5,:))
xlabel('RT - Yabuki et al.')
ylabel('# locis')
axis tight

subplot(2,2,3)
hold off
clearvars nRep
nRep(1,:)=histcounts(binMeta.tor(~isnan(repTime)&ismember([1:size(repTime,1)]',oriBin)),[16.5:41.5]);
nRep(2,:)=histcounts(binMeta.tor(isnan(repTime)&ismember([1:size(repTime,1)]',oriBin)),[16.5:41.5]);
%bar([17:41]',cumsum(nRep'),'stacked')
nPlot=cumsum(cumsum(nRep,2),1)
plot([17:41]',nPlot(2,:),'Linewidth',2,'DisplayName','all annotated','Color',cMap(90,:))
hold on
plot([17:41]',nPlot(1,:),'Linewidth',2,'DisplayName','replicated','Color',cMap(5,:))
xlabel('RT - Yabuki et al.')
ylabel('# ORIs')
axis tight

subplot(1,2,2)
hold off
plotBin=~isnan(repTime)
ds=dscatter(binMeta.tor(plotBin),repTime(plotBin));
ds.Children.Marker='o';
ds.Children.SizeData=20;
ds.Children.DisplayName=sprintf('loci cr:%.2f',corr(binMeta.tor(plotBin),repTime(plotBin)));
hold on
scatter(binMeta.tor(oriBin),repTime(oriBin),[20],.2*[1 1 1],'filled','DisplayName',sprintf('oris cr:%.2f',corr(binMeta.tor(oriBin),repTime(oriBin),'rows','pairwise')))
axis tight
xlabel('RT - Yabuki et al.')
ylabel('RT - wt in HU')
colormap(cMap)
save_gf(gcf,sprintf('FigS1WTHURep'),'type',{'fig'},'paper','gilad22','size',[])