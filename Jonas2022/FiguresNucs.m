%% Figures nucleosomes this script is used to create  Figure 3D, S2C, S3D, S4F and S5A
%% nucleosmes
clear all
[~,~,smeta,nucData,nucMeta]=standardImport(26,'save',true);
%load('forFigure26.mat','nucData','smeta','nucMeta')
goodExps={'s39.mat';'s43.mat';'s44.mat';'s47.mat';'x35.mat';'x36.mat';'s51.mat';'t51.mat';'c51.mat'}
smeta.exp(cellfun('prodofsize',smeta.folder)==54)=strrep(smeta.exp(cellfun('prodofsize',smeta.folder)==54),'x','s');
smeta.exp(ismember(smeta.exp,'s51.mat')&contains(smeta.tag,'h3x2b'))={'t51.mat'};
smeta.gt(ismember(smeta.gt,'rtt109'))={'rtt'}

keep=contains(smeta.exp,goodExps);
smeta=smeta(keep,:);
nucData=nucData(:,keep);
clearvars keep
smeta.tp=str2double(regexp(smeta.con,'(?<=HU)\d+','match','once'));
smeta.tp(contains(smeta.con,'async'))=-10;
smeta.tp(contains(smeta.con,'alpha'))=-5;
c=0;
clear mSmeta meanBin
for ab={'ha','DNA'}
    for gt=unique(smeta.gt(ismember(smeta.ab,ab)&smeta.bad==0&ismember(smeta.gt,[{'mecsml'}    {'rtt'}    {'vps75'}    {'wt'}])))'
        for tp=unique(smeta.tp(ismember(smeta.gt,gt)&ismember(smeta.ab,ab)&smeta.bad==0))'
            selSmp=find(smeta.tp==tp&ismember(smeta.gt,gt)&ismember(smeta.ab,ab)&smeta.bad==0)
            c=c+1;
            meanNuc(:,c)=mean(nucData(:,selSmp),2);
            mSmeta(c,:)=smeta(selSmp(1),:);
        end
    end
end
mSmeta.exp(:)={'c51.mat'}
nucData=[nucData,meanNuc];
smeta=[smeta;mSmeta];
clearvars -except nucData nucMeta smeta goodExps

[tcs,~,smeta.tid]=unique(smeta(:,{'gt','exp','ab'}),'stable')
[smeta,idx]=sortrows(smeta,{'tid','tp'});
nucData=nucData(:,idx);
GP=load('/home/labs/barkailab/felixj/Documents/MATLAB/scripts/gene/group_imp.mat');


clearvars -except nucData smeta nucMeta goodExps GP tcs

goodHa=mean(nucData(:,contains(smeta.ab,'ha')&~smeta.bad)>5,2)>0.9; % above 5 in 90% of the nucs
dNucs=nan(size(nucData));
for t=find(contains(tcs.exp,goodExps))'
    selSmp=find(smeta.tid==t);
    if ismember(tcs.ab(t),{'ha','DNA'})
        alphaSmp=find(ismember(smeta.ab,tcs.ab(t))&ismember(smeta.exp,goodExps)&smeta.tp==-5);
    else
        alphaSmp=find(ismember(smeta.gt,tcs.gt(t))&ismember(smeta.ab,tcs.ab(t))&smeta.tp==-5)
    end
    for s=selSmp'
        [p,r]=robustfit(log2(mean(nucData(goodHa,alphaSmp),2)+0.1),log2(0.1+mean(nucData(goodHa,s),2)));
        dNucs(goodHa,s)=r.resid;
    end
end

xVals=log2(nucData(:,smeta.tid==39)+1)%-log2(nucData(:,smeta.tid==16)+1);
yVals=log2(nucData(:,smeta.tid==35)+1)%-log2(nucData(:,smeta.tid==11)+1);
goodNucs=all(~isnan([xVals,yVals]),2);



figure;
subplot(1,2,1)
plot(smeta.tp(smeta.tid==39),diag(corr(xVals(nucMeta.gene>0,:),yVals(nucMeta.gene>0,:),'rows','pairwise')),'Linewidth',2)
xlabel('time in HU')
ylabel('corr(wt-rtt-TO)')
subplot(1,2,2)
hold off
linortfit2(mean(xVals(nucMeta.ToR<23,1),2,'omitnan'),mean(yVals(nucMeta.ToR<23,1),2,'omitnan'),[],nucMeta.ToR(nucMeta.ToR<23),'filled')
hold on
scatter(mean(xVals(ismember(nucMeta.gene,GP.groups{23}{2}{45})&nucMeta.order<0,2:3),2,'omitnan'),mean(yVals(ismember(nucMeta.gene,GP.groups{23}{2}{45})&nucMeta.order<0,2:3),2,'omitnan'),'filled')

asyncHa=find(smeta.tp==-10&contains(smeta.ab,'ha'))
[~,asyncMyc]=ismember(regexprep(smeta.name(asyncHa),'^ha','myc'),smeta.name);
asyncHa(asyncMyc==0)=[];
asyncMyc(asyncMyc==0)=[];


asyncTo=log2(nucData(:,asyncMyc)+0.1)-log2(nucData(:,asyncHa)+0.1);
plot(corr(nucMeta.ToR,asyncTo,'rows','pairwise'))
hold on
text(1:numel(asyncHa),corr(nucMeta.ToR,asyncTo,'rows','pairwise'),smeta.name(asyncHa))
ylabel('corr(asyncTO vs. ToR)')
axis tight
xlabel('samples')

%% cluster option 1
dnaCluster=load('DNAcluster.mat');
[~,nucMeta.binIdx]=min(abs(nucMeta.pos-dnaCluster.meta(:,1)'),[],2);
nucMeta.binIdx=dnaCluster.finalClu(nucMeta.binIdx);
nucMeta.binIdx(ismember(nucMeta.binIdx,find(accumarray(nucMeta.binIdx(nucMeta.binIdx>0),1)<100)))=NaN;
clear dnaCluster

%% absolute levels from clusters
nClu=max(nucMeta.binIdx);
cluSize=accumarray(nucMeta.binIdx(nucMeta.binIdx>0),1,[max(nucMeta.binIdx),1],@sum,0)
goodHa=mean(nucData(:,contains(smeta.ab,'ha')&~smeta.bad)>5,2)>0.9; % above 5 in 90% of the nucs
cluL=log2(cell2mat(arrayfun(@(x)median(nucData(nucMeta.binIdx==x&goodHa,:),'omitnan'),[1:nClu]','UniformOutput',false)))-log2(median(nucData(goodHa&nucMeta.binIdx>0,:),'omitnan'));
cluA=cell2mat(arrayfun(@(x)median(nucData(nucMeta.binIdx==x&goodHa,:),'omitnan'),[1:nClu]','UniformOutput',false));
%cluD=cell2mat(arrayfun(@(x)median(ducData(nucMeta.binIdx==x&goodHa,:),'omitnan'),[1:nClu]','UniformOutput',false));
cluTor=accumarray(nucMeta.binIdx(nucMeta.binIdx>0&goodHa),nucMeta.ToR(nucMeta.binIdx>0&goodHa),[max(nucMeta.binIdx),1],@(x)mean(x,'omitnan'),nan);

dTot=zeros(1,size(cluL,2))
binForTot={[1;3],[nClu+[-1:0]]}
sTot=zeros(1,size(cluL,2))
close all
c=0;
for g=find(strcmp(tcs.ab,'ha')&accumarray(smeta.tid,smeta.tp,[],@max)>5)'
    idx=find(smeta.tid==g&smeta.tp>-5 & ~smeta.bad);      
    [~,binForTot{2}]=mink(sum(cluL(:,idx),2),3);
    intMeans=[0 mean(cluL(binForTot{1},idx),1);0 mean(cluL(binForTot{2},idx),1)];
    dMeans=diff(intMeans,1,2)';
    [~,turningPoint]=min(movsum(dMeans(:),size(dMeans,1),'Endpoints','discard'));
     dTot(1,idx)=cumsum([max(-dMeans([numel(idx)+1:numel(idx)+turningPoint-1,turningPoint:numel(idx)]),0)]);       
    [temp,~,~]=fitAbsHA([0,smeta.tp(idx)'],[ones(1,1) 2.^dTot(:,idx)],'figure',false);
    sTot(1,idx)=log2(temp(2:end));
end
smeta.dTot=dTot';
clearvars dIdx intMeans dMeans turningPoint idx g c ans binForTot
absClu=2.^(cluL+smeta.dTot');
figure
hold off
for t=find(contains(tcs.gt,{'wt','rtt','mec'})&contains(tcs.ab,'ha')&~contains(tcs.exp,{'30','s35','38'}))'
    selSmp=find(smeta.tid==t)
    subplot(1,2,1)
    plot(smeta.tp(selSmp),sum(cluSize.*absClu(:,selSmp),'omitnan')./sum(cluSize),'linewidth',2,'DisplayName',strjoin(tcs{t,:}))
    hold on
end
legend()

%% calculate replicationr rate ine ach bin 
absCluF=NaN(nClu,size(smeta,1));
absCludT=NaN(nClu,size(smeta,1));
fit=0
for g=find(strcmp(tcs.ab,'ha')&accumarray(smeta.tid,smeta.tp,[],@max)>5)'
    selSmp=find(smeta.tid==g&smeta.tp>-5&~smeta.bad);    
    if fit==0
        temp=diff([ones(nClu,1) absClu(:,selSmp)],1,2)./diff([0;smeta.tp(selSmp)])'*100;
        absCludT(:,selSmp)=temp/2+temp(:,[2:end,end])/2; 
        absCluF(:,selSmp)=absClu(:,selSmp);
    else
            [temp,tempdT,~]=fitAbsHA([0,smeta.tp(selSmp)'],[ones(nClu,1) absClu(:,selSmp)],'figure',false);
            absCludT(:,selSmp)=tempdT(:,2:end)*100;
            absCluF(:,selSmp)=temp(:,2:end);
    end    
end
clear selSmp temp tempdT

%% cmpare TO
allHa=find(contains(smeta.ab,'ha'));
[~,allMyc]=ismember(regexprep(smeta.name(allHa),{'^.*?-','\*'},{'myc-',''}),smeta.name)
allHa=allHa(allMyc>0);
allMyc=allMyc(allMyc>0);
nucTo(:,allHa)=log2(nucData(:,allMyc)+1)-log2(nucData(:,allHa)+1);
for b=1:max(nucMeta.binIdx)
    medTo(b,:)=median(nucTo(nucMeta.binIdx==b,:),'omitnan');
end

%% Figure 2B
close all
hisNucs=find(ismember(nucMeta.gene,GP.groups{23}{2}{45})&nucMeta.order<0);
for t=13
    figure;
    selSmp=find(smeta.tid==t&smeta.tp>10&smeta.bad==0)
    [~,sMyc]=ismember(regexprep(smeta.name(selSmp),{'^ha','*'},{'myc',''}),smeta.name);
    selSmp=selSmp(sMyc>0);
    sMyc=sMyc(sMyc>0);
    c=0;
    c=c+1;
    subplot(3,5,c)
    legend()
    for s=selSmp'        
        sMyc=find(contains(smeta.name,regexprep(smeta.name{s},{'^ha','*'},{'myc',''})));
        c=c+1;
        p=robustfit(absCludT(50:end,s),log2(cluA(50:end,sMyc))-log2(cluA(50:end,s)));
        [~,topBin]=max(absCludT(:,s))
        selBins=cluA(:,s)<cluA(topBin,s);
        adjFit(s,:)=robustfit(absCludT(selBins,s),cluA(selBins,sMyc)-cluA(selBins,s).*2.^p(1));
        subplot(3,5,c)
        scatter(absCludT(:,s),cluA(:,sMyc)-cluA(:,s).*2.^p(1),[],cluTor,'filled')
        hold on
        plot(xlim,xlim.*adjFit(s,2)+adjFit(s,1),'k-')
        title(sprintf('%.2fx%+.2f; %.2f',adjFit(s,2),adjFit(s,1),mean(nucData(hisNucs([4,5,8]),sMyc)./absClu(nucMeta.binIdx(hisNucs([4,5,8])),s)) ))
        sgtitle(strjoin(tcs{t,:}))
    end    
end
figure
for t=[13,14,16]
    subplot(2,3,1)
    selSmp=find(smeta.tid==t&smeta.tp>10&smeta.bad==0)
    [~,sMyc]=ismember(regexprep(smeta.name(selSmp),{'^ha','*'},{'myc',''}),smeta.name);
    selSmp=selSmp(sMyc>0);
    sMyc=sMyc(sMyc>0);
    plot(smeta.tp(selSmp),diag(corr(absClu(:,selSmp),cluL(:,sMyc)-cluL(:,selSmp),'rows','pairwise')),'Linewidth',2,'DisplayName','cr with Level')
    hold on
    plot(smeta.tp(selSmp),diag(corr(absCludT(:,selSmp),cluL(:,sMyc)-cluL(:,selSmp),'rows','pairwise')),'Linewidth',2,'DisplayName','cr with RepRate')
    xlabel('time after HU')
    ylabel('correlation w. turnover')
    legend()
    subplot(2,3,2)
    plot(smeta.tp(selSmp),adjFit(selSmp,2),'Linewidth',2)
    hold on
    xlabel('time after HU')
    ylabel('myc-rep-slope')
    subplot(2,3,3)
    plot(smeta.tp(selSmp),sum(max(absCludT(:,selSmp).*adjFit(selSmp,2)',0).*cluSize,'omitnan')./sum(nucData(nucMeta.binIdx>0,sMyc),'omitnan'),'linewidth',2)
    hold on
end

dNuc=nan(size(nucData));
for t=find(contains(tcs.ab,{'9ac','k56ac','myc','ha'})&contains(tcs.exp,{'35','36','39','43','44'}))'
    selSmp=find(smeta.tid==t);
    alphaSmp=find(smeta.tid==t&smeta.tp==-5)
    for s=selSmp'
        [~,r]=robustfit(log2(mean(nucData(goodHa,alphaSmp),2)+1),log2(mean(nucData(goodHa,s),2)+1));
        dNuc(goodHa,s)=r.resid;
    end
end
close all
t1=37
for t2=[13,24,27]
    plot(smeta.tp(smeta.tid==t1),diag(corr(movmean(dNuc(:,smeta.tid==t1),25),movmean(dNuc(:,smeta.tid==t2),25),'rows','pairwise')),'linewidth',2)
    hold on
end

figure
for t=[8,10,11,13,14,16,12]
    selSmp=find(smeta.tid==t&smeta.tp>-10);
    [~,selMycs]=ismember(regexprep(smeta.name(selSmp),{'*','^ha'},{'','myc'}),smeta.name);
    subplot(4,1,1)
    plot(smeta.tp(selSmp),adjFit(selSmp,2),'Linewidth',2,'DisplayName',strjoin(tcs{t,:}))
    hold on
    subplot(4,1,2)
    plot(smeta.tp(selSmp),mean(nucData(hisNucs([8]),selMycs),1,'omitnan'),'Linewidth',2,'DisplayName',strjoin(tcs{t,:}))
    hold on
    subplot(4,1,3)
    plot(smeta.tp(selSmp),adjFit(selSmp,2)'./mean(nucData(hisNucs([8]),selMycs),1,'omitnan'),'Linewidth',2,'DisplayName',strjoin(tcs{t,:}))   
    hold on
    subplot(4,1,4)
    plot(smeta.tp(selSmp),sum(max(absCludT(:,selSmp).*adjFit(selSmp,2)',0).*cluSize,'omitnan')./sum(nucData(binIdx>0,selMycs),'omitnan'))   
    hold on
end

%% FIgure 3
%% 3A - TO, K56ac, K9ac
ha=find(contains(smeta.ab,'ha')&contains(smeta.gt,'wt')&contains(smeta.tag2,'h3x2')&contains(smeta.exp,'43')&smeta.tp==-5)
myc=find(contains(smeta.ab,'myc')&contains(smeta.gt,'wt')&contains(smeta.tag2,'h3x2')&contains(smeta.exp,'43')&smeta.tp==-5)
h3k56=find(contains(smeta.ab,'56')&contains(smeta.gt,'wt')&contains(smeta.tag2,'h3x2')&smeta.tp==-5)
h3k9=find(contains(smeta.ab,'k9')&contains(smeta.gt,'wt')&contains(smeta.tag2,'h3x2')&smeta.tp==-5)

haVal=movmean(mean(nucData(:,ha),2),1)
mycVal=movmean(mean(nucData(:,myc),2),1)
k56Val=movmean(mean(nucData(:,h3k56),2),1)
k9Val=movmean(mean(nucData(:,h3k9),2),1)

gNucs=(haVal>1)&(mycVal>1)&(k56Val>1);

subplot(1,2,1)
dscatter(log2(k56Val(gNucs)+1)-log2(haVal(gNucs)+1),log2(mycVal(gNucs)+1)-log2(haVal(gNucs)+1))%,[],log2(k9Val+1)-log2(haVal+1),'filled','Linewidth',0.5,'MarkerEdgeColor',[0 0 0])
xlabel('k56 enr')
ylabel('TO')
%ylabel(colorbar(),'k9 enr.')
ylabel(colorbar(),'data density')
title(sprintf('K56ac/HA and myc//HA in alpha cells (r=%.2f)',corr(log2(k56Val+1)-log2(haVal+1),log2(mycVal+1)-log2(haVal+1),'rows','pairwise')))
colormap(gca,flipud(brewermap(128,'purples')))
p=linortfit2(log2(k56Val(gNucs)+1)-log2(haVal(gNucs)+1),log2(mycVal(gNucs)+1)-log2(haVal(gNucs)+1));
hold on
plot([-1 2],[-1 2]*p(1)+p(2),'k--')
axis tight

set(get(gca,'Children'),'Marker','o')

save_gf(gcf,sprintf('Fig2TOvsK56inG1'),'type',{'fig'},'paper','gilad22','size',[])

%
intK56ac=find(contains(tcs.ab,'56')&contains(tcs.gt,'wt')&contains(tcs.exp,'s43'))
intMyc=find(contains(tcs.ab,'myc')&contains(tcs.gt,'wt')&contains(tcs.exp,'43'))
intHa=find(contains(tcs.ab,'ha')&contains(tcs.gt,'wt')&contains(tcs.exp,'43'))
intK9ac=find(contains(tcs.ab,'k9ac')&contains(tcs.gt,'wt')&contains(tcs.exp,'s43'))

xVals1=log2(movmean(nucData(:,smeta.tid==intK56ac),3)+.1)-log2(movmean(nucData(:,smeta.tid==intHa),3)+.1);
xVals2=movmean(dNucs(:,smeta.tid==intK9ac),3)-movmean(dNucs(:,smeta.tid==intHa),3);
xVals3=log2(movmean(nucData(:,smeta.tid==intK9ac),3)+.1)-log2(movmean(nucData(:,smeta.tid==intHa),3)+.1);
xVals4=movmean(dNucs(:,smeta.tid==intK56ac),3)-movmean(dNucs(:,smeta.tid==intHa),3);

yVals=log2(movmean(nucData(:,smeta.tid==intMyc),3)+.1)-log2(movmean(nucData(:,smeta.tid==intHa),3)+.1);

subplot(1,2,2)
hold off
plot(smeta.tp(smeta.tid==intHa),diag(corr(xVals1,yVals,'rows','pairwise')),'linewidth',2,'displayname','k56enr-TO')
hold on
plot(smeta.tp(smeta.tid==intHa),diag(corr(xVals3,yVals,'rows','pairwise')),'linewidth',2,'displayname','k9enr-TO')
plot(smeta.tp(smeta.tid==intHa),diag(corr(xVals2,yVals,'rows','pairwise')),'linewidth',2,'displayname','dk9enr-TO')
plot(smeta.tp(smeta.tid==intHa),diag(corr(xVals4,yVals,'rows','pairwise')),'linewidth',2,'displayname','dk56enr-TO')

axis tight
ylabel('correlation')
xlabel('time after hU')
%%
rttMyc=find(contains(smeta.ab,'myc')&contains(smeta.gt,'rtt')&smeta.tp==-5)
wtMyc=find(contains(smeta.ab,'myc')&contains(smeta.gt,'wt')&smeta.tp==-5&contains(smeta.tag2,'h3x2'))

figure
scatter(log2(1+mean(nucData(:,wtMyc),2)),log2(1+mean(nucData(:,rttMyc),2)),[],nucData(:,contains(smeta.ab,'k56')&smeta.tp==-5&contains(smeta.gt,'wt')&contains(smeta.exp,'s43')),'filled')
xlabel(smeta.name(wtMyc(1)))
ylabel(smeta.name(rttMyc(1)))
ylabel(colorbar(),'k56ac')
% turnover of replicating (for cluster and ind. nucleoosmes)
for t=find(contains(tcs.ab,'ha')&contains(tcs.gt,{'rtt','wt'}))'
    selSmp=find(smeta.tid==t);
    if numel(selSmp)>5
        figure;
        c=0;
        for s=selSmp'
            c=c+1;
            subplot(3,5,c)
            scatter(absCludT(:,s),medTo(:,s),[],absClu(:,s),'filled')
        end
    end
end
close all
for t=find(contains(tcs.ab,'ha')&contains(tcs.gt,{'rtt','wt','vps'}))'
    selSmp=find(smeta.tid==t);
    if numel(selSmp)>5
        figure;
        c=0;
        for s=selSmp'
            mycS=allMyc(allHa==s);
            if numel(mycS)>0
                c=c+1;
                subplot(3,5,c)
                scatter(absCludT(:,s),diff(log2(cluA(:,[s mycS])),1,2),[],absClu(:,s),'filled')
                xlabel('repRate')
                ylabel('log2(Myc_c/HA_c)')
                title(extractBetween(smeta.name{s},'ha-','.mat'))
                caxis([1 1.7])
            end
        end
        ylabel(colorbar(),'DNA content')
    end
end

% similiar turnover of not replicating
GP=load('/home/labs/barkailab/felixj/Documents/MATLAB/scripts/gene/group_imp.mat');
intExp={'x35','s39'}
hisNucs=ismember(nucMeta.gene,GP.groups{23}{2}{45})&nucMeta.order<0;
close all
repTh=0.1
cMap=brewermap(9,'purples');
for e=1:numel(intExp)
    figure('Units','normalized','Position',[0 0 0.5 1], 'color','w')
    c=0;    
    alphaMyc=log2(mean(nucData(:,(contains(smeta.gt,{'wt','rtt'})&contains(smeta.ab,'myc')&smeta.tp==-5&contains(smeta.exp,intExp{e}))),2)+1);
    for tp=unique(smeta.tp(contains(smeta.exp,intExp{e})&smeta.tp>60))'
        c=c+1;
        xMyc=find(contains(smeta.gt,'wt')&contains(smeta.ab,'myc')&smeta.tp==tp&contains(smeta.exp,intExp{e}),1);
        xHa=find(contains(smeta.gt,'wt')&contains(smeta.ab,'ha')&smeta.tp==tp&contains(smeta.exp,intExp{e}),1);        
        yMyc=find(contains(smeta.gt,'rtt')&contains(smeta.ab,'myc')&smeta.tp==tp&contains(smeta.exp,intExp{e}),1);
        yHa=find(contains(smeta.gt,'rtt')&contains(smeta.ab,'ha')&smeta.tp==tp&contains(smeta.exp,intExp{e}),1);
        if numel([xMyc,yMyc,xHa,yHa])==4
            xVals=log2(nucData(:,[xMyc])+1);%diff((dNucs(:,[xHa xMyc])),1,2);%
            yVals=log2(nucData(:,[yMyc])+1);%diff((dNucs(:,[yHa yMyc])),1,2);%
            selBin=find(all(absCludT(:,[xHa,yHa])<repTh,2)&all(absClu(:,[xHa,yHa])<1.15,2));
            %repLvl=nan(size(xVals,1),2);
            %repLvl(nucMeta.binIdx>0,1:2)=absClu(nucMeta.binIdx(nucMeta.binIdx>0),[xHa yHa]);
            subplot(3,3,c)
            selNucs=ismember(nucMeta.binIdx,selBin)&all(~isnan([xVals yVals]),2)&~hisNucs;            
            scatter(xVals(selNucs),yVals(selNucs),20,.8.*[1 1 1],'filled','MarkerFaceAlpha',0.5)
            hold on
            scatter(xVals(selNucs&alphaMyc>4),yVals(selNucs&alphaMyc>4),[],cMap(end-4,:),'filled','DisplayName','noRep+high','MarkerEdgeColor',cMap(end,:))
            %subplot(4,6,c)
            p=linortfit2(xVals(selNucs&alphaMyc>4),yVals(selNucs&alphaMyc>4));
            %scatter(xVals(hisNucs),yVals(hisNucs),'r.')
            axis tight
            title(sprintf('Fit eq: %.2fx%+.2f,%d',p(1),p(2),sum(selNucs)))
            hold on
            plot(xlim,xlim,'k-','DisplayName','1:1')
            plot(xlim,xlim*p(1)+p(2),'k--','DisplayName','linFit')
            axis tight
            xlabel(sprintf('%s-log2',smeta.name{xMyc}))
            ylabel(sprintf('%s-log2',smeta.name{yMyc}))
            legend('show','Location','best')
        end
    end
end

alphaMyc=log2(mean(nucData(:,(contains(smeta.gt,{'wt','rtt','vps'})&contains(smeta.ab,'myc')&smeta.tp==-5&contains(smeta.exp,goodExps))),2)+1);
alphaHa=log2(mean(nucData(:,(contains(smeta.gt,{'wt','rtt','vps'})&contains(smeta.ab,'ha')&smeta.tp==-5&contains(smeta.exp,goodExps))),2)+1);

endLvl=mean(absClu(:,contains(smeta.ab,'ha')&contains(smeta.exp,goodExps)&smeta.tp>150&contains(smeta.exp,goodExps)&contains(smeta.gt,{'vps','rtt','wt'})),2)
selBin=find(endLvl<1.11)
intTps=unique(smeta.tp(contains(smeta.exp,goodExps)&contains(smeta.gt,{'wt','rtt'}) &smeta.tp>=70))
cmpGt={'rtt','vps'}
mycTh=4.5;
close all
for g=1:numel(cmpGt)
    figure;
    c=0;
    for tp=intTps'
        c=c+1;
        xMyc=find(contains(smeta.gt,'wt')&contains(smeta.ab,'myc')&smeta.tp==tp&contains(smeta.exp,goodExps)&smeta.bad==0);
        xHa=find(contains(smeta.gt,'wt')&contains(smeta.ab,'ha')&smeta.tp==tp&contains(smeta.exp,goodExps)&smeta.bad==0);
        yMyc=find(contains(smeta.gt,cmpGt{g})&contains(smeta.ab,'myc')&smeta.tp==tp&contains(smeta.exp,goodExps)&smeta.bad==0);
        yHa=find(contains(smeta.gt,cmpGt{g})&contains(smeta.ab,'ha')&smeta.tp==tp&contains(smeta.exp,goodExps)&smeta.bad==0);
        xVals=log2(mean(nucData(:,[xMyc]),2)+1)-log2(mean(nucData(:,[xHa]),2)+1);%diff((dNucs(:,[xHa xMyc])),1,2);%
        yVals=log2(mean(nucData(:,[yMyc]),2)+1)-log2(mean(nucData(:,[yHa]),2)+1);%diff((dNucs(:,[yHa yMyc])),1,2);%
        %xVals=2.^xVals;
        %yVals=2.^yVals;
        %repLvl=nan(size(xVals,1),2);
        %repLvl(nucMeta.binIdx>0,1:2)=absClu(nucMeta.binIdx(nucMeta.binIdx>0),[xHa yHa]);
        subplot(2,3,c)
        hold off
        selNucs=ismember(nucMeta.binIdx,selBin)&all(~isnan([xVals yVals]),2)&~hisNucs;
        scatter(xVals(selNucs),yVals(selNucs),20,.8.*[1 1 1],'filled','MarkerFaceAlpha',0.5)
        hold on
        hlNucs=selNucs&((alphaMyc-alphaHa)>1.5);
        scatter(xVals(hlNucs),yVals(hlNucs),[],cMap(end-4,:),'filled','DisplayName','noRep+high','MarkerEdgeColor',cMap(end,:))
        %subplot(4,6,c)
        p=linortfit2(xVals(hlNucs),yVals(hlNucs));
        %scatter(xVals(hisNucs),yVals(hisNucs),'r.')
        axis tight
        title(sprintf('Fit eq: %.2fx%+.2f,%d',p(1),p(2),sum(selNucs)))
        hold on
        plot(xlim,xlim,'k-','DisplayName','1:1')
        plot(xlim,xlim*p(1)+p(2),'k--','DisplayName','linFit')
        axis tight
        xlabel(sprintf('%s-log2',smeta.name{xMyc(1)}))
        ylabel(sprintf('%s-log2',smeta.name{yMyc(1)}))
        legend('show','Location','best')
        xlim([-2 4])
        ylim([-2 4])
        %xlim([-3 6])
        %ylim([-3 6])
    end
    %save_gf(gcf,sprintf('%s_nonRepCmp',cmpGt{g}),'type',{'fig'},'paper','gilad22','size',[])
end

%% compare all Nucs rtt109 wt
selTP=-5;

xMyc=find(contains(smeta.gt,'wt')&smeta.tp==selTP&contains(smeta.ab,'myc')&ismember(smeta.exp,setdiff(goodExps,'c51.mat'))&smeta.bad==0)
xHa=find(contains(smeta.gt,'wt')&smeta.tp==selTP&contains(smeta.ab,'ha')&ismember(smeta.exp,setdiff(goodExps,'c51.mat'))&smeta.bad==0)
xVal=log2(mean(nucData(:,xMyc),2)+1)-log2(mean(nucData(:,xHa),2)+1);

yMyc=find(contains(smeta.gt,'rtt')&smeta.tp==selTP&contains(smeta.ab,'myc')&ismember(smeta.exp,setdiff(goodExps,'c51.mat'))&smeta.bad==0)
yHa=find(contains(smeta.gt,'rtt')&smeta.tp==selTP&contains(smeta.ab,'ha')&ismember(smeta.exp,setdiff(goodExps,'c51.mat'))&smeta.bad==0)
yVal=log2(mean(nucData(:,yMyc),2)+1)-log2(mean(nucData(:,yHa),2)+1);
goodNucs=all(nucData(:,[xHa])>1,2)
figure
subplot(1,2,1)
hold off
%dscatter(xVal(goodNucs),yVal(goodNucs))
%set(get(gca,'Children'),'Marker','o')
%colormap(gca,flipud(brewermap(128,'buPu')))
%ylabel(colorbar(),'Data density')

scatter(xVal(goodNucs),yVal(goodNucs),[],yVal(goodNucs)-xVal(goodNucs),'filled')
colormap([flipud(brewermap(128,'purples'));brewermap(128,'blues')])
caxis(1.5*[-1 1])
ylabel(colorbar(),'\Delta Turnover')

% scatter(xVal(goodNucs),yVal(goodNucs),[],mean([yVal(goodNucs),xVal(goodNucs)],2,'omitnan'),'filled')
% colormap(gca,(brewermap(128,'buPu')))
% caxis([-1 4.5])
% ylabel(colorbar(),'average Turnover')

hold on

highNucs=(xVal>3)&all(nucData(:,[xHa])>1,2)
scatter(xVal(highNucs),yVal(highNucs),[],[0 0 0],'DisplayName','high turnover')


xlabel(sprintf('wt turnover',selTP))
ylabel(sprintf('rtt109 turnover',selTP))
title(sprintf('%d %.2f',selTP,corr(xVal,yVal,'rows','pairwise')))
plot(xlim,xlim,'k--')

subplot(1,3,2)
dscatter(log2(1+mean(nucData(goodNucs,xMyc),2,'omitnan')),log2(1+mean(nucData(goodNucs,yMyc),2,'omitnan')))
set(get(gca,'Children'),'Marker','o')
colormap(gca,flipud(brewermap(128,'buPu')))
xlabel(sprintf('wt Myc at %d',selTP))
ylabel(sprintf('rtt109 Myc at %d',selTP))
hold on
plot(xlim,xlim,'k--')
subplot(1,3,3)
dscatter(log2(1+mean(nucData(goodNucs,xHa),2,'omitnan')),log2(1+mean(nucData(goodNucs,yHa),2,'omitnan')))
set(get(gca,'Children'),'Marker','o')
colormap(gca,flipud(brewermap(128,'buPu')))
xlabel(sprintf('wt HA at %d',selTP))
ylabel(sprintf('rtt109 HA at %d',selTP))
hold on
plot(xlim,xlim,'k--')
sgtitle(sprintf('%d',selTP))

save_gf(gcf,sprintf('FigS3rttVswtaFactor'),'type',{'fig'},'paper','gilad22','size',[])
intGts=unique(sampleTable.gt)'

plotVlines(nucData,smeta,nucMeta,[xMyc;yMyc])

close all
for g={'wt','vps75','rtt'}
    figure
    c=0
for ab ={'ha','myc','k9ac','k56ac'}
    for tp=[40,50,70,90,110,130,150]
        c=c+1;
        subplot(4,7,c)
        hold off
        selSmp=find(contains(smeta.ab,ab)&ismember(smeta.gt,g)&smeta.tp==tp&ismember(smeta.exp,goodExps))
        meanLvl=mean(nucData(:,selSmp),2);
        selDna=find(contains(smeta.ab,'DNA')&ismember(smeta.gt,g)&smeta.tp==tp&ismember(smeta.exp,goodExps))
        meanDna=mean(nucData(:,selDna),2);
        
        selB=nucMeta.binIdx>0&nucMeta.order>0;
        scatter(log2(accumarray(nucMeta.binIdx(selB),meanDna(selB),[],@(x)median(x,'omitnan'))),log2(accumarray(nucMeta.binIdx(selB),meanLvl(selB),[],@(x)median(x,'omitnan'))),'filled')
        %        scatter((accumarray(nucMeta.binIdx(selB),meanDna(selB),[],@(x)median(x,'omitnan'))),(accumarray(nucMeta.binIdx(selB),meanLvl(selB),[],@(x)median(x,'omitnan'))),'filled')

        hold on
        
        selB=nucMeta.binIdx>0&nucMeta.order<0;
        scatter(log2(accumarray(nucMeta.binIdx(selB),meanDna(selB),[],@(x)median(x,'omitnan'))),log2(accumarray(nucMeta.binIdx(selB),meanLvl(selB),[],@(x)median(x,'omitnan'))),'filled')
        %scatter((accumarray(nucMeta.binIdx(selB),meanDna(selB),[],@(x)median(x,'omitnan'))),(accumarray(nucMeta.binIdx(selB),meanLvl(selB),[],@(x)median(x,'omitnan'))),'filled')
    end
end
end