function figure4_ver2()
clear all
%% load nucleosome data
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/nucMeta015AddFried.mat')
[data,smeta]=getNaamaData('figure',4);
h3B={'fast1b-wt-async-x18','fast2b-wt-async-x18','fast2a-wt-async-x22','fast2b-wt-async-x22'};
smeta.tag2(contains(smeta.name,h3B))={'hb3'};
smeta.con(strcmp(smeta.con,'async') & strcmp(smeta.gt,'wt') & contains(smeta.exp,'x22') & contains(smeta.tag,'fast1'),:)={'alpha0hr'}
nBase=size(data,1);
nucExt=50;
nucPos=zeros(nBase,1);
for i=1:size(nucMeta,1)
    posi=nucMeta.pos(i)-nucExt:nucMeta.pos(i)+nucExt;
    if i>1
        posi=posi(posi>0 & posi>mean(nucMeta.pos([i-1,i])) & posi<=nBase);
    end
    if i<size(nucMeta,1)
        posi=posi(posi<mean(nucMeta.pos([i+1,i])) & posi>0);
    end
    nucPos(posi)=i;
end
nucData=zeros(size(nucMeta,1),size(data,2));
for i=1:size(nucData,2)
    i
    nucData(:,i)=accumarray(nucPos(nucPos>0),data(nucPos>0,i),[],@mean);
end
clear data i posi h3B nucPos nucExt
highOcc=find(any(movmean(nucData(:,strcmp(smeta.ab,'ha')& strcmp(smeta.gt,'wt')),5)>30,2));
nucData(unique(acol(highOcc+[-2:2])),:)=NaN;
lowOcc=find(any(movmean(nucData(:,strcmp(smeta.ab,'ha') & strcmp(smeta.gt,'wt')),3)<0.01,2));
nucData(unique(acol(lowOcc+[-1:1])),:)=NaN;
clearvars  highOcc lowOcc
[smeta,idx]=sortrows(smeta,[7,14,15,9,10]);
nucData=nucData(:,idx);
clear idx
%% check deleltions
GP=load('group_imp.mat');
loci=movmean(ismember(nucMeta.gene,[GP.gene_table.RTT109,GP.gene_table.HIR1,GP.gene_table.ASF1,GP.gene_table.BAR1]),2)>0.25;
figure;imagesc(nucData(loci,contains(smeta.exp,'x26')))
% imagesc(corr(nucData,'rows','pairwise'))
% yticks(1:size(smeta,1))
% yticklabels(smeta.name) 
%% load additional data
regev=load ('/home/labs/barkailab/LAB/data/DataExternal/EXPRESSION/regev_pnas_2008_rnaseq/data');
geneLvl=median(regev.regev_rnaseq,2,'omitnan');
load('holstege.mat')
geneStdH=std(holstege.TsXMut,[],2,'omitnan');
meanChg=median(holstege.TsXMut,2,'omitnan');

dirStd=geneStdH.*(2*(meanChg>0)-1);
dirStd(abs(meanChg)<2*geneStdH./sqrt(sum(~isnan(holstege.TsXMut),2)))=NaN;

GP=load('group_imp.mat');
ToR = readtable('/home/labs/barkailab/yuliago/Documents/MATLAB/ToR_raz.csv');
ToR.idx=GP.chrIdx(ToR.chr)+ToR.start+500;
[val,idx]=min(abs(nucMeta.pos-ToR.idx'),[],2);
nucMeta.ToR(val<501)=ToR.score(idx(val<501));
nucMeta.ToR(nucMeta.ToR==0)=NaN;
load('spell.mat')
geneStdS=std(spell.data(:,spell.abs==1),[],2);
clear holstege regev ToR idx val spell
[expTypes,~,smeta.expId]=unique(smeta(:,[7,9,10,14]),'stable')
[~,~,smeta.sid]=unique(smeta(:,[8:11]),'stable')
expTypes.n=accumarray(smeta.expId,smeta.bad==0)
medNuc=nan(size(nucData,1),max(smeta.expId));
for i=unique(smeta.expId(smeta.bad==0 ))'   
    selRpt=find(smeta.expId==i & ~smeta.bad);
    selNucs=all(nucData(:,selRpt)>1,2);
    [~,bestRpt]=max(sum(corr(log(nucData(selNucs,selRpt))),2));
    bestRpt=selRpt(bestRpt);
    for j=selRpt'
        pj=median(diff(log(nucData(selNucs,[j bestRpt])),1,2));
        fac(j)=exp(pj);
    end
    medNuc(:,i)=median(nucData(:,selRpt).*fac(selRpt),2); 
    %stdNuc(:,i)=std(log(nucData(:,selRpt).*fac(selRpt)+1),[],2);  
end
clear crMat selRpt i j fac bestRpt pj selNucs loci 
save('forFigure4add.mat')

clear all
load('forFigure4add.mat')
close all

%% scatter plots
intNuc=[-2,-1,2]
intExp=find(strcmp(expTypes.ab,'myc')&strcmp(expTypes.gt,'wt')&strcmp(expTypes.con,'async') & startsWith(expTypes.tag2,{'h2a','h3'}))
colorLim=[2 17;3 18]
figure
for n=1:numel(intNuc)
    selNucs=ismember(nucMeta.order,intNuc(n))& all(~isnan(medNuc(:,intExp)),2) & ismember(nucMeta.gene,find(all(~isnan([geneLvl,geneStdH]),2))) & ~ismember(nucMeta.gene,GP.groups{23}{2}{45});
    cData=medNuc(selNucs,intExp);
    [cSmooth,cNum]=smooth2ParaPlot([geneLvl(nucMeta.gene(selNucs)) geneStdH(nucMeta.gene(selNucs))],cData,'posLim',[4,0;16,0.7])
    for e=1:numel(intExp)
        subplot(numel(intExp),numel(intNuc),(e-1)*numel(intNuc)+n)
        scatter(geneLvl(nucMeta.gene(selNucs)),geneStdH(nucMeta.gene(selNucs)),[],cSmooth(:,e),'.')
        caxis(colorLim(e,:))
        set(gca,'Color',.7*[1 1 1])        
        ylabel('geneStdH')      
        if n==1
            colorbar('west')
        end
        axis tight
        xlim([5 15])   
        text(5,max(ylim),sprintf('%s,nuc:%+d',strjoin(expTypes{intExp(e),[1:4]}),intNuc(n)),'HorizontalAlignment','left','VerticalAlignment','top')
        if e==numel(intExp)
            xlabel('geneLvl')
        else
            xticks([])
        end
    end      
end


%% compare geneLvl, gneweStdH and meanChg to define dirStd
figure;
selGenes=abs(meanChg)>=2*geneStdH./sqrt(sum(~isnan(holstege.TsXMut),2));
scatter(geneStdH(selGenes),meanChg(selGenes),[],dirStd(selGenes),'.','DisplayName','flexibel genes')
hold on
scatter(geneStdH(~selGenes),meanChg(~selGenes),[],[.7 .7 .7],'.','MarkerEdgeAlpha',.3,'DisplayName','not significant')

hold on;plot([fliplr(xlim) xlim],[fliplr(xlim)*2./sqrt(sum(~isnan(holstege.TsXMut(1,:)),2)) -xlim*2./sqrt(sum(~isnan(holstege.TsXMut(1,:)),2))],'k--','DisplayName','2X StdError')
%plot(xlim,xlim*2./sqrt(sum(~isnan(holstege.TsXMut),2)),'k--','DisplayName','2X StdError')
xlim([0 .25])
ylim([-.06 .06])
ylabel(colorbar(),'flexibility')
caxis([-0.2 +0.2])
xlabel('geneStdH')
ylabel('median Chg in mutants')
colormap(brewermap(64,'Spectral'))
save_gf(gcf,'Fig4_classifyGenes2')


%% vplot
load('holstege.mat')
% load('spell.mat')
% geneStdS=std(spell.data(:,spell.abs==1),[],2);
% geneChgS=mean(spell.data(:,spell.abs==1),2,'omitnan');

dirStd=geneStdH.*sign(meanChg);
dirStd(abs(meanChg)<2*geneStdH./sqrt(sum(~isnan(holstege.TsXMut),2)))=NaN;
%dirStd(abs(dirStd)>.3)=NaN;
dirStdClean=dirStd;
%dirStdClean(dirStd<quantile(dirStd,0.025) | dirStd>quantile(dirStd,0.975))=NaN;
dirStdClean(abs(dirStd)>quantile(abs(dirStd),.95))=NaN
%dirStd2=geneStdS.*(2*(meanChg>0)-1);
%dirStd2(abs(meanChg)<2*geneStdH./sqrt(sum(~isnan(holstege.TsXMut),2)))=NaN;

dirStdClean=-dirStdClean;
figure;
intNuc=[-2];
%intExp=find(contains(expTypes.ab,{'myc','ha'})&strcmp(expTypes.gt,'wt')&strcmp(expTypes.con,'async') & startsWith(expTypes.tag2,{'h3'}))
%intExp=find(contains(expTypes.ab,{'myc','ha'})&strcmp(expTypes.gt,'wt')&strcmp(expTypes.con,'async') & startsWith(expTypes.tag2,{'h2a'}))
%intExp=find(contains(expTypes.ab,{'myc','ha'})&strcmp(expTypes.gt,'hir1')&strcmp(expTypes.con,'async') & startsWith(expTypes.tag2,{'h2a'}))
intExp=[74,170,78,174]
%intExp=find(contains(expTypes.ab,{'myc','ha'})&strcmp(expTypes.gt,'hir1')&strcmp(expTypes.con,'async') & startsWith(expTypes.tag2,{'h3'}))
%intExp=find(contains(expTypes.ab,{'myc','ha'})&strcmp(expTypes.gt,'asf')&strcmp(expTypes.con,'async') & startsWith(expTypes.tag2,{'h2a'}))
intExp=[2,14]
intExp=[3,15]
intPar={'dirStdClean'}
movPar=[0.05]
seperateCr=[1,0];

intPar={'geneLvl'}
movPar=[0.25]
seperateCr=0
intNuc=[2]

%lineColor=[lines(2);0.3.*[1 1 1];0.6*[1 1 1]]

blackList=ismember(nucMeta.gene,GP.groups{23}{2}{45})% & ismember(nucMeta.gene,find(geneLvl>12));
step=5;
selPlot=1:12
%close all
%c=0;
if contains(smeta.tag2(intSmp(1)),{'2'})
    lineColor=[100,100,100;127,63,152]/256;
else
    lineColor=[100,100,100;0,174,239]/256;
end
for p=1:numel(intPar)
    parData=eval(intPar{p});
    for n=1:numel(intNuc)
        %subplot(numel(intPar),numel(intNuc),(p-1)*numel(intNuc)+n)
        c=c+1;
        subplot(2,2,c)
        hold off
        selNucs=find(nucMeta.order==intNuc(n) & ismember(nucMeta.gene,find(~isnan(parData))) &~blackList & all(~isnan(medNuc(:,intExp)),2) );
        hisNucs=find(nucMeta.order==intNuc(n) & ismember(nucMeta.gene,GP.groups{23}{2}{45}));
        [xVal,idx]=sort(parData(nucMeta.gene(selNucs)));
        selNucs=selNucs(idx);
        clear lineObj
        for e=1:numel(intExp) 
            expVal=medNuc(selNucs,intExp(e));
            %expVal=expVal./mean(expVal(abs(xVal-quantile(parData,0.05))<=movPar(p)/2));
            yVal=movmean(expVal,movPar(p),'omitnan','SamplePoints',makeUnique(xVal,4));
            nVal=movsum(~isnan(expVal),movPar(p),'omitnan','SamplePoints',makeUnique(xVal,4))
            yValE=movstd(expVal,movPar(p),'omitnan','SamplePoints',makeUnique(xVal,4))./sqrt(movsum(~isnan(expVal),movPar(p),'omitnan','SamplePoints',makeUnique(xVal,4)));
            if seperateCr(p)==1
                cr1=corr(xVal(xVal<=0),expVal(xVal<=0),'rows','pairwise');
                cr2=corr(xVal(xVal>=0),expVal(xVal>=0),'rows','pairwise');
                lineName=sprintf('%s,%.2f,%.2f',strjoin(expTypes{intExp(e),[2,1,4]}),cr1,cr2);
                idxN=find((xVal<max(xVal(xVal<0))-movPar(p)/2) & mod([1:numel(xVal)]',step)==1,1,'last');
                idxP=find(xVal>min(xVal(xVal>0))+movPar(p)/2,1,'first')+step;
                lineObj(e)=plot(xVal(1:step:idxN),yVal(1:step:idxN),'-','Linewidth',2,'Color',brighten(lineColor(e,:),0),'DisplayName',lineName);
                hold on                        
                fill([xVal(1:step:idxN);flipud(xVal(1:step:idxN))],[yVal(1:step:idxN)-yValE(1:step:idxN);flipud(yVal(1:step:idxN)+yValE(1:step:idxN))],brighten(lineColor(e,:),0),'LineStyle','none','FaceAlpha',.3)
                plot(xVal(idxP:step:end),yVal(idxP:step:end),'-','Linewidth',2,'Color',brighten(lineColor(e,:),0),'DisplayName',lineName)
                hold on                        
                fill([xVal(idxP:step:end);flipud(xVal(idxP:step:end))],[yVal(idxP:step:end)-yValE(idxP:step:end);flipud(yVal(idxP:step:end)+yValE(idxP:step:end))],brighten(lineColor(e,:),0),'LineStyle','none','FaceAlpha',.3)
                plot(xVal([idxN,idxP]),yVal([idxN,idxP]),'--','Linewidth',2,'Color',brighten(lineColor(e,:),0))
            else
                cr=corr(xVal,expVal,'rows','pairwise');
                lineName=sprintf('%s,%.2f',strjoin(expTypes{intExp(e),[2,1,4]}),cr);
                lineObj(e)=plot(xVal(1:step:end),yVal(1:step:end),'-','Linewidth',2,'Color',brighten(lineColor(e,:),0),'DisplayName',lineName)
                hold on                        
                fill([xVal(1:step:end);flipud(xVal(1:step:end))],[yVal(1:step:end)-yValE(1:step:end);flipud(yVal(1:step:end)+yValE(1:step:end))],brighten(lineColor(e,:),0),'LineStyle','none','FaceAlpha',.3)
            end            
            %scatter(parData(nucMeta.gene(hisNucs)),medNuc(hisNucs,intExp(e)),[],lineColor(e,:),'o','filled')
        end
        xlabel(intPar{p})
        axis tight
        %xlim([-.3 .3])%(quantile(parData,[-0.3 0.3]))
        %axis('tight')        
        title(intNuc(n))
        xlim(quantile(xVal,[0 1])+movPar(p).*[+.5 -.5])
        %plot([quantile(xVal(xVal<0),.05) quantile(xVal(xVal>0),.95)].*[1;1],ylim()'.*[1 1],'k-')
        %if n==1 & p==1
            legend(lineObj)
        %end
%         yyaxis right
%         hold off
%         plot(xVal(xVal<0),sum(xVal(xVal<0)'>(xVal+movPar(p)/2)),'k-')        
%         hold on
%         plot(xVal(xVal>0),sum(xVal(xVal>0)'<(xVal-movPar(p)/2)),'k-')
    end
    ylabel(sprintf('%.3f - movAcg',movPar(p)))
end
save_gf(gcf,sprintf('Fig4_%sLines','H2A'))
save_gf(gcf,sprintf('Fig4_%sVLines','All'))
intNuc=[-2 -1 1 2]
intExp=find(contains(expTypes.ab,{'ha','myc'})&contains(expTypes.con,{'async'})&startsWith(expTypes.tag2,{'h2a','h3'})&strcmp(expTypes.gt,'wt'))
parData=dirStd;
crMat=cell(2,1)
for n=1:numel(intNuc)
    selNucs=find(nucMeta.order==intNuc(n) & ismember(nucMeta.gene,find(~isnan(parData))) &~blackList & all(~isnan(medNuc(:,intExp)),2));
    parVal=parData(nucMeta.gene(selNucs));
    expVal=medNuc(selNucs,intExp);
    crMat{1}(n,:)=corr(parVal(parVal<0),expVal(parVal<0,:));
    crMat{2}(n,:)=corr(parVal(parVal>0),expVal(parVal>0,:));
end
crMat=cat(3,crMat{:})
subplot(2,4,3)
hold off
imagesc([squeeze(crMat(:,4,:)) squeeze(crMat(:,2,:))],[-.25 .25])
ylabel('nuc position')
yticks(1:numel(intNuc))
yticklabels(intNuc)
xticks([1,3])
xticklabels(expTypes.ab(intExp([3,1])))
colormap(gca,'redblue')
colorbar()
title(expTypes.tag2(intExp(4)))

subplot(2,4,7)
hold off
imagesc([squeeze(crMat(:,3,:)) squeeze(crMat(:,1,:))],[-.25 .25])
ylabel('nuc position')
yticks(1:numel(intNuc))
yticklabels(intNuc)
colormap(gca,'redblue')
xticks([1,3])
xticklabels(expTypes.ab(intExp([3,1])))
colorbar()
title(expTypes.tag2(intExp(3)))

%% histone plots
figure;
intExp=[
    find(contains(expTypes.ab,{'ha'})&contains(expTypes.con,{'async'})&startsWith(expTypes.tag2,{'h3'})&contains(expTypes.gt,{'asf','wt'}))';
    find(contains(expTypes.ab,{'myc'})&contains(expTypes.con,{'async'})&startsWith(expTypes.tag2,{'h3'})&contains(expTypes.gt,{'asf','wt'}))'
    find(contains(expTypes.ab,{'ha'})&contains(expTypes.con,{'async'})&startsWith(expTypes.tag2,{'h3'})&contains(expTypes.gt,{'hir1','wt'}))';
    find(contains(expTypes.ab,{'myc'})&contains(expTypes.con,{'async'})&startsWith(expTypes.tag2,{'h3'})&contains(expTypes.gt,{'hir1','wt'}))';    
    find(contains(expTypes.ab,{'ha'})&contains(expTypes.con,{'async'})&startsWith(expTypes.tag2,{'h2a'})&contains(expTypes.gt,{'hir1','wt'}))';
    find(contains(expTypes.ab,{'myc'})&contains(expTypes.con,{'async'})&startsWith(expTypes.tag2,{'h2a'})&contains(expTypes.gt,{'hir1','wt'}))';
    ];
hisPromNucs=ismember(nucMeta.gene,GP.groups{23}{2}{45}) & nucMeta.order<0;
hisORFNucs=ismember(nucMeta.gene,GP.groups{23}{2}{45}) & nucMeta.order>0;

for e=1:size(intExp,1)
    subplot(3,2,e)
    hold off
    xVal=(medNuc(:,intExp(e,2)));
    yVal=(medNuc(:,intExp(e,1)));
    xVal(medNuc(:,intExp(e,2))<0.2)=NaN;   
    yVal(medNuc(:,intExp(e,1))<0.2)=NaN;
    dscatter(xVal(nucMeta.order~=0 & ~isnan(xVal)& ~isnan(yVal)),yVal(nucMeta.order~=0& ~isnan(xVal)& ~isnan(yVal)))
    hold on
    scatter(xVal(hisPromNucs),yVal(hisPromNucs),'ro','filled','DisplayName','hisPromNucs')
    scatter(xVal(hisORFNucs),yVal(hisORFNucs),[],.5*[1 1 1],'o','filled','DisplayName','hisORFNucs')    
    xlabel(strjoin(expTypes{intExp(e,2),[1:4]}))
    ylabel(strjoin(expTypes{intExp(e,1),[1:4]}))
    axis tight
    text(xVal(hisPromNucs),yVal(hisPromNucs),GP.gene_infoR64.name(nucMeta.gene(hisPromNucs)))
end
save_gf(gcf,'Fig4_mutDensScatterWH2A')

%% profiels
clear all
[data,smeta]=getNaamaData('figure',4);
intData=contains(smeta.ab,{'myc','ha'}) & contains(smeta.gt,{'asf','wt','hir1','hir2','spt6oe'}) & contains(smeta.tag2,{'h3','h2a'})& ~contains(smeta.tag2,{'nc_h2a'})& contains(smeta.con,{'async'});
data=data(:,intData);
smeta=smeta(intData,:);
[~,idx]=sortrows(smeta(:,[7,14,9,10]));
smeta=smeta(idx,:);
data=data(:,idx);
clear intData ans
[expTypes,~,smeta.expId]=unique(smeta(:,[14,7,9,10]),'stable')
for i=unique(smeta.expId(smeta.bad==0))'
    selSmp=smeta.expId==i &smeta.bad==0;
    meanProfile(:,i)=mean(data(:,selSmp),2);
end
GP=load('group_imp.mat');
tssIdx=nan(6701,1);
hasTss=~isnan(GP.gene_infoR64.felixTss(:,1));
tssIdx(hasTss)=GP.chrIdx(GP.gene_infoR64.felixTss(hasTss,1))+GP.gene_infoR64.felixTss(hasTss,2);

intGene={[102,103],[221 222],[1146 1147],[4636 4637]}
close all
%intCmb={[36,26,28,24]}%{[14,10:12]}
%lineColor=lines(4);


%intCmb={[36,31:35]}
%lineColor=[0.2,.2,.2;   brewermap(9,'blues')];
%lineColor=lineColor([1,5:9],:)

%intCmb={[32,20,25,27]};
intCmb={[5,6,7,9]}
lineColor=lines(4);

intCmb={[4,2,13,11]}
lineColor=lines(4);
%
c=0;
for ce=1:numel(intCmb)
    intExp=intCmb{ce}
    figure
    for g=1:numel(intGene)
        c=c+1;
        subplot(2,2,c)
        hold off
        selRegion=quantile(acol(tssIdx(intGene{g})+[-1000,1000]),[0 1]);
        for e=1:numel(intExp)
            plot(selRegion(1):selRegion(2),conv2(meanProfile(selRegion(1):selRegion(2),intExp(e)),[1:50 49:-1:1]'.^2./sum([1:50 49:-1:1].^2),'same'),'LineWidth',2,'DisplayName',strjoin(expTypes{intExp(e),[1:4]}),'Color',lineColor(e,:))
            hold on
        end
        for gg=1:numel(intGene{g})
            posg=GP.chrIdx(GP.gene_infoR64.position(intGene{g}(gg),1))+GP.gene_infoR64.position(intGene{g}(gg),2:3);
            plot(posg,1*[1 1],'k-','Linewidth',5)
            text(mean(posg),1.5,GP.gene_infoR64.name(intGene{g}(gg)),'HorizontalAlignment','center','VerticalAlignment','bottom')
        end
        %selNucs=nucMeta.pos>selRegion(1) &nucMeta.pos<selRegion(2);
        %scatter(nucMeta.pos(selNucs),.5*ones(sum(selNucs),1),'ro','filled')
        %set(gca,'YScale','log')
        axis('tight')
        title(strjoin(GP.gene_infoR64.name(intGene{g})'))
        ylabel('signal')
        %ylim([1 40])
        if g==1
            legend()
        end
        
    end    
end
save_gf(gcf,'SF4_mutHisAsyncHA')

save_gf(gcf,'Fig4_mtAalphaHis')
save_gf(gcf,'Fig4_mtH2aHis')

[~,~,smeta.sid]=unique(smeta(:,[8,9,10,11]),'stable')
for i=1:max(smeta.sid)        
    myci=find(strcmp(smeta.ab,'myc') & smeta.sid==i);
    hai=find(strcmp(smeta.ab,'ha') & smeta.sid==i);
    sampleTo(:,i)=log(data(:,myci)+1)-log(data(:,hai)+1);
end


%% promoter scatter
metaProfile=meta_profile(meanProfile,'promoter',600,'useORF',false,'afterTSS',600,'scaled',false,'disCDS',false);
meanProm=squeeze(mean(metaProfile.cube(:,1:600,:),2))
hasTss=~isnan(meanProm(:,1));
GP=load('group_imp.mat');
%myc lvl
xExp=[18,18,13]
yExp=[14,15,11]
figName='Fig4_promoterMut'
%ha lvls

xExp=[9,9,4]
yExp=[5,6,2]
figName='SF4_promterMutHA'

for i=1:numel(xExp)
    subplot(1,numel(xExp),i)
    hold off
    dscatter(meanProm(hasTss,xExp(i)),meanProm(hasTss,yExp(i)))
    hold on
    scatter(meanProm(GP.groups{23}{2}{45},xExp(i)),meanProm(GP.groups{23}{2}{45},yExp(i)),'o','filled')
    text(meanProm(GP.groups{23}{2}{45},xExp(i)),meanProm(GP.groups{23}{2}{45},yExp(i)),GP.gene_infoR64.name(GP.groups{23}{2}{45}))
    xlabel(strjoin(expTypes{xExp(i),:}))
    ylabel(strjoin(expTypes{yExp(i),:}))
    title(corr(meanProm(hasTss,xExp(i)),meanProm(hasTss,yExp(i))))
end
save_gf(gcf,figName)

intExp=[13,10,11,18,14,15]
figure
for e=1:numel(intExp)
    subplot(2,3,e)
    intSmp=smeta.expId==intExp(e) & smeta.bad==0;
    imagesc(squeeze(median(metaProfile.cube(:,:,intSmp),1,'omitnan'))','XData',[-600 600])
    title(strjoin(expTypes{intExp(e),:}))
    yticks(1:sum(intSmp))
    yticklabels(smeta.name(intSmp))
end
save_gf(gcf,'Fig4_metaProfile')


end