function figure1_new()
clear all
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/nucMeta015AddFried.mat')
[data,smeta]=getNaamaData();
h3B={'fast1b-wt-async-x18','fast2b-wt-async-x18','fast2a-wt-async-x22','fast2b-wt-async-x22'};
smeta.tag2(contains(smeta.name,h3B))={'hb3'};
nucExt=50;
nucPos=zeros(size(data,1),1);
for i=1:size(nucMeta,1)
    posi=nucMeta.pos(i)-nucExt:nucMeta.pos(i)+nucExt;
    if i>1
        posi=posi(posi>0 & posi>mean(nucMeta.pos([i-1,i])) & posi<=size(data,1));
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
clear data i nucExt nucPos posi
highOcc=find(any(movmean(nucData(:,strcmp(smeta.ab,'ha')& strcmp(smeta.gt,'wt')),5)>30,2));
nucData(unique(acol(highOcc+[-2:2])),:)=NaN;
lowOcc=find(any(movmean(nucData(:,strcmp(smeta.ab,'ha') & strcmp(smeta.gt,'wt')),3)<0.01,2));
nucData(unique(acol(lowOcc+[-1:1])),:)=NaN;
[smeta,idx]=sortrows(smeta,[7,14,15]);
nucData=nucData(:,idx);
clearvars  highOcc lowOcc idx h3B


imagesc(corr(nucData,'rows','pairwise'))
title('corr nucleosomes')
colorbar()
yticks(1:numel(smeta.tag2))
yticklabels(smeta.name)
smeta.tag2(contains(smeta.tag2,'h4')&contains(smeta.exp,'28'))=strcat(smeta.tag2(contains(smeta.tag2,'h4')&contains(smeta.exp,'28')),'new')
clear data posi
%% averaging nucDat
%% find rpts
[expTypes,~,smeta.expId]=unique(smeta(:,[7,9,10,14]));
for i=unique(smeta.expId)'   
    selRpt=find(smeta.expId==i & smeta.bad==0);
    selNucs=all(nucData(:,selRpt)>1,2);
    [~,bestRpt]=max(sum(corr(log(nucData(selNucs,selRpt))),2));
    bestRpt=selRpt(bestRpt)
    for j=selRpt'
        pj=median(diff(log(nucData(selNucs,[j bestRpt])),1,2))
        fac(j)=exp(pj)
    end
    medNuc(:,i)=median(nucData(:,selRpt).*fac(selRpt),2);
end
clearvars i j pj fac h3b bestRpt selNucs selRpt
save('forFigure1plus.mat')
clear all
load('forFigure1plus.mat')
load('forFigure15.mat')
%%comapre Ha / myc samples
haExp=[11,10,8,9,12,4,3,7,5,1,2];
mycExp=haExp+12;
sampleOrder=[56,55,53,54,58,37,34,48,47,15,16,5,110,109,107,108,112,91,88,102,101,69,70,59]
% c=0;
% for tag={'haExp','mycExp'}
%     c=c+1;
%     subplot(1,2,c)
%     tagSmp=eval(tag{1});
%     imagesc(corr(medNuc(:,tagSmp),'rows','pairwise'),[0 1])
%     yticks(1:numel(tagSmp))
%     ylabel(tag{1})
%     yticklabels(expTypes.tag2(tagSmp))
%     colorbar()
% end
% close all
figure;
imagesc(corr(medNuc(:,sampleOrder),'rows','pairwise'))
yticks(1:numel(sampleOrder))
yticklabels(strcat(expTypes.ab(sampleOrder),'-',strrep(expTypes.tag2(sampleOrder),'_',''),'-',num2str(expTypes.n(sampleOrder))))
figName='SF1_HAmycMedCorr'
save_gf(gcf,figName)

mycExp=
clear bestRpt selRpt pj i j ans selnucs
%xExp=[7 3] ;yExp=[15,11] for Figure 1
xExp=5;yExp=16 % for Figure 1S
for i=1:numel(xExp)
    subplot(1,numel(xExp),i)
    selNucs=all(~isnan(medNuc(:,[xExp(i),yExp(i)])),2);
    dscatter(medNuc(selNucs,xExp(i)),medNuc(selNucs,yExp(i)))
    xlabel(strjoin(expTypes{xExp(i),[1,3,4]}))
    ylabel(strjoin(expTypes{yExp(i),[1,3,4]}))
    title(corr(medNuc(selNucs,xExp(i)),medNuc(selNucs,yExp(i))))
end
save_gf(gcf,'SF1_slowH3')
%% import external datasets
%% import external data
GP=load('group_imp.mat');
dion=readtable('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/extData/dionTO.bed','FileType','text','ReadVariableNames',false);
[chrNames,~,dion.chrId]=unique(dion.Var1,'stable');
dion.pos=GP.chr_idx(dion.chrId)'+mean([dion.Var2 dion.Var3],2);

rufiange=readtable('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/extData/rufiange2007TOScer3.bed','FileType','text','ReadVariableNames',false);
[~,rufiange.chrID]=ismember(rufiange.Var1,chrNames);
rufiange.chrID(rufiange.chrID==0)=17;
rufiange.pos=mean([rufiange.Var2 rufiange.Var3],2)+GP.chrIdx(rufiange.chrID);

[dist,met2dion]=min(abs(nucMeta.pos-dion.pos'),[],2);
[dist2,dion2met]=min(abs(nucMeta.pos'-dion.pos),[],2);

met2dion(dion2met(met2dion)~=[1:numel(met2dion)]'|dist>80)=nan;
nucDion=nan(size(nucMeta,1),1);
nucDion(~isnan(met2dion))=dion.Var4(met2dion(~isnan(met2dion)));
clear met2dion dist dist2 dion2met dion

nucExt=50;
nucBorders=[nucMeta.pos(1)-nucExt;movmean(nucMeta.pos,2,'Endpoints','discard');nucMeta.pos(end)+nucExt];
for i=1:ceil(size(rufiange,1)/1000);
    seli=i*1000-999:i*1000;
    seli=seli(seli<=size(rufiange,1));
    rufiange.nuc(seli)=sum(rufiange.pos(seli)>nucBorders',2);
end
rufiange.dist(rufiange.nuc<size(nucMeta,1))=(rufiange.pos(rufiange.nuc<size(nucMeta,1))-nucMeta.pos(rufiange.nuc(rufiange.nuc<size(nucMeta,1))));
rufiange.dist(rufiange.nuc==0 | rufiange.nuc>size(nucMeta,1))=NaN;
clear seli nucBorders pj nucPos chrNames
nucRufiange=accumarray(rufiange.nuc(abs(rufiange.dist)<nucExt),rufiange.Var4(abs(rufiange.dist)<nucExt),[size(nucMeta,1) 1],@mean,NaN);

[~,~,expTypes.tagID]=unique(expTypes.tag2,'stable')
nExp=size(expTypes,1)
for i=1:max(expTypes.tagID)
    myci=find(strcmp(expTypes.gt,'wt') &strcmp(expTypes.ab,'myc') & expTypes.tagID==i);
    hai=find(strcmp(expTypes.gt,'wt') &strcmp(expTypes.ab,'ha') & expTypes.tagID==i);
    expTypes(i+nExp,:)=expTypes(hai,:);
    expTypes.ab{i+nExp}='logto';
    medNuc(:,i+nExp)=log(medNuc(:,myci)+1)-log(medNuc(:,hai)+1);
end
subplot(2,3,4)
toCmp=19;
selNucs=~isnan(nucDion) & ~isnan(medNuc(:,toCmp));
x=dscatter(nucDion(selNucs),medNuc(selNucs,toCmp))
hold on
scatter(nucDion(ismember(nucMeta.gene,GP.groups{23}{2}{45}) & nucMeta.order<0),medNuc(ismember(nucMeta.gene,GP.groups{23}{2}{45}) & nucMeta.order<0,toCmp),'ro','filled','DisplayName','histone promoter nucs')
xlabel('Dion TO')
ylabel(strjoin(expTypes{toCmp,[4,1]}));
axis tight

selNucs=~isnan(nucRufiange) & ~isnan(medNuc(:,toCmp));
x=dscatter(nucRufiange(selNucs),medNuc(selNucs,toCmp))
hold on
scatter(nucRufiange(ismember(nucMeta.gene,GP.groups{23}{2}{45}) & nucMeta.order<0),medNuc(ismember(nucMeta.gene,GP.groups{23}{2}{45}) & nucMeta.order<0,toCmp),'ro','filled','DisplayName','histone promoter nucs')
xlabel('Rufiange TO')
ylabel(strjoin(expTypes{toCmp,[4,1]}));
axis tight
title(corr(nucRufiange(:,1),medNuc(:,toCmp),'rows','pairwise'))
save_gf(gcf,'SF1_CmpToRufiange')
%% compare HA between samples
figure
xExp=[3,3,3];
yExp=[1,2,4];
nucType(nucMeta.order==0,1)=1;
nucType(nucMeta.order>0,1)=2;
nucType(nucMeta.order<0,1)=3;
typeColor=[.3 .3 .3;lines(2)];
typeAlpha=[.3 1 1];
typeName={'N/A','ORF','promoter'}
for i=1:numel(xExp)
    subplot(1,numel(xExp),i)
    hold off
    ok = all(~isnan(medNuc(:,[xExp(i) yExp(i)])),2);
    dscatter(medNuc(ok,xExp(i)),medNuc(ok,yExp(i)))
    %ol = ok & any(medNuc(:,[xExp(i) yExp(i)])>20,2) & ((medNuc(:,xExp(i))>2*medNuc(:,yExp(i)))|(medNuc(:,yExp(i))>2*medNuc(:,xExp(i))));
    %hold on
%     for j=2:max(nucType)
%         subplot(2,3,(i-1)*3+j)
%         %scatter(medNuc(nucType==j&o,xExp(i)),medNuc(nucType==j&ol,yExp(i)),[],typeColor(j,:),'o','filled','DisplayName',typeName{j})        
%         dscatter(medNuc(nucType==j&ok,xExp(i)),medNuc(nucType==j&ok,yExp(i)),'plottype','scatter')
%         hold on
%     end
    xlabel(strjoin(expTypes{xExp(i),[1,3,4]}))
    ylabel(strjoin(expTypes{yExp(i),[1,3,4]}))
    title(corr(medNuc(ok,xExp(i)),medNuc(ok,yExp(i))))
    hold on
    plot(xlim,xlim,'k--')
end
save_gf(gcf,'SF1_HAcmpAll')
%% turnover correlations between individual sampels
[~,~,smeta.sid]=unique(smeta(:,[8,9,10,11]),'stable')
sampleTo=nan(size(nucMeta,1),max(smeta.sid));
for i=1:max(smeta.sid)        
    myci=find(strcmp(smeta.ab,'myc') & smeta.sid==i);
    hai=find(strcmp(smeta.ab,'ha') & smeta.sid==i);
    sampleTo(:,i)=log(nucData(:,myci)+1)-log(nucData(:,hai)+1);
end
%% allc orrelations for supplement
subplot(1,2,1)
imagesc(corr(medNuc,'rows','pairwise'))
yticks(1:24)
yticklabels(strrep(strcat(expTypes.ab,'-',expTypes.tag2),'_',' '))
[~,idx]=sortrows(smeta(:,[7,14]))
subplot(1,2,2)
imagesc(corr([nucData(:,idx),sampleTo(:,idx(1:44))],'rows','pairwise'))
save_gf(gcf,'SF1_allCrs')
%%
[tags,~,smeta.tagid]=unique(smeta.tag2,'stable');
scale=max(accumarray(smeta.tagid(smeta.bad==0 &strcmp(smeta.ab,'ha')),1));
intTags=[2,4,1,6];
for i=1:max(smeta.tagid)
    selTO=smeta.sid(smeta.tagid==i & smeta.bad==0 & strcmp(smeta.ab,'ha'));
    selTO=repmat(selTO,ceil(scale./numel(selTO)),1);
    collTO{i}=sort(selTO(1:scale));    
end
crMat=corr(sampleTo(:,cat(1,collTO{intExp})),'rows','pairwise');
subplot(1,3,3)
hold off
imagesc(crMat)
yticks(scale/2:scale:size(crMat,1))
yticklabels(tags(intTags))
xticks(scale/2:scale:size(crMat,1))
xticklabels(tags(intTags))
ylabel(colorbar('west'),'all nuc TO correlation')

%% turnover correlations median levels
intExp=[21,22,23,24,19,20]
crMat=corr(medNuc(:,intExp),'rows','complete');
subplot(1,3,3)
hold off
imagesc(crMat)
yticks(1:numel(intExp))
yticklabels(expTypes.tag2(intExp))
xticks(1:numel(intExp))
xticklabels(expTypes.tag2(intExp))
ylabel(colorbar('west'),'all nuc TO correlation')
save_gf(gcf,'/Fig1_medTOcorrelations')

%% compare H2A and H2B
subplot(1,2,1)
hold off
intExp={[1,2],[9,10],[3,4],[11,12]}
intNuc=[-1,1,2]
for e=1:numel(intExp)
    subplot(2,2,e)
    for n=1:numel(intNuc)
        selNucs=nucMeta.order==intNuc(n)& all(~isnan(medNuc(:,intExp{e})),2);
        p=linortfit2(medNuc(selNucs,intExp{e}(1)),medNuc(selNucs,intExp{e}(2)))
        scatterName=sprintf('%.2fx+%.2f',p(1),p(2))
        scatter(medNuc(selNucs,intExp{e}(1)),medNuc(selNucs,intExp{e}(2)),'.','DisplayName',scatterName);
        hold on        
    end
    legend()
end

hold on
scatter(medNuc(nucMeta.order==-1,1),medNuc(nucMeta.order==-1,2),'o','filled')
subplot(1,2,2)
hold off
scatter(medNuc(nucMeta.order==1,9),medNuc(nucMeta.order==1,10),'.')
hold on
scatter(medNuc(nucMeta.order==-1,9),medNuc(nucMeta.order==-1,10),'o','filled')


%% profiles

clear all
[data,smeta]=getNaamaData();
[~,~,smeta.sid]=unique(smeta(:,[8,9,10,11]),'stable')
for i=1:max(smeta.sid)        
    myci=find(strcmp(smeta.ab,'myc') & smeta.sid==i);
    hai=find(strcmp(smeta.ab,'ha') & smeta.sid==i);
    sampleTo(:,i)=log(data(:,myci)+1)-log(data(:,hai)+1);
end
GP=load('group_imp.mat');
load('cyclebase.mat','cyclebase')
metaProfile=meta_profile(sampleTo,'promoter',600,'useORF',false,'afterTSS',700,'scaled',false,'disCDS',true);
geneGroups={1:6701,GP.groups{6}{2}{1},find(cyclebase.data(:,1)<200),GP.groups{7}{2}{66}}
groupNames={'All gene*','ESR ind','CC Top 200','RiBi'}
intTag={'h3','h4','h2a','h2b'}
cLim=[0,.9;0 .5;0 .5;0 .5];
close all
figure;
for t=1:numel(intTag)
    selSmp=smeta.sid(smeta.bad==0 & strcmp(smeta.tag2,intTag{t}) & strcmp(smeta.ab,'ha'));
    for i=1:numel(geneGroups)
        subplot(numel(intTag),numel(geneGroups),(t-1)*numel(geneGroups)+i)
        imMat=permute(median(metaProfile.cube(geneGroups{i},:,selSmp),1,'omitnan'),[3,2,1]);
        imagesc(imMat,'XData',[-600 700])
        yticks([])
        xticks([-500,0,500])
        if t==1
            title(groupNames{i})
        end
        ylabel([intTag{t} '-logTO'] )
        colorbar()
        caxis(cLim(t,:))
    end
end
save_gf(gcf,'Fig1_groupProfileAll')

%% create mean Profiles
[expTypes,~,smeta.expId]=unique(smeta(:,[14,7,9,10]),'stable')
for i=unique(smeta.expId(smeta.bad==0))'
    selSmp=smeta.expId==i &smeta.bad==0;
    meanProfile(:,i)=mean(data(:,selSmp),2);
end
clear data
for i=find(strcmp(expTypes.ab,'ha'))'
    myci=find(strcmp(expTypes.ab,'myc')& strcmp(expTypes.tag2,expTypes.tag2(i)) & strcmp(expTypes.con,expTypes.con(i)) & strcmp(expTypes.gt,expTypes.gt(i)))
    meanTo(:,i)=diff(log(meanProfile(:,[i myci])+1),1,2);
end

metaProfile=meta_profile(meanTo,'promoter',600,'useORF',false,'afterTSS',700,'scaled',false,'disCDS',false);
metaProfile=meta_profile(meanProfile,'promoter',600,'useORF',false,'afterTSS',700,'scaled',false,'disCDS',false);

hasProfile=all(~isnan(metaProfile.cube(:,:,1)),2)
regev=load ('/home/labs/barkailab/LAB/data/DataExternal/EXPRESSION/regev_pnas_2008_rnaseq/data');
geneLvl=median(regev.regev_rnaseq,2,'omitnan');
load('holstege.mat')
geneStdH=std(holstege.TsXMut,[],2,'omitnan');
clear holstege regev
lvlIdx=sum(geneLvl>=quantile(geneLvl(hasProfile),[0:1/6:5/6]),2);
stdIdx=sum(geneStdH>=quantile(geneStdH(hasProfile),[0:1/6:5/6]),2);
intIdx={'lvlIdx','stdIdx'};
figure;
for t=1:numel(intTag)
    %selSmp=smeta.sid(smeta.bad==0 & strcmp(smeta.tag2,intTag{t}) & strcmp(smeta.ab,'ha'));
    %selExp=find(strcmp(expTypes.tag2,intTag{t}) & strcmp(expTypes.ab,'ha'));
    selExp=find(strcmp(expTypes.tag2,intTag{t}) & strcmp(expTypes.ab,'myc'));
    for i=1:numel(intIdx)
        subplot(numel(intIdx),numel(intTag),(i-1)*numel(intTag)+t)
        selIdx=eval(intIdx{i});
        imMat=nan(max(selIdx),size(metaProfile.cube,2));
        for j=1:max(selIdx)
            %imMat(j,:)=mean(mean(metaProfile.cube(selIdx==j,:,selSmp),'omitnan'),3,'omitnan');
            imMat(j,:)=mean(metaProfile.cube(selIdx==j,:,selExp),'omitnan');
        end
        imagesc(imMat,'XData',[-600 700])
        ylabel(intIdx{i})
        if i==1
            title(intTag{t})
        end
        colorbar()
        %caxis(cLim(t,:))
    end
end
save_gf(gcf,'Fig1_MeanLvlStdProfileAll')
close all

%% compare Epimarks
clear all
load('forFigure4add.mat')
load('epMarks.mat')
nucMark=nucMark/100;
nucMark(nucMark<.5)=NaN;
fig12=load('forFigure12.mat')
asyncGilad=find(contains(expTypes.ab,{'k56','k4','k79','k9','myc'}) & contains(expTypes.tag2,'h3')& contains(expTypes.con,'async') & contains(expTypes.gt,'wt'));
asyncGilad=asyncGilad([2,4,1,3])
intFried=[14,17,13,16]%13,14,16,17];
figure
imagesc(corr(medNuc(:,asyncGilad),nucMark(:,intFried),'rows','pairwise'))
xticks(1:numel(intFried))
xticklabels(markNames.Var1(intFried))
xticklabels(markNames.Var1(intFried))
yticks(1:numel(asyncGilad))
yticklabels(expTypes.ab(asyncGilad))
ylabel('Wt async Marks')
xlabel('Friedmann marks')
figName='SF1_FriedmanMarks';
colorbar()
save_gf(gcf,figName)
clear all
%%markrofiles
clear all
[data,smeta]=getNaamaData('figure',4);
tooLow=sum(movmean(data(:,contains(smeta.gt,'wt')&contains(smeta.ab,'ha')),100)<0.05,2)>sum(contains(smeta.gt,'wt')&contains(smeta.ab,'ha'))/3;
data(tooLow,:)=NaN;
[expTypes,~,smeta.expId]=unique(smeta(:,[14,7,9,10]),'stable')
for i=unique(smeta.expId(smeta.bad==0))'
    selSmp=smeta.expId==i &smeta.bad==0;
    meanProfile(:,i)=mean(data(:,selSmp),2,'omitnan');
end
clear data smeta tooLow selSmp intData i 
metaProfile=meta_profile(meanProfile,'promoter',600,'useORF',false,'afterTSS',700,'scaled',false,'disCDS',false);
figName='SF1_MarksMeta'
intAb={'k56','k9ac','k4m3','k79m3'}
for i=1:numel(intAb)
    subplot(numel(intAb),1,i)
    hold off
    for j= find(contains(expTypes.ab,intAb{i}) & contains(expTypes.con,{'async'}) & contains(expTypes.gt,{'wt'}))'
        plot(-600:700,mean(metaProfile.cube(:,:,j),'omitnan'),'DisplayName',strjoin(expTypes{j,:}),'Linewidth',2)
        hold on
    end
    title(intAb{i})
    axis tight
    legend('show','Location','best')
    ylabel('Intensity')
end
save_gf(gcf,figName)
end