%% review()
%[expTypes,medNuc,smeta,nucData,nucMeta]=standardImport(16)
clear all
load('forFigure15.mat')
GP=load('group_imp.mat');
ToR = readtable('/home/labs/barkailab/yuliago/Documents/MATLAB/ToR_raz.csv');
ToR.idx=GP.chrIdx(ToR.chr)+ToR.start+500;
[val,idx]=min(abs(nucMeta.pos-ToR.idx'),[],2);
nucMeta.ToR(val<501)=ToR.score(idx(val<501));
nucMeta.ToR(nucMeta.ToR==0)=NaN;
clear ToR
figure;
c=0;
clv=30;
for i=find(contains(expTypes.con,{'async'})&contains(expTypes.ab,'ha','IgnoreCase',true)&contains(expTypes.gt,{'wt'},'IgnoreCase',true))'
    c=c+1;
    subplot(4,4,c)
    hold off
    myci=find(strcmp(expTypes.con,expTypes.con{i})&contains(expTypes.ab,'myc','IgnoreCase',true)&strcmp(expTypes.tag2,expTypes.tag2{i})&strcmp(expTypes.gt,expTypes.gt{i}))
    scatter(medNuc(:,i),medNuc(:,myci),[],nucMeta.ToR,'.')
    selNucs=all(~isnan(medNuc(:,[i,myci])),2) & nucMeta.ToR>30;
    %dscatter(medNuc(selNucs,i),medNuc(selNucs,myci))
    title(strjoin(expTypes{i,[3,4]}))
    p=robustfit(medNuc(selNucs,i),medNuc(selNucs,myci))
    hold on    
    plot([0:30],[0:30]*p(2)+p(1),'r-','Linewidth',2)
    title(sprintf('%s (%d/%.0f)',strjoin(expTypes{i,[2,3,4]}),expTypes.cc(i),clv))
    %title(sprintf('%s,%.2f%+.2f',strjoin(expTypes{i,[3,4]}),p(2),p(1)))
    text(20,50,sprintf('%.2fx%+.2f',p(2),p(1)),'Color',[1 0 0])
    [~,topNucs]=maxk(diff(log(medNuc(:,[i myci])+2),1,2).*ismember(nucMeta.gene,GP.groups{23}{2}{45}),7)
    
    %p=robustfit(medNuc(topNucs,i),medNuc(topNucs,myci))
    hold on    
    if expTypes.cc(i)>0
        plot([0:30],[0:30]*p(2)*expTypes.cc(i)/clv+p(1),'k--','Linewidth',1) %theoretical 100% line   
        frac=(medNuc(topNucs,myci)-p(1))./(p(2)*medNuc(topNucs,i))*clv/expTypes.cc(i);
        text(medNuc(topNucs,i),medNuc(topNucs,myci),num2str(frac*100,'%.0f'))
    end    
    scatter(medNuc(topNucs,i),medNuc(topNucs,myci),'ro','filled')
    %text(10,120,sprintf('%.2f%+.2f',p(2),p(1)),'Color',[0 0 0])
    xlim([0 30])
    ylim([0 40])
end

xCmp=4
figure
c=0
for i=setdiff(find(strcmp(expTypes.gt,'wt')&strcmp(expTypes.con,'async')&contains(expTypes.ab,{'ha','h3'})),xCmp)'
    c=c+1;
    subplot(4,5,c)
    scatter(medNuc(:,xCmp),medNuc(:,i),[],nucMeta.ToR,'.')
    %selNucs=all(~isnan(medNuc(:,[i,myci])),2);
    hold on
    scatter(medNuc(topNucs,xCmp),medNuc(topNucs,i),'ro','filled')    
    xlabel(strjoin(expTypes{xCmp,[1,2,3,4]}))
    ylabel(strjoin(expTypes{i,[1,2,3,4]}))
    xlim([0 30])
    ylim([0 30])
    plot([0 30],[0 30],'k--')
end
suptitle(strjoin(expTypes{xCmp,[1,2,3,4]}))
close all
figure;
xId=85;
yId=108;
scatter(medNuc(:,xId),medNuc(:,yId),[],diff(log(medNuc(:,[14 66])+1),1,2),'o','filled')
scatter(medNuc(:,xId),medNuc(:,yId),[],medNuc(:,[14]),'o','filled')

caxis([-1 1])
selnucs=medNuc(:,xId)>0&medNuc(:,yId)>0
p=robustfit(medNuc(selNucs,yId),medNuc(selNucs,xId))
hold on
plot([0 60]*p(2)+p(1),[0 60],'LineWidth',2)
scatter(medNuc(topNucs,xId),medNuc(topNucs,yId),'ro','Linewidth',2)
xlabel(strjoin(expTypes{xId,[1,2,3,4]}))
ylabel(strjoin(expTypes{yId,[1,2,3,4]}))
ylabel(colorbar(),'log(H2Ato)')
p=robustfit(medNuc(topNucs,xId),medNuc(topNucs,yId),[],[],0)%robustfit(medNuc(topNucs,xId),medNuc(topNucs,yId),'cauchy)

selNucs=(medNuc(:,85)>5)&(medNuc(:,33)>5)
plot([0 30],[0 30]*p(1),'LineWidth',2)

scatter(medNuc(selNucs,14),medNuc(selNucs,66),[],diff(log(medNuc(selNucs,[33 85])+1),1,2),'o','filled')
median(diff(log(medNuc(selNucs,[33 85])+1),1,2))
hold on
figure
%xId=49;yId=43;cId=81;%
%xId=102;yId=97;cId=81;
%xId=101;yId=95;cId=81;
xId=46+52;yId=41+52;cId=75;
xId=45+52;yId=39+52;cId=75;
xId=41;yId=38;cId=75;
xId=29;yId=17;cId=75;
%xId=25;yId=14;cId=75;


subplot(1,2,2)
hold off
scatter((medNuc(:,xId)+1),(medNuc(:,yId)+1),[],medNuc(:,cId),'o','filled')
ylabel(colorbar(),strjoin(expTypes{cId,1:4}))
hold on
%scatter(medNuc(topNucs,xId),medNuc(topNucs,yId),'ro','filled')
hold on
%text(medNuc(topNucs,xId),medNuc(topNucs,yId),GP.gene_infoR64.name(nucMeta.gene(topNucs)))
ylabel(strjoin(expTypes{yId,1:4}))
xlabel(strjoin(expTypes{xId,1:4}))
caxis([0 20])
set(gca,'Color',[.7 .7 .7])
selNucs=all(medNuc(:,[xId,yId])>0,2)
p=linortfit2(medNuc(selNucs,xId),medNuc(selNucs,yId))
hold on
%plot(xlim,xlim*p(1)+p(2),'k--')
%text(mean(xlim),mean(ylim),sprintf('%.2fx%+.2f',p(1),p(2)))


%% vplots
clear all

% figure
% imagesc(corr(medNuc,'rows','pairwise'))
% set(gca,'YTick',1:size(expTypes,1),'YTickLabel',strrep(strcat(expTypes.ab,'-',expTypes.tag2),'_',' '),'XTick',1:size(expTypes,1),'XTickLabel',strrep(strcat(expTypes.ab,'-',expTypes.tag2),'_',' '),'XTickLabelRotation',90)
% title('Genome-wide correlaiton between all WT async systems (myc and HA)')
% ylabel(colorbar(),'Pearson correlaiton')
load('forFigure15.mat')
load('forFigure2add.mat')
load('forFigure3add.mat','nucData','smeta','nucMeta')
medNuc=nucData;
expTypes=smeta(:,[7,9,10,14])
expTypes.tp=str2double(regexp(expTypes.con,'(?<=T)\d+','match','once'))
expTypes.tp(isnan(expTypes.tp))=-1
[expTypes,idx]=sortrows(expTypes,[1,4,5]);
medNuc=medNuc(:,idx);
clear idx
regev=load ('/home/labs/barkailab/LAB/data/DataExternal/EXPRESSION/regev_pnas_2008_rnaseq/data');
geneLvl=median(regev.regev_rnaseq,2,'omitnan');
load('holstege.mat')a
geneStdH=std(holstege.TsXMut,[],2,'omitnan');
meanChg=median(holstege.TsXMut,2,'omitnan');
GP=load('group_imp.mat');
ToR = readtable('/home/labs/barkailab/yuliago/Documents/MATLAB/ToR_raz.csv');
ToR.idx=GP.chrIdx(ToR.chr)+ToR.start+500;
[val,idx]=min(abs(nucMeta.pos-ToR.idx'),[],2);
nucMeta.ToR(val<501)=ToR.score(idx(val<501));
nucMeta.ToR(nucMeta.ToR==0)=NaN;
load('spell.mat')
%geneStdS=std(spell.data(:,spell.abs==1),[],2);
clear holstege regev ToR idx val spell
load('holstege.mat')
dirStd=geneStdH.*sign(meanChg);
dirStd(abs(meanChg)<2*geneStdH./sqrt(sum(~isnan(holstege.TsXMut),2)))=NaN;
dirStdClean=dirStd;
dirStdClean(abs(dirStd)>quantile(abs(dirStd),.95))=NaN
load('addData.mat','tsScore')
asLvl=log2(tsScore(:,2)+.1);
clearvars -except asLvl dirStdClean geneLvl dirStd medNuc expTypes smeta nucData nucMeta
GP=load('group_imp.mat');

figureType(1).intNuc={-2,-1,1,2};
figureType(1).intPar={'dirStdClean'}
figureType(1).movPar=[0.05]
figureType(1).seperateCr=[1];
figureType(1).xLim=[-.25 .25]
figureType(1).normLvl=false;
figureType(1).normRange=[-0.1 0.1]
figureType(1).badNucs=false(size(nucMeta.ToR));

figureType(2).intPar={'geneLvl'}
figureType(2).movPar=[0.75]
figureType(2).seperateCr=0
figureType(2).intNuc={-2,-1,1,2}
figureType(2).xLim=[5 13]
figureType(2).normLvl=false;
figureType(2).normRange=[4.5 5.5]
figureType(2).badNucs=false(size(nucMeta.ToR));%ismember(nucMeta.gene,find(geneLvl>10 | abs(dirStd)>0.15));

figureType(3).intNuc={-2,-1,1,2};
figureType(3).intNuc={min(nucMeta.order):max(nucMeta.order)}
%figureType(3).intNuc={-2,-1,1,2}
figureType(3).intPar={'nucMeta.ToR'}
figureType(3).movPar=1
figureType(3).seperateCr=0
figureType(3).xLim=[17 40]
figureType(3).normLvl=false;
figureType(3).normRange=[-0.1 0.1]
figureType(3).badNucs=min(abs(nucMeta.pos-GP.oris.loc'),[],2)<2000 ;
%figureType(3).badNucs=~ismember(1:size(nucMeta,1),top500)';false(size(nucMeta.ToR));

figureType(4).intPar={'asLvl'}
figureType(4).movPar=[0.5]
figureType(4).seperateCr=0
figureType(4).intNuc={-3,-2,-1,1,2,3}
figureType(4).xLim=[-3.3 8]
figureType(4).normLvl=false;
figureType(4).normRange=[4.5 5.5]
figureType(4).badNucs=false(size(nucMeta.ToR));

figureType(5).intPar={'geneLvl'}
figureType(5).movPar=[0.75]
figureType(5).seperateCr=2
figureType(5).intNuc={-2,-1,1,2,3};%{-1 1 2};
figureType(5).xLim=[5 13]
figureType(5).normLvl=false;
figureType(5).normRange=[4.5 5.5]
figureType(5).badNucs=false(size(nucMeta.ToR));%ismember(nucMeta.gene,find(geneLvl>10 | abs(dirStd)>0.15));
figureType(5).yLim=[]

lineColor=[lines(7);0.3.*[1 1 1];0.6*[1 1 1]]
lineColor=[.3 .3 .3;brewermap(10,'Set1')]
globalBlackList=ismember(nucMeta.gene,GP.groups{23}{2}{45})% & ismember(nucMeta.gene,find(geneLvl>12));
hisNucs=[2733,2734,12699,12700,54681,54682,54683];
step=5;
selPlot=1:12
%rlf2 vs. wt
%[0.4980 0.2471 0.5961;0.3922 0.3922 0.3922;0 0.6824 0.9373;0.3922 0.3922 0.3922]
intCmb={[1 12],[2 13]}
intCmb={[15 15+54],[16 16+54 5 5+54],[34 34+54 48 48+54],[47 47+54]}% {[25,29,11],[25,29,11]+52}
intCmb={[15 10 11 12 7],[15 10 11 12 7]+93,[74 59 60 61 32],[75 59 60 61 33]+95};figName='geneExpAll_4RowsDark'
intCmb={[74 78 38],[75 78 37]+95}
%intCmb={[37+54 36+54],[37 36]}
%intCmb{1}=find(contains(expTypes.ab,'ha','IgnoreCase',true)& contains(expTypes.gt,'rpb','IgnoreCase',true) & expTypes.n>0)
%intCmb{2}=find(ismember(expTypes(:,2:4),expTypes(intCmb{1},2:4)) & contains(expTypes.ab,'myc'));
intCmb={[1,38,30,51],[2,39,34,52]}
%intCmb={[39,34,52],[39,34,52]+54}
intCmb={[2,3,39,34,52,50,47]}
intCmb={find(contains(expTypes.ab,'ha')&contains(expTypes.con,'async')&contains(expTypes.tag2,'h3x2')) find(contains(expTypes.ab,'myc')&contains(expTypes.con,'async')&contains(expTypes.tag2,'h3x2'))}
% lines
close all 
for cmb=1:numel(intCmb)
    intExp=intCmb{cmb}(:)';
    for f=1:3
        intNuc=figureType(f).intNuc;
        intPar=figureType(f).intPar;
        movPar=figureType(f).movPar;
        seperateCr=figureType(f).seperateCr;
        xLim=figureType(f).xLim;
        normLvl=figureType(f).normLvl;
        normRange=figureType(f).normRange;
        blackList=globalBlackList | figureType(f).badNucs;
        if cmb==1
            figLink(f)=figure('Color',[1 1 1]);
        else
            figure(figLink(f))
        end
        c=0;
        for p=1:numel(intPar)
            parData=eval(intPar{p});
            for n=1:numel(intNuc)
                if contains(expTypes.tag2(intExp(1)),'h2')
                    lineColor=[ 0.3922 0.3922 0.3922;
                        0.4980 0.2471 0.5961;
                        brighten([0.3922 0.3922 0.3922],0.6);
                      brighten([0.4980 0.2471 0.5961],0.6);];
                     %yLim=[4 18]
                     lineColor=brewermap(8,'Set2');                     
                else
                    lineColor=[
                        0.3922 0.3922 0.3922;
                        0 0.6824 0.9373
                        brighten([0.3922 0.3922 0.3922],0.6);
                        brighten([0 0.6824 0.9373],-.6)] ;
                    %lineColor=[0 0 0;brewermap(4,'Set2')]
                    %yLim=[4 16] 
                    lineColor=brewermap(16,'BuPu');
                    %lineColor=lineColor(2:end,:);
                    lineColor=brewermap(8,'Set2')      ;              

                end
                %subplot(numel(intPar),numel(intNuc),(p-1)*numel(intNuc)+n)
                c=c+1;
                subplot(numel(intCmb),numel(intNuc),c+(cmb-1)*numel(intNuc))
                hold off
                hisNucs=find(ismember(nucMeta.order,intNuc{n}) & ismember(nucMeta.gene,GP.groups{23}{2}{45}));
                if size(parData,1)==6701
                    selNucs=find(ismember(nucMeta.order,intNuc{n}) & ismember(nucMeta.gene,find(~isnan(parData))) &~blackList & all(~isnan(medNuc(:,intExp)),2) );
                    [xVal,idx]=sort(parData(nucMeta.gene(selNucs)));
                else
                    selNucs=find(ismember(nucMeta.order,intNuc{n}) & ~isnan(parData) &~blackList & all(~isnan(medNuc(:,intExp)),2) );
                    [xVal,idx]=sort(parData(selNucs));
                end
                selNucs=selNucs(idx);
                clear lineObj idx
                for e=1:numel(intExp)
                    expVal=medNuc(selNucs,intExp(e));%./exp(mean(log(medNuc(hisNucs,intExp(e)))));
                    if normLvl
                        expVal=expVal./mean(expVal(xVal>normRange(1) & xVal<normRange(2)));
                    end
                    %expVal=expVal./mean(expVal(abs(xVal-quantile(parData,0.05))<=movPar(p)/2));
                    yVal=movmean(expVal,movPar(p),'omitnan','SamplePoints',makeUnique(xVal,4));
                    nVal=movsum(~isnan(expVal),movPar(p),'omitnan','SamplePoints',makeUnique(xVal,4));
                    yValE=movstd(expVal,movPar(p),'omitnan','SamplePoints',makeUnique(xVal,4))./sqrt(movsum(~isnan(expVal),movPar(p),'omitnan','SamplePoints',makeUnique(xVal,4)));
                    if seperateCr(p)==1
                        cr1=corr(xVal(xVal<=0),expVal(xVal<=0),'rows','pairwise');
                        cr2=corr(xVal(xVal>=0),expVal(xVal>=0),'rows','pairwise');
                        lineName=sprintf('%s,%.2f,%.2f',strjoin(expTypes{intExp(e),[1,4,2,3]}),cr1,cr2);
                        idxN=find((xVal<max(xVal(xVal<0))-movPar(p)/2) & mod([1:numel(xVal)]',step)==1,1,'last');
                        idxP=find(xVal>min(xVal(xVal>0))+movPar(p)/2,1,'first')+step;
                        lineObj(e)=plotWithStd(xVal(1:step:idxN),yVal(1:step:idxN),yValE(1:step:idxN),'Linewidth',2,'Color',brighten(lineColor(e,:),0),'DisplayName',lineName);
                        hold on
                        [~]=plotWithStd(xVal(idxP:step:end),yVal(idxP:step:end),yValE(idxP:step:end),'Linewidth',2,'Color',brighten(lineColor(e,:),0),'DisplayName',lineName);
                        plot(xVal([idxN,idxP]),yVal([idxN,idxP]),'--','Linewidth',2,'Color',brighten(lineColor(e,:),0))
                    elseif seperateCr(p)==2
                        cr=corr(xVal,expVal,'rows','pairwise');
                        lineName=sprintf('%s,%.2f',strjoin(expTypes{intExp(e),[1,4,2,3]}),cr);
                        xNuc=1:numel(xVal);
                        lineObj(e)=plotWithStd(xNuc(1:step:end),yVal(1:step:end),yValE(1:step:end),'Linewidth',2,'Color',brighten(lineColor(e,:),0),'DisplayName',lineName);
                        hold on
                    else
                        cr=corr(xVal,expVal,'rows','pairwise');
                        lineName=sprintf('%s,%.2f',strjoin(expTypes{intExp(e),[1,4,2,3]}),cr);
                        if false
                            xEnd=find(movmean(nVal>median(nVal)/4,4)>0.5,1,'last');
                            xEnd2=find(movmean(nVal>median(nVal)/8,4)>0.5,1,'last');
                            plotWithStd(xVal(xEnd:step:xEnd2),yVal(xEnd:step:xEnd2),yValE(xEnd:step:xEnd2),'LineStyle','--','Linewidth',2,'Color',brighten(lineColor(e,:),0.5),'DisplayName',lineName);
                            hold on                             
                            plotWithStd(xVal(xEnd2:step:end),yVal(xEnd2:step:end),yValE(xEnd2:step:end),'LineStyle','--','Linewidth',2,'Color',brighten(lineColor(e,:),0.9),'DisplayName',lineName,'fillLine','-','fillAlpha',0);
                        else
                            xEnd=numel(xVal);
                        end
                        lineObj(e)=plotWithStd(xVal(1:step:xEnd),yVal(1:step:xEnd),yValE(1:step:xEnd),'Linewidth',2,'Color',brighten(lineColor(e,:),0),'DisplayName',lineName);                                       
                        hold on
                    end
                    %scatter(parData(nucMeta.gene(hisNucs)),medNuc(hisNucs,intExp(e)),[],lineColor(e,:),'o','filled')
                end
                axis tight
                if seperateCr(p)==2
                    xlabel(sprintf('Nucleosomes orderd by %s',intPar{p}))
                    %xticks(round([0.1:0.2:0.9].*numel(xVal)))
                    xticks(sum(xVal<[5.8,6.7,7.4,8.3,10.1],1))
                    xticklabels(round(xVal(xticks),1))
                    xlim([0.03 0.99]*numel(xVal))
                    %ylim(yLim)
                    %yticks(5:3:max(yLim))
                else
                    xlabel(intPar{p})
                    xlim(xLim)
                end                
                
                if numel(figureType(f).yLim)==2
                    ylim(figureType(f).yLim)
                end
                if numel(intNuc{n})==1
                    title(intNuc(n))
                end
                if 2==2
                    legend(lineObj,'Location','best')
                end
            end
            ylabel(sprintf('%.3f - movAcg',movPar(p)))
        end
    end
end
save_gf(gcf,figName)
%% correaltion low-imntermed
selGenes=setdiff(find((geneLvl<10) ),GP.groups{23}{2}{45})
intNucs=[-2,-1,1,2]
for n=1:numel(intNucs)
    selNucs=ismember(nucMeta.gene,selGenes) & nucMeta.order==intNucs(n);
    allCr(n,:)=corr(medNuc(selNucs,:),geneLvl(nucMeta.gene(selNucs)),'rows','pairwise')
end
allCr(:,[29 81])
%% qPCRprobes()
clear all
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/nucMeta015AddFried.mat')
GP=load('/home/labs/barkailab/felixj/Documents/MATLAB/scripts/gene/group_imp.mat');
load('forFigure15.mat')
qPCR=readtable('extData/qPCRStrubin.csv');
[probes,~,qPCR.id]=unique(qPCR.Var1);
for i=1:numel(probes)
    selProbes=find(qPCR.id==i & qPCR.Var3==100 & startsWith(qPCR.Var2,'B'));
    [best,idx]=max(qPCR.Var4(selProbes));
    if sum(qPCR.Var4(selProbes)==best)==1
        qPCR.best(selProbes(idx))=1;
    else
       selProbes(qPCR.Var4(selProbes)==best)
    end
end
qPCR=qPCR(qPCR.best==1,:);
chrNames={'BK006935.2','BK006936.2','BK006937.2', ...
'BK006938.2','BK006939.2','BK006940.2','BK006941.2', ...
'BK006934.2','BK006942.2','BK006943.2','BK006944.2', ...
'BK006945.2','BK006946.2','BK006947.3','BK006948.2', ...
'BK006949.2'}
[~,qPCR.chr]=ismember(qPCR.Var2,chrNames)
qPCR.mid=(qPCR.Var9+qPCR.Var10)/2+GP.chrIdx(qPCR.chr);
clear targets
[targets.name,~,qPCR.targetId]=unique(regexp(qPCR.Var1,'.*(?=-Fw)|.*(?=-Rev)','match','once'))
targets=struct2table(targets);
for i=find(accumarray(qPCR.targetId,1)==2)'
    targets.chr(i)=min(qPCR.chr(qPCR.targetId==i));
    targets.up(i)=GP.chrIdx(targets.chr(i))+min([qPCR.Var9(qPCR.targetId==i);qPCR.Var10(qPCR.targetId==i)]);
    targets.down(i)=GP.chrIdx(targets.chr(i))+max([qPCR.Var9(qPCR.targetId==i);qPCR.Var10(qPCR.targetId==i)]);
    targets.pos(i)=mean(qPCR.mid(qPCR.targetId==i))
    [targets.dis(i),targets.nuc(i)]=min(abs(nucMeta.pos-targets.pos(i)))
    targets.nucPos(i)=nucMeta.pos(targets.nuc(i));
    if nucMeta.gene(targets.nuc(i))>0
        targets.gene(i)=GP.gene_infoR64.name(nucMeta.gene(targets.nuc(i)))
    end
end
regev=load ('/home/labs/barkailab/LAB/data/DataExternal/EXPRESSION/regev_pnas_2008_rnaseq/data');
geneLvl=median(regev.regev_rnaseq,2,'omitnan');
load('qPCRtargets.mat')

orderedOligos={'HTA1','HTA2','HHF1','HHT2','HTB2','KOG1','HYP2','RPL3','NIT2','AUS1','ACS1'}
targets.gene(cellfun('isempty',targets.gene))={''}
targets.ordered=ismember(targets.gene,orderedOligos)
figure
c=0;
for xId=[11,31,49]
    c=c+1;
    subplot(1,3,c)
    yId=xId+52;
    xData=medNuc(:,xId);
    yData=medNuc(:,yId);
    hold off
    selNucs=nucMeta.gene>0;
    scatter(xData(selNucs),yData(selNucs),[],geneLvl(nucMeta.gene(selNucs)),'.')
    hold on   
    ylabel(colorbar(),strjoin(expTypes{3,1:4}))
    xlabel(strjoin(expTypes{xId,1:4}))
    ylabel(strjoin(expTypes{yId,1:4}))
    p=robustfit(xData,yData)
    hold on    
    plot([0:30],[0:30]*p(2)+p(1),'k-','Linewidth',2)
    
    %hlNucs=lowH2
    %scatter(xData(hlNucs),yData(hlNucs),'ro','filled')
    %text(xData(hlNucs),yData(hlNucs),GP.gene_infoR64.name(nucMeta.gene(hlNucs)))
    
    scatter(xData(targets.nuc(targets.ordered)),yData(targets.nuc(targets.ordered)),'bo','filled')
    text(xData(targets.nuc(targets.ordered)),yData(targets.nuc(targets.ordered)),targets.gene(targets.ordered))
    hold on
    xlim([0 25])    
end


clear all
[data,smeta]=getNaamaData('figure',15);

%% create mean Profiles
[expTypes,~,smeta.expId]=unique(smeta(:,[14,7,9,10]),'stable')
for i=unique(smeta.expId(smeta.bad==0))'
    selSmp=smeta.expId==i &smeta.bad==0;
    meanProfile(:,i)=mean(data(:,selSmp),2);
end
GP=load('/home/labs/barkailab/felixj/Documents/MATLAB/scripts/gene/group_imp.mat');
clear data
load('qPCRtargets.mat')
c=0
figure
for i=find(targets.nucPos>0)'
    c=c+1;
    subplot(5,5,c)
    hold off
    if ~contains(targets.name(i),'h2','IgnoreCase',true)
        plot(targets.nucPos(i)+[-500:500]-GP.chrIdx(targets.chr(i)),smoothdata(meanProfile(targets.nucPos(i)+[-500:500],2),1,'gaussian',100))
        hold on
        plot(targets.nucPos(i)+[-500:500]-GP.chrIdx(targets.chr(i)),smoothdata(meanProfile(targets.nucPos(i)+[-500:500],20),1,'gaussian',100))
    else
        plot(targets.nucPos(i)+[-500:500]-GP.chrIdx(targets.chr(i)),smoothdata(meanProfile(targets.nucPos(i)+[-500:500],67),1,'gaussian',100))
        %plot([targets.up(i)-200:targets.down(i)+200]-GP.chrIdx(targets.chr(i)),smoothdata(meanProfile(targets.up(i)-200:targets.down(i)+200,1),1,'gaussian',100))
        hold on
        plot(targets.nucPos(i)+[-500:500]-GP.chrIdx(targets.chr(i)),smoothdata(meanProfile(targets.nucPos(i)+[-500:500],78),1,'gaussian',100))
    end
    %plot([targets.up(i)-200:targets.down(i)+200]-GP.chrIdx(targets.chr(i)),smoothdata(meanProfile(targets.up(i)-200:targets.down(i)+200,10),1,'gaussian',100))
    %xlim([targets.up(i)-100 targets.down(i)+100]-GP.chrIdx(targets.chr(i)))
    title(strjoin(targets{i,[1,9]}))
    axis tight
    if targets.up(i)>0
        xticks(sort([targets.down(i),targets.nucPos(i),targets.up(i)])-GP.chrIdx(targets.chr(i)))
    else
        xticks(targets.nucPos(i)-GP.chrIdx(targets.chr(i)))
    end
    xtickangle(-45)
    %yLim=ylim()
    %plot((targets.up(i)-GP.chrIdx(targets.chr(i)))*[1 1],yLim,'k--')
    %plot((targets.down(i)-GP.chrIdx(targets.chr(i)))*[1 1],yLim,'k--')
    %ylim(yLim)
   %plot((targets.nucPos(i)-GP.chrIdx(targets.chr(i)))*[1 1],yLim,'k-')
   %xlabel(sprintf('Chr %d',targets.chr(i)))
   if targets.peakPos(i)==0
       [x{i} y{i}]=getpts();
   else
       plot(targets.peakPos(i)+[-100 100],mean(ylim)*[1 1],'r-','LineWidth',2)
   end
end
for i=find(targets.peakPos>0)'
    targets.amplicon{i}=SC_genome(targets.chr(i)).Sequence(round(targets.peakPos(i))+[-75:75]);
end
%%
GP.gff(find(contains(GP.gff.type,'centromer')),:)
intExp=[1,15,67,78,7,18]
c=0
for i=find(contains(GP.gff.type,'centromer')&contains(GP.gff.details,'Dbxref'))'
    c=c+1;
    subplot(4,4,c)
    hold off
    pos=GP.gff.gc(i,:);
    for e=1:numel(intExp)
        plot(pos(1)-200:pos(2)+200,meanProfile(pos(1)-200:pos(2)+200,intExp(e)),'-','LineWidth',2,'DisplayName',strjoin(expTypes{intExp(e),:}))
        hold on
    end    
    axis tight
    title(regexp(GP.gff.details(i),'CEN\d+','match','once'))
end
end

%% oligoStand()
clear all
xlsFiles=dir('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/qPCR/*oligo*.xls');
intVars={'SampleName','TargetName','C_','runId'}
for i=1:numel(xlsFiles)
    impOpts=detectImportOptions('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/qPCR/Chipoligos1to6.xls','Sheet','Results');
    temp=readtable([xlsFiles(i).folder,'/',xlsFiles(i).name],impOpts);
    temp.runId=repmat(i,size(temp,1),1)
    intCols=ismember(temp.Properties.VariableNames,intVars)
    samples=~cellfun('isempty',regexp(temp{:,1},'^[A-Z]\d$|^[A-Z]1\d$'));
    resultTable{i}=temp(samples,intCols)    
end
resultTable=cat(1,resultTable{:});
resultTable.n=str2double(regexp(resultTable.SampleName,'\d+','match','once'));
resultTable.targetId=str2double(resultTable.TargetName);
targetNames={'HHF1','HHT2','HTA1','HTA2','HTB2','HYP2','AUS1','ACS1','KOG1','NIT2','RPL3'};

for i=1:max(resultTable.targetId)
    selSmp=find(resultTable.targetId==i);
    subplot(3,4,i)
    hold off
    scatter(log2(resultTable.n(selSmp)),resultTable.C_(selSmp),'o','filled','DisplayName','Data')
    hold on
    title(targetNames{i})
    xlabel('log2(Dilution)')    
    ylabel('Ct Value')
    p=robustfit(log2(resultTable.n(selSmp)),resultTable.C_(selSmp));
    hold on;
    plot([0 11],[0 11]*p(2)+p(1),'k-','DisplayName',sprintf('%.2fx+.%.2f',p(2),p(1)))
    legend('show','Location','best');
end
end

%%  analyseqPCR
clear all
xlsFiles=[dir('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/qPCR/qX*.xls');dir('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/qPCR/qx*.xls')];
intVars={'SampleName','TargetName','C_','runId','Well','Tm1','Tm2','Tm3'}
dnaConc=[0.175600000000000,1.17200000000000,0.194400000000000,0.0216000000000000,0.195200000000000,0.0164000000000000;0.151600000000000,0.184800000000000,0.157600000000000,0,0.179800000000000,0;0.120400000000000,0.0796000000000000,0.116400000000000,0,0.140200000000000,0;0.0504000000000000,0.0344000000000000,0.0646000000000000,0,0.0898000000000000,0;0.0194000000000000,0.0116000000000000,0.0316000000000000,0,0.0234000000000000,0;0,0,0.0124000000000000,0,0.0114000000000000,0;0,0,0,0,0,0;0,0,0,0,0,0];
clear resultTable
for i=1:numel(xlsFiles)
    impOpts=detectImportOptions('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/qPCR/Chipoligos1to6.xls','Sheet','Results');
    temp=readtable([xlsFiles(i).folder,'/',xlsFiles(i).name],impOpts);
    temp.runId=repmat(i,size(temp,1),1)
    intCols=ismember(temp.Properties.VariableNames,intVars)
    samples=~cellfun('isempty',regexp(temp{:,1},'^[A-Z]\d$|^[A-Z]1\d$'));
    resultTable{i}=temp(samples,intCols)
    if contains(xlsFiles(i).name,{'HHF1','HTA1'}) & contains(xlsFiles(i).name,{'qX1'})
        sampleConc=readtable('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/qPCR/samplConc.xlsx','Sheet','HTA1','ReadVariableNames',false);
        sampleConc=sampleConc{:,:};
        [~,rowN]=ismember(regexp(resultTable{i}.Well,'^[A-H]','match','once'),{'A','B','C','D','E','F','G','H'});
        colN=str2double(regexp(resultTable{i}.Well,'\d+','match','once'));
        indN=sub2ind([8,12],rowN,colN);
        resultTable{i}.n=sampleConc(indN);
        resultTable{i}.qX(:)=1;
    elseif contains(xlsFiles(i).name,{'qX2','qX3'},'IgnoreCase',true)
        if contains(xlsFiles(i).name,{'qX2'})
            sampleConc=readtable('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/qPCR/samplConc.xlsx','Sheet','qx2','ReadVariableNames',false);
            if contains(xlsFiles(i).name,'rep2','IgnoreCase',true)
                resultTable{i}.qX(:)=2.5;
            else
                resultTable{i}.qX(:)=2;
            end
        elseif contains(xlsFiles(i).name,{'qX3'},'IgnoreCase',true)
            if contains(xlsFiles(i).name,'qx3f')
                resultTable{i}.qX(:)=3.75;
                sampleConc=readtable('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/qPCR/samplConc.xlsx','Sheet','qx3f','ReadVariableNames',false);               

            elseif contains(xlsFiles(i).name,'asf1')
                resultTable{i}.qX(:)=3.5;
                sampleConc=readtable('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/qPCR/samplConc.xlsx','Sheet','qx3asf1','ReadVariableNames',false);               
            else
                resultTable{i}.qX(:)=3;
                sampleConc=readtable('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/qPCR/samplConc.xlsx','Sheet','qx3','ReadVariableNames',false);
            end
        end
        sampleConc=sampleConc{:,:};
        [~,rowN]=ismember(regexp(resultTable{i}.Well,'^[A-H]','match','once'),{'A','B','C','D','E','F','G','H'});
        colN=str2double(regexp(resultTable{i}.Well,'\d+','match','once'));
        indN=sub2ind([8,12],rowN,colN);
        resultTable{i}.n=sampleConc(indN);
    elseif contains(xlsFiles(i).name,{'qX4'},'IgnoreCase',true)
        resultTable{i}.qX(:)=4;
        sampleConc=readtable('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/qPCR/samplConc.xlsx','Sheet','qx4','ReadVariableNames',false);
        sampleConc=sampleConc{:,:};
        [~,rowN]=ismember(regexp(resultTable{i}.Well,'^[A-H]','match','once'),{'A','B','C','D','E','F','G','H'});
        colN=str2double(regexp(resultTable{i}.Well,'\d+','match','once'));
        indN=sub2ind([8,12],rowN,colN);
        resultTable{i}.n=sampleConc(indN);
    elseif contains(xlsFiles(i).name,{'qX5'},'IgnoreCase',true)
        resultTable{i}.qX(:)=5;
        sampleConc=readtable([xlsFiles(i).folder,'/',xlsFiles(i).name],'Sheet','dilutions','ReadVariableNames',false);
        sampleConc=sampleConc{:,:}; 
        [~,rowN]=ismember(regexp(resultTable{i}.Well,'^[A-H]','match','once'),{'A','B','C','D','E','F','G','H'});
        colN=str2double(regexp(resultTable{i}.Well,'\d+','match','once'));
        indN=sub2ind([8,12],rowN,colN);
        resultTable{i}.n=sampleConc(indN);
        genoType=table2array(readtable([xlsFiles(i).folder,'/',xlsFiles(i).name],'Sheet','Genotype','ReadVariableNames',false));
        abIp=table2array(readtable([xlsFiles(i).folder,'/',xlsFiles(i).name],'Sheet','antibody','ReadVariableNames',false));
        resultTable{i}.SampleName= strcat(abIp(indN),'-',genoType(indN),'-',resultTable{i}.SampleName)
      else
        [~,resultTable{i}.n]=ismember(regexp(resultTable{i}.Well,'^[A-H]','match','once'),{'A','B','C','D','E','F','G','H'})
        resultTable{i}.qX(:)=1;
        sampleConc=readtable('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/qPCR/samplConc.xlsx','Sheet','qx4','ReadVariableNames',false);
    end
    resultTable{i}.TargetName=strcat(resultTable{i}.TargetName,sprintf('_%02d',i));
end
clear impOpts temp intVars indN intCols dnaConc colN i rowN sampleConc samples
resultTable=cat(1,resultTable{:});
resultTable.SampleName=strrep(resultTable.SampleName,'clean','')
allResults=resultTable;

%resultTable=allResults(allResults.qX>=3&~contains(allResults.SampleName,{'DDW','template','Sample 1'}),:)
resultTable=allResults(allResults.qX==5 & contains(allResults.SampleName,'wt'),:)%&contains(allResults.SampleName,{'asf1'}),:)
%resultTable.TargetName(contains(resultTable.Well,{'E','F','G','H'}))=strcat(resultTable.TargetName(contains(resultTable.Well,{'E','F','G','H'})),'B')
resultTable.SampleName(contains(resultTable.Well,{'E','F','G','H'})&resultTable.qX<5)=strcat('Rep2',resultTable.SampleName(contains(resultTable.Well,{'E','F','G','H'})&resultTable.qX<5))

[sample,~,resultTable.sid]=unique(resultTable.SampleName);
chipSample=regexp(resultTable.SampleName,'.*(?=_myc)|.*(?=_HA)|(?<=myc-).*|(?<=HA-).*','match','once');
[chipSample,~,resultTable.cid]=unique(cat(1,chipSample));

dropN={5:7,1:4};
concCol=[5,6,1,2,3,4,1,2];
[targets,~,resultTable.tid]=unique(resultTable.TargetName)
targetColor=lines(numel(targets));
selRun=unique(resultTable.runId);

close all
figure('units','normalized','outerposition',[0 0 1 1])
for i=1:max(resultTable.sid)    
    subplot(3,max(resultTable.cid),resultTable.cid(find(resultTable.sid==i,1))+contains(sample{i},'myc')*numel(chipSample))
    dropi=dropN{1+(contains(sample{i},'myc')&~contains(sample{i},'NC'))};
    hold off
    c=0;
    clear linePlots
    for j=unique(resultTable.tid(resultTable.sid==i&ismember(resultTable.runId,selRun)))'%intersect([1:5],unique(resultTable.tid(resultTable.sid==i)))'
        selWell=resultTable.sid==i & resultTable.tid==j;
        scatter(resultTable.n(selWell),resultTable.C_(selWell),[],targetColor(j,:),'o','filled','DisplayName',targets{j})
        fitWell=selWell & ismember(resultTable.n,dropi);
        if sum(~isnan(resultTable.C_(fitWell)))>=2
            c=c+1;
            if sum(~isnan(resultTable.C_(fitWell)))>2
                p=robustfit(resultTable.n(fitWell),resultTable.C_(fitWell));
            else
                p=fliplr(polyfit(resultTable.n(fitWell&~isnan(resultTable.C_)),resultTable.C_(fitWell&~isnan(resultTable.C_)),1));
            end
            hold on
            linePlots(c)=plot(dropi,dropi.*p(2)+p(1),'k-','Color',targetColor(j,:),'DisplayName',sprintf('%s:%.2fx%+.1f',extractBefore(targets{j},'_'),p(2),p(1)),'LineWidth',2);
            bestFit{i,j}=p;
        else
            bestFit{i,j}=[NaN NaN]
        end
    end
    title(strrep(sample{i},'_',' '))
    rangeY=range(ylim());
    xlabel('dilution step')
    ylabel('Ct')
%     yyaxis right
%     hold off
%     scatter(1:8,-log2(dnaConc(:,concCol(:,i))),'ko','filled','DisplayName','Qubit')
%     ylim(ylim()+[0 rangeY])
%     ylabel('-log2(qubit)')
    if exist('linePlots','var')
        legend(linePlots,'Location','best')
    end
end

for i=1:numel(chipSample)
    subplot(3,numel(chipSample),2*numel(chipSample)+i)
    hold off
    mycId=find(strcmp(sample,[chipSample{i},'_myc']));
    haId=find(strcmp(sample,[chipSample{i},'_HA']));
    for j=unique(resultTable.tid(resultTable.sid==mycId&ismember(resultTable.runId,selRun)))'
        haLevel=[1:7]*bestFit{haId,j}(2)+bestFit{haId,j}(1);        
        mycLevel=[1:7]*bestFit{mycId,j}(2)+bestFit{mycId,j}(1);
        plot([1:7],mycLevel-haLevel,'-o','Color',targetColor(j,:),'DisplayName',targets{j})
        hold on
        text(4,mycLevel(4)-haLevel(4),sprintf('%.2f',mycLevel(4)-haLevel(4)),'Color',targetColor(j,:),'VerticalAlignment','cap','HorizontalAlignment','center')
    end
    title(sprintf('%s vs. %s',sample{mycId},sample{haId}))
    xlabel('Dilution')
    ylabel('\DeltaCt')
end

%% process allResult table
allResults=allResults(~cellfun('isempty',extractBefore(allResults.SampleName,'_')),:)
allResults=allResults(~cellfun('isempty',extractBefore(allResults.SampleName,'_')),:)
allResults.SampleName(contains(allResults.Well,{'E','F','G','H'})&allResults.qX>=3)=strcat('Rep2',allResults.SampleName(contains(allResults.Well,{'E','F','G','H'})&allResults.qX>=3))

%[chipSample,~,allResults.cid]=unique(table(extractBefore(allResults.SampleName,'_'),allResults.qX));
%[sample,~,allResults.sid]=unique(table(allResults.SampleName,allResults.qX));
%[targetR,~,allResults.tid]=unique(table(allResults.TargetName,allResults.runId));
[fitSample,~,allResults.fid]=unique(table(allResults.runId,allResults.TargetName,allResults.SampleName,allResults.qX));
fitSample.Var2=regexp(upper(fitSample.Var2),'[A-Z]+','match','once')
for i=1:size(fitSample,1)
    dropi=dropN{1+(contains(fitSample.Var3{i},'myc')&~contains(fitSample.Var3{i},'NC'))};
    selSample=find(allResults.fid==i & ismember(allResults.n,dropi))
    if sum(~isnan(allResults.C_(selSample)))>2
        pFit{i}=robustfit(allResults.n(selSample),allResults.C_(selSample))
    else
        pFit{i}=[NaN NaN]
    end
    fitSample.Ct(i)=pFit{i}(1)+4*pFit{i}(2);
end
[oligo,~,fitSample.oid]=unique(fitSample.Var2)
intOligos=accumarray(fitSample.oid,1)>2
expColor=lines(3)
c=0
figure
for i=find(intOligos)'
    c=c+1;
    selSample=find(fitSample.oid==i);
    subplot(3,2,c)    
    for j=unique(floor(fitSample.Var4(selSample)))'
        bar(find(floor(fitSample.Var4(selSample))==j),fitSample.Ct(selSample(floor(fitSample.Var4(selSample))==j)),'CData',expColor(j,:));
        hold on
    end
    title(oligo{i})
    ylim([20 35])
    xticks(1:numel(selSample))
    set(gca,'xticklabel',strcat(strrep(fitSample.Var3(selSample),'_',''),'_{',num2str(fitSample.Var1(selSample),'%02g'),'}'),'XTickLabelRotation',90)
end
chipNames=regexp(fitSample.Var3,'.*(?=_myc)|.*(?=_HA)','match');
chipNames=cat(1,chipNames{:});
[deltaCt,~,fitSample.cid]=unique([fitSample(:,[1,2,4,6]),table(chipNames,'VariableNames',{'Chip'})],'stable')
clear chipNames
for i=1:size(deltaCt,1)
    myci=find(fitSample.cid==i&contains(fitSample.Var3,'_myc','IgnoreCase',true));
    hai=find(fitSample.cid==i&contains(fitSample.Var3,'_HA','IgnoreCase',true));
    deltaCt.dCt(i)=fitSample.Ct(myci)-fitSample.Ct(hai);
end
figure
c=0;
for i=find(intOligos)'
    c=c+1;
    selSample=find(deltaCt.oid==i);
    subplot(3,2,c)
    for j=unique(floor(deltaCt.Var4(selSample)))'
        bar(find(floor(deltaCt.Var4(selSample))==j),deltaCt.dCt(selSample(floor(deltaCt.Var4(selSample))==j)),'CData',expColor(j,:));
        hold on
    end
    %bar(1:numel(selSample),deltaCt.dCt(selSample));
    title(oligo{i})
    xticks(1:numel(selSample))
    set(gca,'xticklabel',strcat(strrep(deltaCt.Chip(selSample),'_',''),'_{',num2str(deltaCt.Var1(selSample),'%02g'),'}'),'XTickLabelRotation',90)
    text(1:numel(selSample),deltaCt.dCt(selSample),num2str(deltaCt.dCt(selSample),'%.1f'),'VerticalAlignment','bottom','HorizontalAlignment','center')
end

figure;
c=0;
for i=[1,2,2.5]
    c=c+1;
    subplot(1,3,c)
    selSample=find(deltaCt.Var4==i & contains(deltaCt.Chip,'NC'));
    bar(1:numel(selSample),deltaCt.dCt(selSample))
    text(1:numel(selSample),deltaCt.dCt(selSample),num2str(deltaCt.dCt(selSample),'%.1f'),'VerticalAlignment','bottom','HorizontalAlignment','center')
    set(gca,'XTick',1:numel(selSample),'xticklabel',strcat(deltaCt.Chip(selSample),'-',deltaCt.Var2(selSample)),'XTickLabelRotation',90)
end

end

%% qPCRvsSeq()
clear all
load('./qPCR/resultSummary.mat','deltaCt')
deltaCt.outlier([14,19,20])=1;
load('forFigure15.mat')
load('qPCRtargets.mat')
targets(targets.Var13==0,:)=[];
for i=1:max(floor(deltaCt.Var4))
    deltaCt.adj(floor(deltaCt.Var4)==i)= deltaCt.dCt(floor(deltaCt.Var4)==i) - mean(deltaCt.dCt(floor(deltaCt.Var4)==i & contains(deltaCt.Chip,'NC')))
end
[~,deltaCt.tidx]=ismember(deltaCt.Var2,extractBefore(targets.gene,4));
deltaCt.nuc=targets.nuc(deltaCt.tidx);
sampleConv=readtable('./qPCR/sampleConversion.xlsx','ReadVariableNames',false)
[~,idx]=ismember(strrep(deltaCt.Chip,'Rep2',''),sampleConv.Var1);
deltaCt.string=sampleConv.Var2(idx);
[chipSample,~,expTypes.sid]=unique(expTypes(:,2:4))
[~,deltaCt.sid]=ismember(deltaCt.string,strcat(chipSample.gt,'-',chipSample.tag2,'-',chipSample.con))
for i=1:size(chipSample,1)
    hai=find(expTypes.sid==i & strcmp(expTypes.ab,'ha'));
    myci=find(expTypes.sid==i & strcmp(expTypes.ab,'myc'));
    log2To(:,i)=log2(medNuc(:,myci)+1)-log2(medNuc(:,hai)+1);
end
close all
c=0;
for i=[44,33,45]%unique(deltaCt.sid(deltaCt.sid>0))'%[40,42,33,43,10]%[3,6,7]%unique(deltaCt.sid(deltaCt.sid>0))'
    c=c+1;
    subplot(1,3,c)
    hold off
    selSmp=deltaCt.sid==i &deltaCt.outlier==0;
    for t=unique(deltaCt.tidx(selSmp))'
        scatter(mean(log2To(deltaCt.nuc(selSmp&deltaCt.tidx==t),i)),mean(-deltaCt.adj(selSmp&deltaCt.tidx==t)),[],'o','filled','DisplayName',targets.gene{t})        
        hold on
        errorbar(mean(log2To(deltaCt.nuc(selSmp&deltaCt.tidx==t),i)),mean(-deltaCt.adj(selSmp&deltaCt.tidx==t)),std(-deltaCt.adj(selSmp&deltaCt.tidx==t))/sqrt(sum(selSmp&deltaCt.tidx==t)))
    end
    %text(accumarray(deltaCt.tidx(selSmp),log2To(deltaCt.nuc(selSmp),i),[9,1],@mean,nan),accumarray(deltaCt.tidx(selSmp),-deltaCt.adj(selSmp),[9,1],@mean,nan),targets.gene)
    crVal=corr(log2To(deltaCt.nuc(selSmp),i),-deltaCt.adj(selSmp),'rows','pairwise')
    title(sprintf('%s-%d (%.2f)',strrep(strjoin(chipSample{i,:}),'_',' '),i,crVal))
    %xlabel('log2TO in ChipSeq')
    %ylabel('adj. \Delta\DeltaCt_{qPCR} (ha vs. myc) ~ log2(abs. myc./ha)')
    xlabel('log2(relative HA/myc) (ChIP-Seq)')
    ylabel('log2(absolute HA/myc) (qPCR)')
    p=robustfit(log2To(deltaCt.nuc(selSmp),i),-deltaCt.adj(selSmp));
    Rsquared=corr(-deltaCt.adj(selSmp),p(1)+p(2)*log2To(deltaCt.nuc(selSmp),i))^2
    hold on
    plot(xlim,xlim*p(2)+p(1),'k-','LineWidth',2,'DisplayName',sprintf('%.2fx %+.2f;R^2=%.2f',p(2),p(1),Rsquared))
    hold on;
    if i==44
        plot(xlim,xlim+log2(1.46/316.74),'DisplayName','ChIPSeq+SpikeIn norm')
    end
    axis tight
    legend('show','Location','best')    
    colormap(gca,lines)
end
suptitle('turnover: qPCR vs. Seq')
save_gf(gcf,'qPCRvsSeqPlusSpikeIn')
%% aFactor turnover
clear all
load('forFigure15.mat')
close all
xIdHa=34;
yIdHa=30;
GP=load('group_imp.mat');
load('./extData/addData.mat','h3AdjMid')
xIdMyc=find(ismember(expTypes(:,2:4),expTypes(xIdHa,2:4)) &strcmp(expTypes.ab,'myc'))
yIdMyc=find(ismember(expTypes(:,2:4),expTypes(yIdHa,2:4)) &strcmp(expTypes.ab,'myc'))

xTo=diff(log(medNuc(:,[xIdHa xIdMyc])+1),1,2);
yTo=diff(log(medNuc(:,[yIdHa yIdMyc])+1),1,2);
hold off
nonNan=all(~isnan([xTo,yTo]),2)
%scatter(xTo(fitNucs),yTo(fitNucs),[],0.9*[1 1 1],'.','MarkerEdgeAlpha',.5)
%scatter(xTo(fitNucs),yTo(fitNucs),[],0.9*[1 1 1],'.','MarkerEdgeAlpha',.5)
scatter(xTo(nonNan),yTo(nonNan),[],diff(h3AdjMid(nonNan,:),1,2),'.')
caxis([0 5])
hold on
colormap(gca,parula)
highTo=xTo>log(2);
%scatter(xTo(highTo),yTo(highTo),[],0.7*[1 1 1],'.')

xlabel(['log TurnOver ' strjoin(expTypes{xIdHa,2:4})])
ylabel(['log TurnOver ' strjoin(expTypes{yIdHa,2:4})])
p=polyfit(xTo(highTo),yTo(highTo),1);
hold on
plot(xlim,xlim()*p(1)+p(2),'-','LineWidth',2,'DisplayName',sprintf('%.2f x%+.2f : %.2f',p(1),p(2),corr(xTo(highTo),yTo(highTo))))
%plot(xlim,xlim,'k--','LineWidth',2,'DisplayName','1:1')
matProms=ismember(nucMeta.gene,GP.groups{8}{2}{9}) & nucMeta.order<0;
hisNucs=find(nucMeta.order<0 & ismember(nucMeta.gene,GP.groups{23}{2}{45}));
scatter(xTo(matProms),yTo(matProms),'o','filled','DisplayName','Mating gene promoters')
%scatter(xTo(hisNucs),yTo(hisNucs),'o','filled','DisplayName','Histone gene promoters')
crVal=corr(xTo(nonNan),yTo(nonNan));
title(crVal)
ylabel(colorbar(),'estimated non Sphase exch. events')
save_gf(gcf,'asyncVsAlphaCmap')

normFac=2^(-9.71+7.04)
p=linortfit2(medNuc(highTo,yIdMyc),diff(h3AdjMid(highTo,:),1,2));
normFac=p(1)
close all
figure

haBins=[7.5:1:16.5]%[6.5:1:14.5]%[0.7:0.1:1.6];
xVal=medNuc(:,yIdHa)./median(medNuc(:,yIdHa),'omitnan');
yVal=medNuc(:,yIdMyc).*normFac;

scatter(xVal,yVal,[],diff(h3AdjMid,1,2),'o','filled')
caxis([0.5 5])
ylabel(colorbar(),'estimated non Sphase exch. events')
xlabel(strjoin(expTypes{yIdHa,[1:4]}))
ylabel(strjoin(expTypes{yIdMyc,[1:4]}))
title(corr(xVal,yVal,'rows','pairwise'))
haBins=[7.5:1:16.5]/median(medNuc(:,yIdHa),'omitnan');
nucBin=sum(xVal>(haBins),2);
lineNucs=nucBin>0 & nucBin<numel(haBins);
medLineX=accumarray(nucBin(lineNucs),xVal(lineNucs),[],@median);
medLineY=accumarray(nucBin(lineNucs),yVal(lineNucs),[],@(x)quantile(x,0.475));
pMid=fliplr(polyfit(medLineX,medLineY,1));
hold on
baseLine=plot(xlim(),xlim()*pMid(2)+pMid(1),'k-','LineWidth',2,'DisplayName',sprintf('BaseLine: %.2fx%+.2f',pMid(2),pMid(1)));
hold on
scatter(xVal(intNucs),yVal(intNucs),'ro','filled')
save_gf(gcf,'alphaMcyvsHaWBaseline')

figure
scatter(medNuc(:,yIdHa),medNuc(:,yIdMyc).*normFac,[],diff(h3AdjMid,1,2),'o','filled')
caxis([0.5 5])
ylabel(colorbar(),'estimated non Sphase exch. events')
xlabel(strjoin(expTypes{yIdHa,[1:4]}))
ylabel(strjoin(expTypes{yIdMyc,[1:4]}))
title(corr(medNuc(:,yIdHa),medNuc(:,yIdMyc),'rows','pairwise'))
save_gf(gcf,'alphaMcyvsHa')

%%HHT2 vs
cId=34
cId2=34+54;
c=0;
for xId=[52,39]
    c=c+1;
    subplot(1,2,c)
    yId=find(ismember(expTypes(:,2:4),expTypes(xId,2:4)) & contains(expTypes.ab,'myc'));
    scatter(medNuc(:,xId),medNuc(:,yId),[],medNuc(:,cId),'.')
    caxis([0 20])
    xlabel(sprintf('%s : %0.2f',strjoin(expTypes{xId,[1:4]}),corr(medNuc(:,xId),medNuc(:,cId),'rows','pairwise')))
    ylabel(sprintf('%s : %0.2f',strjoin(expTypes{yId,[1:4]}),corr(medNuc(:,yId),medNuc(:,cId2),'rows','pairwise')))

    ylabel(colorbar(),strjoin(expTypes{cId,[1:4]}))
end
save_gf(gcf,'hht2-h3x2')
%% H2Aint vs. H2B
clearvars()
xIdHa=34;
yIdHa=30;
GP=load('group_imp.mat');

xIdMyc=find(ismember(expTypes(:,2:4),expTypes(xIdHa,2:4)) &strcmp(expTypes.ab,'myc'))
yIdMyc=find(ismember(expTypes(:,2:4),expTypes(yIdHa,2:4)) &strcmp(expTypes.ab,'myc'))

xTo=diff(log2(medNuc(:,[xIdHa xIdMyc])+1),1,2);
yTo=diff(log2(medNuc(:,[yIdHa yIdMyc])+1),1,2);

%% compare Histone occupancy different mutants
clear all
load('forFigure15.mat')
GP=load('group_imp.mat');
load('qPCRtargets.mat')

intHa=[17,16,31,30]
histoneNucs=find(ismember(nucMeta.gene,GP.groups{23}{2}{45}) & nucMeta.order<0);
intHisNucs=intersect(targets.nuc(targets.Var13==1),histoneNucs);
c=0;
for i=intHa
    c=c+1;
    subplot(2,2,c)
    hold off
    [y,x]=histcounts(medNuc(:,i),20)
    plot(movmean(x,2,'Endpoints','discard'),y,'LineWidth',2)
    hold on
    plot(medNuc(histoneNucs,i)'.*[1;1],repmat(ylim',1,numel(histoneNucs)),'DisplayName','Histone Nucs','Color',[1 0 0 ])
    xlabel(strjoin(expTypes{i,1:4}))    
    ylabel('# nucleosomes')
    text(medNuc(intHisNucs,i),repmat(mean(ylim),numel(intHisNucs),1),GP.gene_infoR64.nameNew(nucMeta.gene(intHisNucs)),'HorizontalAlignment','center')
end

%%  test()
xId=29
yId=17
hisNucs=ismember(nucMeta.gene,GP.groups{23}{2}{45}) & nucMeta.order<=0;

subplot(1,2,1)
scatter(medNuc(:,xId+52),medNuc(:,yId+52),'.');
xlabel(strjoin(expTypes{xId+52,1:4}))
ylabel(strjoin(expTypes{yId+52,1:4}))
hold on
plot(xlim,xlim,'k--')
scatter(medNuc(hisNucs,xId+52),medNuc(hisNucs,yId+52),'ro','filled')

subplot(1,2,2)
xTo=log(medNuc(:,xId+52)+1)-log(medNuc(:,xId)+1);%medNuc(:,xId+52)%
yTo=log(medNuc(:,yId+52)+1)-log(medNuc(:,yId)+1);%medNuc(:,yId+52)%
scatter(xTo,yTo,'.');
xlabel(strjoin(expTypes{xId,1:4}))
ylabel(strjoin(expTypes{yId,1:4}))
hold on
plot(xlim,xlim,'k--')
scatter(xTo(hisNucs),yTo(hisNucs),'ro','filled')


%% look at antisense transcription
clear all
plusWig=readtable('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/extData/GSM617027_WT_NC_plus.wig','FileType','text','ReadVariableNames',false);
negWig=readtable('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/extData/GSM617027_WT_NC_minus.wig','FileType','text','ReadVariableNames',false);
GP=load('group_imp.mat');
chrNames={'chrI'	'chrII'	'chrIII'	'chrIV'	'chrV'	'chrVI'	'chrVII'	'chrVIII'	'chrIX'	'chrX'	'chrXI'	'chrXII'	'chrXIII'	'chrXIV'	'chrXV'	'chrXVI' 'chrMt'};

for i=1:numel(chrNames)
    occCol{i}=zeros(GP.chr_len(i),2);
    plusTableStart=find(endsWith(plusWig.Var2,chrNames{i}));
    if numel(plusTableStart)>0
        plusTableEnd=min([find(contains(plusWig.Var2,'chr') & ([1:size(plusWig,1)]'>plusTableStart),1),size(plusWig,1)]);
        occCol{i}(str2double(plusWig.Var1(plusTableStart+1: plusTableEnd-1)),1)=str2double(plusWig.Var2(plusTableStart+1: plusTableEnd-1));
    end
    minusTableStart=find(endsWith(negWig.Var2,chrNames{i}));
    if numel(minusTableStart)>0
        minusTableEnd=min([find(contains(negWig.Var2,'chr') & ([1:size(negWig,1)]'>minusTableStart),1),size(negWig,1)]);
        occCol{i}(str2double(negWig.Var1(minusTableStart+1: minusTableEnd-1)),2)=str2double(negWig.Var2(minusTableStart+1: minusTableEnd-1));
    end
end
regev=load ('/home/labs/barkailab/LAB/data/DataExternal/EXPRESSION/regev_pnas_2008_rnaseq/data');
geneLvl=median(regev.regev_rnaseq,2,'omitnan');
occCol=cat(1,occCol{:});
clearvars -except GP occCol
tssIdx=nan(6701,1);
hasTss=~isnan(GP.gene_infoR64.felixTss(:,1));
tssIdx(hasTss)=GP.gene_infoR64.felixTss(hasTss,2)+GP.chrIdx(GP.gene_infoR64.felixTss(hasTss,1));
tssSur=nan(6701,601,2);
posTss=hasTss&GP.gene_infoR64.dir>0;
tssSur(posTss,:,:)=reshape(occCol(acol(tssIdx(posTss)+[-300:300]),[1,2]),sum(posTss),601,2);
negTss=hasTss&GP.gene_infoR64.dir<0;
tssSur(negTss,:,:)=reshape(occCol(acol(tssIdx(negTss)+[300:-1:-300]),[2,1]),sum(negTss),601,2);
regev=load ('/home/labs/barkailab/LAB/data/DataExternal/EXPRESSION/regev_pnas_2008_rnaseq/data');
geneLvl=median(regev.regev_rnaseq,2,'omitnan');
clear regev
[~,idx]=sort(geneLvl.*(tssIdx./(tssIdx)),'descend','MissingPlacement','last')
figure
subplot(2,2,1)
imagesc(tssSur(idx,:,1),'XData',[-300 300],[0 20])
ylabel('Genes ordered by ExpLevel')
xlabel('Aligned by TSS')
colormap(gca,brewermap(64,'OrRd'))
ylabel(colorbar(),'Occupancy in sense direction')
title('sense transcription')
subplot(2,2,2)
imagesc(tssSur(idx,:,2),'XData',[-300 300],[0 5])
ylabel('Genes ordered by ExpLevel')
xlabel('Aligned by TSS')
colormap(gca,brewermap(64,'Blues'))
ylabel(colorbar(),'Occupancy in sense direction')
title('antisense transcription')
tsScore=squeeze(sum(tssSur(:,301:end,:),2));
subplot(2,2,3)
histogram(log2(tsScore(:,2)+.1),50)
ylabel('# genes')
xlabel('# antisense transcripts')
subplot(2,2,4)
scatter(log2(tsScore(:,1)+1),log2(tsScore(:,2)+0.1),'.')
[~,intGenes]=ismember({'HMS2','PHO5','GAL80','PHO84','GAL1'},GP.gene_infoR64.name);
hold on
scatter(log2(tsScore(intGenes,1)+1),log2(tsScore(intGenes,2)+0.1),'filled')
xlabel('sense transcription')
ylabel('antisense transcription')


%% look at TATA and SAGA genes
clear all
GP=load('group_imp.mat');
saga=detectImportOptions('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/extData/SAGAPugh2004.xls','VariableNamesRange','A2')
saga=readtable('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/extData/SAGAPugh2004.xls',saga);
tata=readtable('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/extData/TATAGenesPugh2004.xls');
[~,tata.barkaiId]=ismember(tata.Gene,GP.gene_infoR64.orf);
[~,saga.barkaiId]=ismember(saga.GeneID,GP.gene_infoR64.orf);

[sagaType,~,saga.sagaId]=unique(saga.SAGA_TFIIDGroup)
load('holstege.mat')
geneStdH=std(holstege.TsXMut,[],2,'omitnan');
meanChg=median(holstege.TsXMut,2,'omitnan');
figure;
[yAll,bins]=histcounts(geneStdH,'Normalization','probability')
subplot(1,2,1)
hold off
plot(movmean(bins,2,'Endpoints','discard'),yAll,'DisplayName','All genes','LineWidth',2)
hold on
tataType([1,2,3,9])={'TATA','Orth','TBP','None'};
for i=unique(tata.Subclassification_seeFig_S1_(tata.Subclassification_seeFig_S1_>0))'
    [yi]=histcounts(geneStdH(tata.barkaiId(tata.Subclassification_seeFig_S1_==i & tata.barkaiId>0)),bins,'Normalization','probability');
    plot(movmean(bins,2,'Endpoints','discard'),movmean(yi,4),'DisplayName',sprintf('SubClass:%s',tataType{i}),'LineWidth',2)
end
title('TATA genes vs plasticity')
xlabel('abs(flexibility)')
ylabel('fraction')
axis tight
subplot(1,2,2)
hold off
plot(movmean(bins,2,'Endpoints','discard'),yAll,'DisplayName','All genes','LineWidth',2)
hold on
for i=2:4 % only assigned ones are shown
    [yi]=histcounts(geneStdH(saga.barkaiId(saga.sagaId==i & saga.barkaiId>0)),bins,'Normalization','probability');
    plot(movmean(bins,2,'Endpoints','discard'),movmean(yi,4),'DisplayName',sagaType{i},'LineWidth',2)
end
title('SAGA genes vs plasticity')
xlabel('abs(flexibility)')
ylabel('fraction')
axis tight
save_gf(gcf,'SAGA-TATAvsFlex')
load('forFigure15.mat')
close all
figure
intExp=[69,15,88,34];
expColor([6,22,7,23],:) =[100,100,100;127,63,152;150,150,150;0,174,239]/256 
load('greytopurple.mat','Colo')
load('greytoblue.mat','Colo1')
c0=0
subPos=[1,5,2,6,3,7,4,8]
for e=intExp;
    clear imageMat
    c0=c0+1;
    c1=0;
    subplot(4,4,subPos(c0))
    for n=[-2,-1,1,2]
        c1=c1+1;
        c2=0;
        for t=[1:3,9]
            c2=c2+1;
            selNucs=ismember(nucMeta.gene,tata.barkaiId(tata.Subclassification_seeFig_S1_==t)) & nucMeta.order==n;
            imageMat(c2,c1)=median(medNuc(selNucs,e),'omitnan');
        end
    end
    imagesc(imageMat)
    ylabel(colorbar(),'median level')
    xlabel('nuc position')
    if c0==1
            ylabel('TATA class')   
            yticks(1:4)
            yticklabels({'TATA','TATA orth','TBP resp','none'})
    end
    xticks(1:4)
    xticklabels([-2,-1,1,2])
    if contains(expTypes.tag2(e),'h2')
        colormap(gca,Colo)
    else
        colormap(gca,Colo1)
    end
    %cMap=arrayfun(@(x)brighten(expColor(e,:),x),[1:-.01:0]','UniformOutput',false);
    %colormap(gca,cat(1,cMap{:}))
    title(strjoin(expTypes{e,1:4}))
end

for e=intExp;
    clear imageMat
    c0=c0+1;
    c1=0;
    subplot(4,4,subPos(c0))
    for n=[-2,-1,1,2]
        c1=c1+1;
        c2=0;
        for t=[2:4]
            c2=c2+1;
            selNucs=ismember(nucMeta.gene,saga.barkaiId(saga.sagaId==t)) & nucMeta.order==n;
            imageMat(c2,c1)=median(medNuc(selNucs,e),'omitnan');
        end
    end
    imagesc(imageMat)
    ylabel(colorbar(),'median level')
    if c0==5
        ylabel('SAGA class')
        yticks(1:4)
        yticklabels(sagaType(2:4))
    end    
    if contains(expTypes.tag2(e),'h2')
        colormap(gca,Colo)
    else
        colormap(gca,Colo1)
    end
    xticks(1:4)
    xticklabels([-2,-1,1,2])
    title(strjoin(expTypes{e,1:4}))
end
save_gf(gcf,'SAGA-TATATurnOver')
clear imageMat
figure
e=81;
c1=8;
for n=[-2,-1,1,2]
    c1=c1+1;
    subplot(4,4,c1)
    hold off
    c2=0;
    for t=[1:3,9]
        c2=c2+1;
        selNucs=ismember(nucMeta.gene,setdiff(tata.barkaiId(tata.Subclassification_seeFig_S1_==t),GP.groups{23}{2}{45})) & nucMeta.order==n;
        %selNucs=ismember(nucMeta.gene,setdiff(saga.barkaiId(saga.sagaId==t),GP.groups{23}{2}{45})) & nucMeta.order==n;
        %selNucs= nucMeta.order==n;
        scatter(repmat(c2,sum(selNucs),1)+randn(sum(selNucs),1)*.1,medNuc(selNucs,e),[],geneStdH(nucMeta.gene(selNucs)),'.');
        hold on
        scatter(c2,median(medNuc(selNucs,e),'omitnan'),'ko','filled')
        %imageMat(c2,c1)=median(medNuc(selNucs,23),'omitnan');
        forXlabel{c2}=sprintf('%s\\newline%.2f',tataType{t}',corr(medNuc(selNucs,e),geneStdH(nucMeta.gene(selNucs)),'rows','pairwise'))
        caxis([0 0.2])
    end
    xticks(1:c2)
    xticklabels(forXlabel)
    set(gca,'Color',[.9 .9 .9])
    title(sprintf('Nuc: %d',n))
    ylabel(strjoin(expTypes{e,1:4}))
end
c1=12;
for n=[-2,-1,1,2]
    c1=c1+1;
    subplot(4,4,c1)
    hold off
    c2=0;
    for t=[2:4]
        c2=c2+1;
        %selNucs=ismember(nucMeta.gene,setdiff(tata.barkaiId(tata.Subclassification_seeFig_S1_==t),GP.groups{23}{2}{45})) & nucMeta.order==n;
        selNucs=ismember(nucMeta.gene,setdiff(saga.barkaiId(saga.sagaId==t),GP.groups{23}{2}{45})) & nucMeta.order==n;
        %selNucs= nucMeta.order==n;
        scatter(repmat(c2,sum(selNucs),1)+randn(sum(selNucs),1)*.1,medNuc(selNucs,e),[],geneStdH(nucMeta.gene(selNucs)),'.');
        hold on
        scatter(c2,median(medNuc(selNucs,e),'omitnan'),'ko','filled')
        %imageMat(c2,c1)=median(medNuc(selNucs,23),'omitnan');
        forXlabel{c2}=sprintf('%s\\newline%.2f',sagaType{t}',corr(medNuc(selNucs,e),geneStdH(nucMeta.gene(selNucs)),'rows','pairwise'))
        caxis([0 0.2])
    end
    xticks(1:c2)
    xticklabels(forXlabel)
    set(gca,'Color',[.9 .9 .9])
    title(sprintf('Nuc: %d',n))
    ylabel(strjoin(expTypes{e,1:4}))
end

%% compare saga to AntiSense

[yAll,bins]=histcounts(log2(tsScore(:,2)),'Normalization','probability')
subplot(1,2,1)
hold off
plot(movmean(bins,2,'Endpoints','discard'),yAll,'DisplayName','All genes','Color',[.5 .5 .5],'LineWidth',2)
hold on
tataType([1,2,3,9])={'TATA','Orth','TBP','None'};
for i=unique(tata.Subclassification_seeFig_S1_(tata.Subclassification_seeFig_S1_>0))'
    [yi]=histcounts(log2(tsScore(tata.barkaiId(tata.Subclassification_seeFig_S1_==i & tata.barkaiId>0),2)),bins,'Normalization','probability');
    plot(movmean(bins,2,'Endpoints','discard'),movmean(yi,4),'DisplayName',sprintf('SubClass:%s',tataType{i}),'LineWidth',2)
end
title('TATA genes vs plasticity')
xlabel('abs(flexibility)')
ylabel('fraction')
axis tight
subplot(1,2,2)
hold off
plot(movmean(bins,2,'Endpoints','discard'),yAll,'DisplayName','All genes','Color',[.5 .5 .5],'LineWidth',2)
hold on
for i=2:4 % only assigned ones are shown
    [yi]=histcounts(log2(tsScore(saga.barkaiId(saga.sagaId==i & saga.barkaiId>0),2)),bins,'Normalization','probability');
    plot(movmean(bins,2,'Endpoints','discard'),movmean(yi,4),'DisplayName',sagaType{i},'LineWidth',2)
end
title('SAGA genes vs plasticity')


%% saga profile
clear all
[data,smeta]=getNaamaData('figure',16);
data(:,~ismember(smeta.ab,{'ha','myc'}))=[];
smeta(~ismember(smeta.ab,{'ha','myc'}),:)=[];
[smeta,idx]=sortrows(smeta,{'ab','tag2','bad'});
data=data(:,idx);
load('addData.mat','saga','sagaType','tata','tataType')

[~,~,smeta.sid]=unique(regexp(smeta.name,'(?<=(^myc-|^ha-)).*(?=\.mat)','match','once'),'stable')
for i=1:max(smeta.sid)        
    myci=find(strcmp(smeta.ab,'myc') & smeta.sid==i);
    hai=find(strcmp(smeta.ab,'ha') & smeta.sid==i);
    if numel([myci,hai])==2
        sampleTo(:,i)=log(data(:,myci)+1)-log(data(:,hai)+1);        
    end
end
load('cyclebase.mat','cyclebase')
GP=load('group_imp.mat');
geneGroups={1:6701,GP.groups{6}{2}{1},find(cyclebase.data(:,1)<200),GP.groups{7}{2}{66},...
    tata.barkaiId(tata.barkaiId>0&tata.Subclassification_seeFig_S1_==1),saga.barkaiId(saga.barkaiId>0 & saga.sagaId==2)}
groupNames={'All gene*','ESR ind','CC Top 200','RiBi','TATA','SAGA'}
intTag={'h3','h4','h2a','h2b','h2Int'}

metaProfile=meta_profile(sampleTo,'promoter',600,'useORF',false,'afterTSS',700,'scaled',false,'disCDS',true);
cLim=[0,.9;0 .5;0 .5;0 .5;0 .5];
close all
figure;
for t=1:numel(intTag)
    selSmp=smeta.sid(smeta.bad==0 & strcmp(smeta.tag2,intTag{t}) & strcmp(smeta.ab,'myc'));
    if t==numel(intTag)
        selSMp=find(selSmp);
        selSmp=flipud(selSmp([4,1,2,3,5,6]));
    end
    for i=1:numel(geneGroups)
        subplot(numel(intTag),numel(geneGroups),(t-1)*numel(geneGroups)+i)
        imMat=permute(median(metaProfile.cube(geneGroups{i},:,selSmp),1,'omitnan'),[3,2,1]);
        imagesc(imMat,'XData',[-600 700])
        yticks([])
        xlim([-600 600])
        xticks([-500,0,500])
        if t==1
            title(groupNames{i})
        end
        ylabel([intTag{t} '-logTO'] )
        colorbar()
        caxis(cLim(t,:))
    end
end
save_gf(gcf,'Fig1_groupProfileAllwH2Int_6SamplesOrderHightoLow')

metaProfile=meta_profile(data,'promoter',600,'useORF',false,'afterTSS',700,'scaled',false,'disCDS',true);
cLim={[0,.9;0 .5;0 .5;0 .5;0 .5];[]}
intAb={'ha','myc'}
for ab=1:numel(intAb)
    figure;
    for t=1:numel(intTag)
        selSmp=find(smeta.bad==0 & strcmp(smeta.tag2,intTag{t}) & strcmp(smeta.ab,intAb{ab}));
        for i=1:numel(geneGroups)
            subplot(numel(intTag),numel(geneGroups),(t-1)*numel(geneGroups)+i)
            imMat=permute(median(metaProfile.cube(geneGroups{i},:,selSmp),1,'omitnan'),[3,2,1]);
            imagesc(imMat,'XData',[-600 700])
            yticks([])
            xticks([-500,0,500])
            if t==1
                title(groupNames{i})
            end
            ylabel([intTag{t} '-' intAb{ab}] )
            colorbar()
            %caxis(cLim{ab}(t,:))
        end
    end
end

%% correlation median turnover
clear all 
[expTypes,medNuc,smeta,nucData,nucMeta]=standardImport(16,'save',true);

load('greytopurple.mat','Colo')
for i=find(strcmp(expTypes.ab,'ha'))'
    myci=find(strcmp(expTypes.ab,'myc')&ismember(expTypes(:,2:4),expTypes(i,2:4)))
    tempLine=expTypes(i,1:5);
    tempLine.ab={'logTo'};
    medNuc(:,end+1)=diff(log(medNuc(:,[i myci])+1),1,2);
    expTypes(end+1,:)=tempLine;
end
intExp=[39,43,37,38,36]
crMat=corr(medNuc(:,intExp),'rows','complete');
figure;
hold off
imagesc(crMat,[.39 1])
yticks(1:numel(intExp))
yticklabels(expTypes.tag2(intExp))
xticks(1:numel(intExp))
xticklabels(expTypes.tag2(intExp))
ylabel(colorbar('west'),'all nuc TO correlation')
title('turnover correlation voe rall nucleosomes')
colormap(gca,Colo)
save_gf(gcf,'/Fig1_medTOcorrelations')


%% compare meta data between all samples
clear all
dataFolders=readtable('dataFolders.xlsx','ReadVariableNames',false)
innerBCs=readtable('/home/labs/barkailab/LAB/scripts/felix_bcl2fastq/dna_barcodes.txt','ReadVariableNames',false);
innerBCs=innerBCs(:,[5 6]);
outerBCs=table({'Enr1','Enr2','Enr3','Enr4'}',{'TCTACTCT','CTCCTTAC','TATGCAGT','TACTCCTT'}')

for i=1:size(dataFolders,1)
    dnaIndxFile=dir([dataFolders.Var1{i},'/*DNA*.xlsx']);
    clear indexList
    for j=1:numel(dnaIndxFile)
        try
            indexList{j}=readtable([dnaIndxFile(j).folder '/' dnaIndxFile(j).name],'ReadVariableNames',false);
        end
    end
    indexList=cat(1,indexList{:});
    metaDataFile=dir([dataFolders.Var1{i},'/*.tsv']);
    [~,selMeta]=max([metaDataFile.bytes]);
    metaData=readtable([metaDataFile(selMeta).folder '/' metaDataFile(selMeta).name],'FileType','text');
    metaData.Properties.VariableNames=strrep(metaData.Properties.VariableNames,'Sample','File');
    if any(contains(metaData.File,'Enr'))
        [~,enrIdx]=ismember(regexp(metaData.File,'Enr\d','match','once'),outerBCs.Var1);
        [~,bcIdx]=ismember(regexp(metaData.File,'(?<=Enr\d_)\d+[A-Z]','match','once'),innerBCs.Var5);
        metaData.barCode=strcat(outerBCs.Var2(enrIdx),'_',innerBCs.Var6(bcIdx)) ;
    else
        metaData.barCode=extractBefore(metaData.File,'_S');
    end
    [~,idx]=ismember(indexList.Var1,metaData.barCode);
    metaColl{i}=[metaData(idx(idx>0),end-1:end),indexList(idx>0,:)];
end

%% X34 and X33
[expTypes,medNuc,smeta,nucData,nucMeta]=standardImport(18)
GP=load('group_imp.mat');
close all
figure
c=0
histoneNucs=find(ismember(nucMeta.gene,GP.groups{23}{2}{45}) & nucMeta.order<0)
for i=[63,58,52:57,59:61]%[51,47,41:46,48:50]%
    c=c+1;
    subplot(3,4,c)
    hold off
    scatter(medNuc(:,62),medNuc(:,i),[],nucMeta.order,'.','DisplayName','All nucs')
    xlabel(strjoin(expTypes{62,1:4}))
    ylabel(strjoin(expTypes{i,1:4}))
    hold on
    scatter(medNuc(histoneNucs,62),medNuc(histoneNucs,i),'ro','filled','DisplayName','histone promoters')
    caxis([-2 2])    
    set(gca,'Color',.8*[1 1 1])
end
ylabel(colorbar(),'nucleosome position')
for i=1:4
    subplot(2,2,i)
    hold off
    scatter(cerpar(i*4-3:i*4,1),cerpar(i*4-3:i*4,2),'o','filled')
    hold on
    p=robustfit(cerpar(i*4-2:i*4,1),cerpar(i*4-2:i*4,2));
    plot([4 20],[4 20]*p(2)+p(1),'k-','DisplayName',sprintf('%.2fx %+.2f',p(2)*100,p(1)*100))
    legend()
    xlabel('Scer:Par Volume')
    ylabel('Scer:Par Sequencing')
    legend()
end


%% heatshock experiment
clear all
%load('forFigure15.mat')
%h3To=diff(log(medNuc(:,[34 34+54])+1),1,2);
%[expTypes,medNuc,smeta,nucData,nucMeta]=standardImport(17,'save',true);
load('forFigure17.mat')
clearvars -except smeta nucData nucMeta
load('./extData/addData.mat','h3To')
smeta.tp=str2double(regexp(smeta.con,'(?<=t|T)\d+','match','once'));
smeta.tp(isnan(smeta.tp))=-1;
[tcs,~,smeta.tcId]=unique(smeta(:,[7,9,14]),'rows');
[smeta,idx]=sortrows(smeta,{'tcId','tp'})
nucData=nucData(:,idx);
clear idx
nucDelta=nan(size(nucData));
for t=1:numel(tcs)
    temp=log(nucData(:,smeta.tcId==t)+1);
    nucDelta(:,smeta.tcId==t)=temp-median(temp,2);
end
GP=load('group_imp.mat');
expHS=load('./extData/HS_wt.mat')
dynExp=expHS.HS_SD.exp-median(expHS.HS_SD.exp,2,'omitnan');
load('extData/addData.mat','tata','geneLvl')
indGenes=find(~(mean(dynExp(:,4:6),2)<-.5))
hsf1=load('ChecProfile.mat','Hsf1_DK');
checProfiles=hsf1.Hsf1_DK;
clear hsf1
metaProfile=meta_profile(checProfiles,'promoter',700,'useORF',false,'afterTSS',100,'scaled',false,'disCDS',false);
sumProm=permute(sum(metaProfile.cube(:,1:700,:),2),[1,3,2]);
hsf1Targets=sumProm>20000;


selGenes={GP.groups{7}{2}{66}(:),GP.groups{6}{2}{1}(:),find(hsf1Targets),GP.groups{23}{2}{45}(:),...
    tata.barkaiId(tata.Subclassification_seeFig_S1_==1&(tata.barkaiId>0)&~ismember(tata.barkaiId,[indGenes;GP.groups{6}{2}{1}(:)]) ),...
    tata.barkaiId(tata.Subclassification_seeFig_S1_==1&(tata.barkaiId>0))}
selNames={'RiBi','Stress','HSF1','Histones','TATA','TATA full'}
intNuc={[-2,-1],[2,3]};

% by time course
figure
for g=1:numel(selGenes)
    for n=1:numel(intNuc)
        selNucs=ismember(nucMeta.order,intNuc{n}) & ismember(nucMeta.gene,selGenes{g});
        for t=[1:4]    
            subplot(4,2,t*2-2+n)
            plot(smeta.tp(smeta.tcId==t),median(nucDelta(selNucs,smeta.tcId==t),'omitnan'),'LineWidth',2,'DisplayName',selNames{g})
            hold on
            title(sprintf('%s;nuc%d',strjoin(tcs{t,:}),intNuc{n}(1)))
            legend()
            xlabel('time after HS')
            ylabel('median log \delta')
        end        
    end
end
close all
figure
subPlot=[2,4];
c=0
for g=1:numel(selGenes)
    for n=1:numel(intNuc)
        selNucs=ismember(nucMeta.order,intNuc{n}) & ismember(nucMeta.gene,selGenes{g});
        c=c+1;
        %subplot(3,2,subPlot(c))
        hold off
        subplot(6,2,g*2-2+n)
        for t=[3,4]
            lineName=sprintf('%s',strjoin(tcs{t,:}));
            temp=median(nucDelta(selNucs,smeta.tcId==t),'omitnan');
            temp=smoothdata(temp,2,'gaussian',3)
            temp=temp-temp(:,1);
            plot(smeta.tp(smeta.tcId==t),temp,'LineWidth',2,'DisplayName',lineName)
            hold on
            legend()
        end  
        title(sprintf('%s; nuc:%d',selNames{g},intNuc{n}(1)))            
        xlabel('time after HS')
        ylabel('median relative \Deltamyc (log)')
        axis tight
    end
end
subplot(1,2,1)
mycIdx=[find(ismember(smeta.tcId,[4]));find(ismember(smeta.tcId,[3]))];
imagesc(corr(nucData(:,mycIdx),'rows','pairwise'))
caxis([0 1])
colorbar()
yticks(1:numel(mycIdx))
yticklabels(strcat(smeta.gt(mycIdx),'-T:',num2str(smeta.tp(mycIdx))))


for t=1:size(tcs,1)
    tcsMed(:,t)=median(nucData(:,smeta.tcId==t&smeta.bad==0),2,'omitnan');
end

close all
for t=[3,4]
    figure;
    c=0
    for s=find(smeta.tcId==t)'
        c=c+1;
        subplot(3,4,c)
        scatter(tcsMed(:,t-2),nucData(:,s),[],h3To,'.')
        xlabel('median HA level in TC')
        ylabel(smeta.name(s))
        caxis([-1 2.5])
        title(corr(tcsMed(:,t-2),nucData(:,s),'rows','pairwise'))
    end
    ylabel(colorbar(),'logTo in H3 wt cells')
    suptitle(tcs.gt{t})
end

%%%;
clear all;
[expTypes,medNuc,smeta,nucData,nucMeta]=standardImport(18);
figure
imagesc(corr(medNuc,'rows','pairwise'))
yticks(1:size(medNuc,2))
yticklabels(strcat(expTypes.tag2,';',expTypes.con))


xticks(1:size(medNuc,2))
xticklabels(strcat(expTypes.tag2,';',expTypes.con))
xtickangle(90)

%% spike ins
bamFiles=dir('/home/labs/barkailab/LAB/data/SEQ/Gilad_X34Spike_210218/spikeLowTh/temp/*.bam');
bamFiles=rmfield(bamFiles,{'isdir','bytes','date','datenum'});
chrNames={'NC_001133.9','NC_001134.8','NC_001135.5','NC_001136.10','NC_001137.3','NC_001138.5',...
    'NC_001139.9','NC_001140.6','NC_001141.2','NC_001142.9','NC_001143.9','NC_001144.5','NC_001145.3',...
    'NC_001146.8','NC_001147.6','NC_001148.4','NC_001224.1','chrI','chrII','chrIII','chrIV','chrV','chrVI','chrVII','chrVIII','chrIX',...
    'chrX','chrXI','chrXII','chrXIII','chrXIV','chrXV','chrXVI','CBS432_rDNA:1-8984:+'}
samCommand= 'samtools view -f2 -F0x400 %s | grep "AS:i:0.*YS:i:0" | grep -v "XS:" | cut -f3 | uniq -c'
% get names
innerBCs=readtable('/home/labs/barkailab/LAB/scripts/felix_bcl2fastq/dna_barcodes.txt','ReadVariableNames',false);
innerBCs=innerBCs(:,[5 6]);
outerBCs=table({'Enr1','Enr2','Enr3','Enr4'}',{'TCTACTCT','CTCCTTAC','TATGCAGT','TACTCCTT'}');
bamFiles=struct2table(bamFiles);
[~,enrIdx]=ismember(regexp(bamFiles.name,'Enr\d','match','once'),outerBCs.Var1);
[~,innerIdx]=ismember(regexp(bamFiles.name,'(?<=_)\d+[A-H](?=)','match','once'),innerBCs.Var5);
bamFiles.barcode=strcat(outerBCs.Var2(enrIdx),'_',innerBCs.Var6(innerIdx));
sampleTable=readtable('/home/labs/barkailab/LAB/data/SEQ/Gilad_X34Spike_210218/X34_DNA_index.xlsx','ReadVariableNames',false);
[~,sampleIdx]=ismember(bamFiles.barcode,sampleTable.Var1);
bamFiles.sample(sampleIdx>0)=sampleTable.Var3(sampleIdx(sampleIdx>0));
badSample=find(cellfun('isempty',bamFiles.sample))
bamFiles(badSample,:)=[];
bamFiles.par=str2double(regexp(bamFiles.sample,'(?<=spike)\d+','match','once'));
bamFiles.par(isnan(bamFiles.par))=0;
bamFiles.ratio=bamFiles.par./(100-bamFiles.par);

readsChr=cell(1,size(bamFiles,1));
parfor i=1:size(bamFiles,1)
    [~,out]=system(sprintf(samCommand,[bamFiles.folder{i} '/' bamFiles.name{i}]));
    out=regexp(out,'(\d+) (\S+)','tokens')';
    out=cat(1,out{:});
    [~,idx]=ismember(out(:,2),chrNames);
    readsChr{i}=zeros(numel(chrNames),1)
    readsChr{i}(idx)=str2double(out(:,1));
end
clearvars -except bamFiles readsChr
readsChr=cat(2,readsChr{:})
bamFiles.forFit(:)=1
[tags,~,bamFiles.tagId]=unique(regexp(bamFiles.sample,'^[A-Za-z]+-.*?(?=-)','match','once'));
intRatio={[1:16],[18:34];[1:11,13:16],[18:28,30:33]}%[1],[18];2,19;3,20;4,21;5,22;6,23;7,24;8,25;9,26;10,27;11,28;12,[29 34];13,30;14,31;15,32;16,33;}
cmpTag=[1,3;2,4]
mixRatio=unique(bamFiles.ratio);
close all
figure
for r=2
    subplot(1,1,1)
    seqRatio=sum(readsChr(intRatio{r,2},:),1)./sum(readsChr(intRatio{r,1},:),1);
    c=0;
    for t=[2,4]
        c=c+1;
        subplot(2,numel(tags),c)
        xVal=-acol(log2(bamFiles.ratio(bamFiles.tagId==t)))
        yVal=-acol(log2(seqRatio(bamFiles.tagId==t)))
        scatter(xVal,yVal,'o','filled','DisplayName',tags{t})
        hold on
        p=fliplr(linortfit2(xVal(bamFiles.forFit(bamFiles.tagId==t)==1),yVal(bamFiles.forFit(bamFiles.tagId==t)==1)));
        plot(sort(xVal),sort(xVal)*p(2)+p(1),'-','DisplayName',sprintf('%.2fx%+.2f',p(2),p(1)))
    end
    legend('show','Location','best')
    xlabel('log2(Pellet ratio (Scer:Spar))')
    ylabel('log2(Seq ratio (Spar:Scer))')
    title(' myc quantification with sike-ins')
end
save_gf(gcf,'spike-in')
subplot(2,1,2)
hold off
for m=2:numel(mixRatio)
    mycLevel(m,:)=[seqRatio(bamFiles.tagId==cmpTag(1,1)&bamFiles.ratio==mixRatio(m))./seqRatio(bamFiles.tagId==cmpTag(1,2)&bamFiles.ratio==mixRatio(m)),...
        seqRatio(bamFiles.tagId==cmpTag(2,1)&bamFiles.ratio==mixRatio(m))./seqRatio(bamFiles.tagId==cmpTag(2,2)&bamFiles.ratio==mixRatio(m))]
end
for t=1:size(mycLevel,2)
    plot(mixRatio,mycLevel(:,t),'DisplayName',sprintf('%s vs %s',tags{cmpTag(t,1)},tags{cmpTag(t,2)}),'LineWidth',2)
    hold on
end
xlabel('Pellet Ratio (Spar:Scer)')
ylabel('absolute myc level')
title('absolute myc level')
%%
clear all
GP=load('group_imp.mat');
load('forFigure15.mat')
hisNucs=find(ismember(nucMeta.order,intNuc{n}) & ismember(nucMeta.gene,GP.groups{23}{2}{45}));

%% HA comparissons
clear all
load('forFigure15.mat')
intSmp=find(contains(expTypes.ab,{'h3','ha'}) & contains(expTypes.con,'async')& contains(expTypes.gt,'wt')& contains(expTypes.tag2,{'hht','hhf','h3','h4'}));
h3To=diff(log(medNuc(:,[37 37+54])+1),1,2);
sampleOrder=[1,3,2,4,6,12,11,10,9,8,5,7];
intSmp=intSmp(sampleOrder)
imagesc(corr(medNuc(:,intSmp),'rows','pairwise'))
yticks(1:12)
yticklabels(strcat(expTypes.ab(intSmp),';',expTypes.tag2(intSmp)))


xticks(1:12)
xticklabels(strcat(expTypes.ab(intSmp),';',expTypes.tag2(intSmp)))
xtickangle(90)
title('genome-wide correlaiton occupancy')
xId=2
specialNuc=ismember(nucMeta.order,[1]);
c=0;
for s=setdiff(intSmp',xId)
    c=c+1;
    subplot(3,4,c)
    nonNan=all(~isnan(medNuc(:,[s,xId])),2);
    %scatter(medNuc(~specialNuc,xId),medNuc(~specialNuc,s),'.')
    hold off
    dscatter(medNuc(~specialNuc&nonNan,xId),medNuc(~specialNuc&nonNan,s));
    hold on
    plot(xlim,xlim,'k--')
    xlabel(strjoin(expTypes{xId,[1,4]}))
    ylabel(strjoin(expTypes{s,[1,4]}))
end

c=0;
xId=2
for s=setdiff(intSmp',xId)
    c=c+1
    subplot(3,4,c)
    hold off
    nonNan=all(~isnan(medNuc(:,[xId,s])),2);
    dscatter(h3To(nonNan),diff(log(medNuc(nonNan,[xId s])+1),1,2))   
    title(sprintf('%s: %.2f',strjoin(expTypes{s,[1,4]}),corr(h3To(nonNan),diff(log(medNuc(nonNan,[xId s])+1),1,2))))
    hold on
    plot(xlim,[0 0],'k-')
    caxis([-1 1.5])
    xlabel('H3 clean turnover')
    ylabel(sprintf('\\Delta Log (Ha vs. Ha_{%s})',strjoin(expTypes{xId,[1 4]})))
    axis tight
end

selSmp=find(contains(smeta.ab,{'ha','h3'})&contains(smeta.con,'async')&contains(smeta.gt,'wt')&contains(smeta.tag2,{'h3','hht2','h4','hhf1'})&~contains(smeta.tag2,{'tev'},'IgnoreCase',true)&smeta.bad==0)
figure
imagesc(corr(nucData(:,selSmp),'rows','pairwise'))
caxis([.5 1])
yticks(1:numel(selSmp))
yticklabels(strcat(smeta.ab(selSmp),'-',smeta.tag2(selSmp),'-',smeta.exp(selSmp)))
%% histone  marks
clear all
%load('forFigure4add.mat')
load('forFigure16.mat')
load('epMarks.mat')
% turnover
[~,~,expTypes.sid]=unique(expTypes(:,[2:4]),'stable');
nExp=size(expTypes,1)
c=0;
for i=unique(expTypes.sid(strcmp(expTypes.ab,'ha')),'stable')'
    c=c+1;
    hai=find(strcmp(expTypes.ab,'ha')&expTypes.sid==i);
    myci=find(strcmp(expTypes.ab,'myc')&expTypes.sid==i);
    medNuc(:,nExp+c)=log(medNuc(:,myci)+1)-log(medNuc(:,hai)+1);
    expTypes(nExp+c,:)=expTypes(hai,:);
    expTypes.ab(nExp+c,:)={'logTO'};
end
%enrmark
enrMark=nan(size(nucMark));
for i=find(~contains(markNames.Var1,'input')')
    inputId=find(contains(markNames.Var1,markNames.Var2(i)),1);
    enrMark(:,i)=log(nucMark(:,i)+1)-log(nucMark(:,inputId)+1);
end
gbNucs=nucMeta.order>1%true(size(nucData,1),1);
promNucs=ismember(nucMeta.order,[-2 -1])
for i=find(~contains(markNames.Var1,'input')')
    crAll(i)=corr(medNuc(:,39),enrMark(:,i),'rows','pairwise');
    crGB(i)=corr(medNuc(gbNucs,39),enrMark(gbNucs,i),'rows','pairwise');
    crProm(i)=corr(medNuc(promNucs,39),enrMark(promNucs,i),'rows','pairwise');
end
intCr=contains(markNames.Var1,{'K79','K56','k9'},'IgnoreCase',true)
figure
subplot(1,2,1)
scatter(crAll,crGB,'filled')
text(crAll(intCr),crGB(intCr),markNames.Var1(intCr))
xlabel('correlation turnover mark enrichment - over all Nucleosomes')
ylabel('correlation turnover mark enrichment - over gene body nucleosomes')
title('all vs. gene body')
subplot(1,2,2)
scatter(crAll,crProm,'filled')
text(crAll(intCr),crProm(intCr),markNames.Var1(intCr))
xlabel('correlation turnover mark enrichment - over all Nucleosomes')
ylabel('correlation turnover mark enrichment - over promoter nucleosomes')
title('all vs. promoter')

%% histone expression during stress
figure;
for i=2:3
    plot(ExpData(i).time,median(dynExp{i}(GP.groups{23}{2}{45},:)),'DisplayName',ExpData(i).Name,'LineWidth',2)
    hold on
end
legend()
xlabel('time after release')
ylabel('median histone expression')

%% median myc/HA
clear all
load('forFigure15.mat')
%%comapre Ha / myc samples [39,52,50,93,106,104]
sampleOrder=[56,55,53,54,58,37,34,48,52,39,50,47,15,16,5,109,107,108,112,91,88,102,106,93,104,101,69,70,59,110]
sampleOrder=[[34,52,39],[34,52,39]+54]
figure;
imagesc(corr(medNuc(:,sampleOrder),'rows','pairwise'))
yticks(1:numel(sampleOrder))
yticklabels(strcat(num2str(sampleOrder'),'-',expTypes.ab(sampleOrder),'-',strrep(expTypes.tag2(sampleOrder),'_',''),'-',num2str(expTypes.n(sampleOrder))))
figName='SF1_H3X2HHT2'
save_gf(gcf,figName)
%% median myc/HA
clear all
load('forFigure17.mat')
clearvars 
expTypes.tp=str2double(regexp(expTypes.con,'(?<=t)\d+','match','once'))
[expTypes,idx]=sortrows(expTypes,{'ab','gt','tp'});
medNuc=medNuc(:,idx);
[tcs,~,expTypes.tcId]=unique(expTypes(:,[1,2]))
hisNucs=[2733,2734,12699];
for t=[3,4]
    hisDyn=median(log(medNuc(hisNucs,expTypes.tcId==t))-log(median(medNuc(hisNucs,expTypes.tcId==t),2)));
    hisDyn=hisDyn-hisDyn(1);
    plot(expTypes.tp(expTypes.tcId==t),movmean(hisDyn,2),'-','LineWidth',2,'DisplayName',strjoin(tcs{t,:}))    
    hold on
    xlabel('time at 37C')
    ylabel('median(rel. change histone myc) (log)')
end
legend()

scatter(expTypes.tp,median(log(medNuc(hisNucs,:))))

%% compare all