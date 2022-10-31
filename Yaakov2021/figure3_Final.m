function figure2_Final
[expTypes,medNuc,smeta,nucData,nucMeta]=standardImport(2)
load('segal_occ.mat','p_all')
nucSegal=p_all(nucMeta.pos);
clear p_all regev holstege ToR val idx 
save('forFigure2add.mat')

%% sequence dependency
clear all
load('forFigure2add.mat')
load('SC_genome.mat')
allSeq=extractfield(SC_genome,'Sequence');
allSeq=cat(2,allSeq{:});
ATcon=ismember(allSeq',{'A','T','a','t'});
ATper=movmean(ATcon,100);
nucAT=ATper(nucMeta.pos);
clearvars allSeq SC_genome ATcon ATper

figure;
subplot(1,2,1)
hold off
intVal=[.5 .75];
intNames={'AT low','AT high'}
intHis={'h3','h2a'}
histBins=[0:0.3:30];
for i=1:numel(intVal)
    selNucs=nucAT>(intVal(i)-.02)&nucAT<(intVal(i)+.02);
    for h=1:numel(intHis)
        selExp=strcmp(expTypes.tag2,intHis(h)) & strcmp(expTypes.ab,'ha');
        nucHist=histcounts(medNuc(selNucs,selExp),histBins,'Normalization','pdf')
        lineName=sprintf('%s nucs %s HA lvl',intNames{i},intHis{h});
        plot(movmean(histBins,2,'Endpoints','discard'),movmean(nucHist,5),'-','DisplayName',lineName,'LineWidth',2)
        hold on
    end
    legend()
end
%xticks((numel(intHis)+1)*(.5+[1:numel(intVal)]))
%xticklabels(intVal)
xlabel('HA signal')
ylabel('PDF')

%% for figure 2B 
intPar={'nucSegal'}%{'nucMeta.ToR'}%;
intNuc={min(nucMeta.order):max(nucMeta.order)}%{setdiff(min(nucMeta.order):max(nucMeta.order),0)};
genePar=0;
blackList=ismember(nucMeta.gene,GP.groups{23}{2}{45});
seperateCr=0
intExp=[5,16,3,14]
movPar=.1
stepSize=40
hisNucs=ismember(nucMeta.gene,GP.groups{23}{2}{45}) & nucMeta.order<0;

intPar={'geneLvl'};
genePar=1;
blackList=ismember(nucMeta.gene,GP.groups{23}{2}{45});
seperateCr=0
movPar=.25
stepSize=5
%intNuc={3,[2],1,[-1],-2};
intNuc={-2,-1,1,2,3}
% intExp=[1,10]
% selPlot=[2,3];
% selPlot=[6,7];

%intPar={'geneStdH'}
%genePar=1;
%blackList=ismember(nucMeta.gene,GP.groups{23}{2}{45});
%seperateCr=1
%movPar=.05
%stepSize=5

allExp={[3,14] [4 15],[5 16] [6 17],[11 22]}%{[3,14] [4 15],[5 16] [6 17],[7 18]}
%allExp={[5 16,7 18]}
%intExp=[5 16]
%intExp=[1,12]
%selPlot=[2];
%intExp=[2,13]
%selPlot=[3];
allExp={[6,22],[1,17],[4,20],[2,18],[7,23]}


lineColor=lines(max(cellfun('prodofsize',allExp)));
lineColor([1,3],:)=[.4 .4 .4;.6 .6 .6]
close all
for a=1:numel(allExp)
    intExp=allExp{a};
    for p=1:numel(intPar)
        parData=eval(intPar{p});
        for n=1:numel(intNuc)
            %subplot(2,4,selPlot(n)) %for wt
            %subplot(1,5,selPlot(p)) % for spt6oe
            subplot(numel(allExp),numel(intNuc),(a-1)*numel(intNuc)+n)
            hold off
            if genePar(p)==1
                selNucs=find(ismember(nucMeta.order,intNuc{n}) & ismember(nucMeta.gene,find(~isnan(parData))) &~blackList & all(~isnan(medNuc(:,intExp)),2) );
                [xVal,idx]=sort(parData(nucMeta.gene(selNucs)));
            else
                selNucs=find(ismember(nucMeta.order,intNuc{n}) & ~isnan(parData) &~blackList & all(~isnan(medNuc(:,intExp)),2) );
                [xVal,idx]=sort(parData(selNucs));
            end
            selNucs=selNucs(idx);
            clear lineObj
            for e=1:numel(intExp)
                if contains(expTypes.tag2(intExp(e)),{'2'})
                    lineColor=[100,100,100;127,63,152];
                    yLimit=[2 20];
                else
                    lineColor=[100,100,100;0,174,239]
                    yLimit=[2 15];
                end
                expVal=medNuc(selNucs,intExp(e));
                yVal=movmean(expVal,movPar(p),'omitnan','SamplePoints',makeUnique(xVal,4))
                yValE=movstd(expVal,movPar(p),'omitnan','SamplePoints',makeUnique(xVal,4))./sqrt(movsum(~isnan(expVal),movPar(p),'omitnan','SamplePoints',makeUnique(xVal,4)));
                if seperateCr(p)==1
                    cr1=corr(xVal(xVal<=0),expVal(xVal<=0),'rows','pairwise');
                    cr2=corr(xVal(xVal>=0),expVal(xVal>=0),'rows','pairwise');
                    lineName=sprintf('%s,%.2f,%.2f',strjoin(expTypes{intExp(e),[2,1,4]}),cr1,cr2);
                else
                    cr=corr(xVal,expVal,'rows','pairwise');
                    lineName=sprintf('%s,%.2f',strjoin(expTypes{intExp(e),[2,1,4]}),cr);
                end
                lineObj(e)=plot(xVal(1:stepSize:end),yVal(1:stepSize:end),'-','Linewidth',2,'Color',brighten(lineColor(e,:),0),'DisplayName',lineName)
                hold on
                fill([xVal(1:stepSize:end);flipud(xVal(1:stepSize:end))],[yVal(1:stepSize:end)-yValE(1:stepSize:end);flipud(yVal(1:stepSize:end)+yValE(1:stepSize:end))],brighten(lineColor(e,:),0),'LineStyle','none','FaceAlpha',.3)
                %scatter(parData(hisNucs),medNuc(hisNucs,intExp(e)),[],lineColor(e,:),'o','filled')
            end
            xlabel(intPar{p})
            if contains(intPar{p},'geneLvl')
                xlim([5 13])
            end
            %ylim(yLimit)
            ylabel(sprintf('moving average %.2f',movPar(p)))
            legend(lineObj)
            title(sprintf('%+d - %+d nuc',intNuc{n}(1),intNuc{n}(end)))
        end
    end
end
%save_gf(gcf,sprintf('Fig2_%s','Spt6-h3'))
%save_gf(gcf,'Fig2_nucScore')
save_gf(gcf,'SF2_allLines')
%% scatterPlots
selNucs=nucMeta.order==2 & ismember(nucMeta.gene,find(~isnan(geneLvl)));
intExp=[3,14;16,14]
selPlot=[1,5]
for e=1:size(intExp,1)
    subplot(2,4,selPlot(e))
    hold off
    scatter(medNuc(selNucs,intExp(e,1)),medNuc(selNucs,intExp(e,2)),[],geneLvl(nucMeta.gene(selNucs)),'.')
    xlabel(strjoin(expTypes{intExp(e,1),:}))
    ylabel(strjoin(expTypes{intExp(e,2),:}))
    title('+2 nucleosome')
    ylabel(colorbar('east'),'genelvl')
    axis tight
end

%% correlation matrices
intNuc=[-3:-1 1:3]
%intExp=[1,10;3,12]; wt adjust
%selPlot=[4,8] % wt adjust

intExp=[12,1;13,2];
selPlot=[4,5]

for e=1:size(intExp,1)
    for n=1:numel(intNuc)
        selNucs=nucMeta.order==intNuc(n) & ~blackList;
        crMat(n,:)=corr(geneLvl(nucMeta.gene(selNucs)),medNuc(selNucs,intExp(e,:)),'rows','pairwise')
    end
    %subplot(2,4,selPlot(e))
    subplot(1,size(intExp,1),e)
    imagesc(crMat')
    yticks(1:size(intExp,2))
    set(gca,'XAxisLocation','top')
    yticklabels(strcat(expTypes.tag2(intExp(e,:)),'-',expTypes.ab(intExp(e,:)),'-',expTypes.gt(intExp(e,:))))
    xticks(1:numel(intNuc));
    xticklabels(intNuc)
    colorbar()
end
save_gf(gcf,'SF2_spt6Corr')
%% scatterPlots - spt6OE

selNucs=nucMeta.order==2 & ismember(nucMeta.gene,find(~isnan(geneLvl)));
intExp=[13,12];
conExp=[16,14];
subplot(1,5,1)
hold off
xVal=log(medNuc(selNucs,intExp(1))+1)-log(medNuc(selNucs,conExp(1))+1);
yVal=log(medNuc(selNucs,intExp(2))+1)-log(medNuc(selNucs,conExp(2))+1);

scatter(xVal,yVal,[],geneLvl(nucMeta.gene(selNucs)),'.')
xlabel([strjoin(expTypes{intExp(1),[1,2,4]}) '-' strjoin(expTypes{conExp(1),[1,2,4]})] )
ylabel([strjoin(expTypes{intExp(2),[1,2,4]}) '-' strjoin(expTypes{conExp(2),[1,2,4]})] )
title('+2 nucleosome')
ylabel(colorbar('east'),'genelvl')
axis tight
save_gf(gcf,'Fig2_spt6OE')
end