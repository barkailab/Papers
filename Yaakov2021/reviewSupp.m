%% calculate Histone expression 
fileLocs=readtable('~/gilad/RNA1_A09_TGGAGCA/FileLocs.xlsx','ReadVariableNames',false);
for i=1:size(fileLocs,1)
    temp=dir([fileLocs.Var3{i},'/**/*' extractAfter(fileLocs.Var1{i},'_') '*R1*.*gz'])
    fileLocs.path{i}=[temp(1).folder '/' temp(1).name]
end
[samples,~,fileLocs.sid]=unique(fileLocs.Var2)
for i=1:numel(samples)
    infiles=sprintf('%s ',fileLocs.path{fileLocs.sid==i})
    cmd=sprintf('cat %s > %s',infiles,['~/gilad/RNA1_A09_TGGAGCA/',samples{i},'.fastq.gz'])
    system(cmd)
end
clear all
outFiles=struct2table(dir('~/gilad/RNA1_A09_TGGAGCA/*.out'))
for i=1:size(outFiles)
    temp=readtable([outFiles.folder{i} '/' outFiles.name{i}],'FileType', 'text','ReadVariableNames',false);
    temp.Var2(temp.Var2==6702)=GP.gene_table.HHT1;
    temp.Var2(temp.Var2==6703)=GP.gene_table.HTA2;
    geneExp(:,i)=accumarray(temp.Var2,temp.Var4,[6701 1],@sum,nan);
    temp(temp.Var2==GP.gene_table.HHT1,:)
end
geneExp(all(geneExp==0,2),:)=nan;
xId=4;
figure
for i=1:3
    subplot(1,3,i)
    hold off
    scatter(log(geneExp(:,xId)+.5),log(geneExp(:,i)+.5),'.')
    hold on
    scatter(log(geneExp(GP.groups{23}{2}{45},xId)+.5),log(geneExp(GP.groups{23}{2}{45},i)+.5),'o','filled')
    xlim([3 10])
    text(log(geneExp(GP.groups{23}{2}{45},xId)+.5),log(geneExp(GP.groups{23}{2}{45},i)+.5),GP.gene_infoR64.name(GP.groups{23}{2}{45}))
    xlabel(strrep(outFiles.name{xId},'_',' '))
    ylabel(outFiles.name{i})
    title('log #reads+.5')
end

%% calculat ehiostoen expression cyclebase
expData=readtable('~/budding_timecourses/budding_experiments.tsv','FileType','text','Delimiter','\t');
metaData=readtable('~/budding_timecourses/budding_metadata.tsv','FileType','text','Delimiter','\t');

for i=1:size(metaData,1)
    fieldName=strrep(strrep(metaData.source{i},'-','_'),' ','');
    dataLines=strcmp(expData.source,metaData.source{i});
    a=strrep(strrep(expData.y_values_expression_(dataLines),'}',''),'{','');
    a=regexp(a,',','split');
    a = cellfun(@(x)str2double(x),a,'UniformOutput',false);
    a=cat(1,a{:});
    [~,barkaiId]=ismember(expData.identifier(dataLines),GP.gene_infoR64.orf);
    expMat=nan(6701,size(a,2));
    expMat(barkaiId(barkaiId>0),:)=a(barkaiId>0,:);
    tps=str2double(regexp(metaData.x_values__cell_cycle_{i},'\d+.\d+','match'));
    cyclebase.(fieldName).tp=tps;
    cyclebase.(fieldName).exp=expMat;
end
clear all

cyclebase=load('cyclebaseRaw.mat')
GP=load('/home/labs/barkailab/felixj/Documents/MATLAB/scripts/gene/group_imp.mat');
dataSets=fieldnames(cyclebase)';
figure
c=0;
for intGene=GP.groups{23}{2}{45}
    c=c+1;
    subplot(2,4,c)
    clear temp
    for ds=1:numel(dataSets)
        plot(cyclebase.(dataSets{ds}).tp,cyclebase.(dataSets{ds}).exp(intGene,:),'DisplayName',dataSets{ds})
        hold on
        temp(ds,:)=interp1(cyclebase.(dataSets{ds}).tp,cyclebase.(dataSets{ds}).exp(intGene,:),100:300);
    end
    plot(100:300,median(temp),'k-','LineWidth',2)
    title(GP.gene_infoR64.nameNew(intGene))
    xlim([100 300])
    ylim([-6 6])
    ylabel('Exprssion Level-log2?')
    xlabel('CycleBaseTime')
end
%% calculate histone expresion Yoav
clear all
temp=combine_rna_data;
GP=load('group_imp.mat');
intRnaData=[1,3,40]
close all
figure
phaseColor=lines(2);
c=0;
for i=intRnaData
    if c==5
        c=c+5
    end
    c=c+1;
    subplot(3,numel(intRnaData),c)
    hold off
    expMat=temp(i).reads;
    expMat(sum(mean(expMat==0,2)>0.5),:)=NaN;    
    normExp=expMat./sum(expMat,1,'omitnan')*100000;
    logExp=log2(normExp+0.1);
    deltaExp=logExp-median(logExp,2,'omitnan');
    
    imagesc(deltaExp(GP.groups{23}{2}{45},:),[-2.5 2.5])
    if c==1
        yticks(1:numel(GP.groups{23}{2}{45}))
        yticklabels(GP.gene_infoR64.name(GP.groups{23}{2}{45}))
    end
    title(strrep(temp(i).name,'_',' '))
    xlabel('timepoints')
    caxis([-2.5 2.5])
    if c==5
        colorbar()
    end
    subplot(3,numel(intRnaData),c+numel(intRnaData))
    hold off
    scatter(temp(i).tps,median(deltaExp(GP.groups{23}{2}{45},:)),[],[.6 .6 .6],'filled','LineWidth',2,'DisplayName','Histones (median)')
    hold on
    xAll=[1:max(temp(i).tps)];
    yAll=interp1(temp(i).tps(temp(i).tps>=0),median(deltaExp(GP.groups{23}{2}{45},temp(i).tps>=0)),[1:max(temp(i).tps)]);
    plot(xAll,yAll,'-','LineWidth',2,'Color',[.6 .6 .6])   
    xlim([0 min([max(xlim),100])])
    hold on
    [~,sphase]=maxk(yAll(1:60),15);
    sphase=sort(sphase);
    [~,notSphase]=mink(yAll(1:60),15);
    notSphase=sort(notSphase);
    scatter(sphase,yAll(sphase),20,phaseColor(1,:),'filled','DisplayName','SPhase')
    scatter(notSphase,yAll(notSphase),20,phaseColor(2,:),'filled','DisplayName','G1/G2/M')
    xlabel('time after release')
    %plot(temp(i).tps,median(deltaExp(GP.groups{8}{2}{9},:),'omitnan'),'LineWidth',2,'DisplayName','mating')
    ylim([-2.5 2.5])
    if c==1
        legend()
    end
    xlabel('time after release')
    ylabel('Median Histone expression dynamics (log2)')
    subplot(3,numel(intRnaData),c+2*numel(intRnaData))
    hold off
    yAbs=2.^(yAll-median(yAll(sphase)))*3000
    plot(xAll,yAbs,'-','Color',[.6 .6 .6],'LineWidth',2)   
    hold on
    scatter(sphase,yAbs(sphase),20,phaseColor(1,:),'filled')
    scatter(notSphase,yAbs(notSphase),20,phaseColor(2,:),'filled')
    xlabel('time after release')
    ylabel('estimated histone production')
    title(sprintf('%d histones/min outside Sphase',round(mean(yAbs(notSphase)))))
    xlim([0 80])
end

set(gcf,'Color',[1 1 1])
figureName='RNA datasets';
save_gf(gcf,figureName)

    %for geneId=GP.groups{23}{2}{45}
    %    plot(temp(i).tps,normExp(geneId,:))
    %    hold on
    %end
    %plot(temp(i).tps,quantile(normExp,.8),'k--')
% for i=[40]
%     figure
%     temp(i).norm=100000.*temp(i).reads./sum(temp(i).reads,'omitnan')+0.1
%     acMat=nan(numel(GP.groups{23}{2}{45}),max(temp(i).tps))
%     c=0
%     for intGene=GP.groups{23}{2}{45}
%         c=c+1;
%         subplot(2,8,c)
%         hold off
%         plot(temp(i).tps,temp(i).norm(intGene,:),'-','LineWidth',2)
%         acMat(c,:)=interp1(temp(i).tps,temp(i).norm(intGene,:),1:max(temp(i).tps));
%         xlabel('time after release')
%         ylabel('norm MRNA - no log')
%     end
%     [x,y]=findpeaks(movmean(xcorr(median(log(acMat),1),'unbiased'),5),[-139:139]);
%     cmbExp=cat(3,acMat(:,1:59),acMat(:,60:118),[acMat(:,119:end) nan(8,37)]);
%     c=8
%     for g=1:size(cmbExp,1)
%         c=c+1;
%         subplot(2,8,c)
%         hold off
%         for cc=1:size(cmbExp,3)
%             plot(cmbExp(g,:,cc))
%             hold on;
%         end
%         plot(median(cmbExp(g,:,:),3,'omitnan'),'k-','LineWidth',2)        
%         title(sprintf('%s %.1f %.1f',GP.gene_infoR64.name{GP.groups{23}{2}{45}(g)},mean(maxk(median(cmbExp(g,:,:),3,'omitnan'),20)),mean(mink(median(cmbExp(g,:,:),3,'omitnan'),39))))
%         xlabel('aligned cell cycles')
%         ylabel('norm mRNA')
%     end
% end

%% esitmate turnover numbers
clear all
%[expTypes,medNuc,smeta,nucData,nucMeta]=standardImport(19,'save',true);
load('../Gilad-CC/forFigure20.mat','smeta','nucData','nucMeta')
%load('forFigure3add.mat');expTypes=smeta(:,[7,9,10,14]);medNuc=nucData;
load('forFigure15.mat')
GP=load('group_imp.mat');
ToR = readtable('/home/labs/barkailab/yuliago/Documents/MATLAB/ToR_raz.csv');
ToR.idx=GP.chrIdx(ToR.chr)+ToR.start+500;
[val,idx]=min(abs(nucMeta.pos-ToR.idx'),[],2);
nucMeta.ToR(val<501)=ToR.score(idx(val<501));
nucMeta.ToR(nucMeta.ToR==0)=NaN;
regev=load ('/home/labs/barkailab/LAB/data/DataExternal/EXPRESSION/regev_pnas_2008_rnaseq/data');
geneLvl=median(regev.regev_rnaseq,2,'omitnan')
clear ToR regev
fitType={'h2a,h2b','h3,h4,hht2,hhf1'}
fitNucsColl{1}= nucMeta.gene==0  | ismember(nucMeta.gene,find(geneLvl<8));
fitNucsColl{2}=nucMeta.ToR>25 & nucMeta.order>-1;%true(size(nucMeta.ToR))%
load('qPCRtargets.mat','qPCRtargets')
intNucs=[qPCRtargets.nuc(ismember(qPCRtargets.gene,{'HTA1','KOG1','AUS1','HHF1','HTB2'}))]
%smpID=[find(contains(expTypes.ab,'ha','IgnoreCase',true)&contains(expTypes.gt,{'wt'},'IgnoreCase',true)&contains(expTypes.con,{'async','h2o2T38','37t50','alpha2hr','nacl33'},'IgnoreCase',true)&contains(expTypes.tag2,{'h3'},'IgnoreCase',true)&expTypes.n>0)]
%smpID=smpID([8,7,2,4,6])
smpID=[find(contains(expTypes.ab,'ha','IgnoreCase',true)&contains(expTypes.gt,{'wt','asf','hir'},'IgnoreCase',true)&contains(expTypes.con,{'async'},'IgnoreCase',true)&contains(expTypes.tag2,{'h3','hht2','h4','hhf'},'IgnoreCase',true)&~contains(expTypes.tag2,{'tev','slow','nc'},'IgnoreCase',true)&expTypes.n>0)]
%smpID=find(contains(expTypes.ab,'ha','IgnoreCase',true) & contains(expTypes.con,'async')& contains(expTypes.gt,'wt')& ~contains(expTypes.tag2,{'nc','none','slow'}))%[29,25,17,15]%[31,11,46,33];
%smpID=find(strcmp(expTypes.ab,'ha'))
smpID=smpID([4,1,2,3,5,11,9,10])
ccTime(smpID)=[102,158,120 , 120 ,102,102, 140 ,143]
load('greytoblue.mat','Colo1')
cMapModel=flipud(bone)
%smpID=[30]
%smpId=[34,58,55]'%[39,52]'
subPlot=[1:ceil(numel(smpID)/2),ceil(numel(smpID)/2)*2+1:ceil(numel(smpID)/2)*3]
ccTime=repmat(100,size(expTypes,1),1)
close all
figure
c=0
fittingChoice='pMid';
haBins=[7.5:1:16.5]%[6.5:1:14.5]%[0.7:0.1:1.6];
expTypes.midMyc(:)=NaN;
for xId=smpID
    yId=find(ismember(expTypes(:,2:4),expTypes(xId,2:4)) & contains(expTypes.ab,'myc'));
    %xId=29
    fitNucs=fitNucsColl{2};%fitNucsColl{contains(fitType,regexp(expTypes.tag2(xId),'h\d','match','once'))};
    
    xVal=medNuc(:,xId);
    %yVal=medNuc(:,yId);
    
    [pRobust,S]=robustfit(xVal(fitNucs),yVal(fitNucs));
    
    medNucs=abs(1-xVal./median(xVal,'omitnan'))<=0.05;
    medX=median(xVal(medNucs),'omitnan');
    medY=median(yVal(medNucs));
    pSimple=[0 medY./medX];
    %nucBin=sum(nucData(:,xId)>(haBins.*medX),2);
    nucBin=sum(xVal>(haBins),2);
    lineNucs=nucBin>0 & nucBin<numel(haBins);
    medLineX=accumarray(nucBin(lineNucs),xVal(lineNucs),[],@median);
    medLineY=accumarray(nucBin(lineNucs),yVal(lineNucs),[],@(x)quantile(x,0.475));
    pMid=fliplr(polyfit(medLineX,medLineY,1));
    pBlack=[0 median(yVal,'omitnan')./median(xVal,'omitnan')];
    pMix=[pMid(1)+median(yVal)- pMid(2)];       
    
    pOrig=eval(fittingChoice);%fliplr(pMid)%;
    pExp(xId,:)=pOrig;
    
    nucFac=median(xVal,'omitnan')./1;
    adjX=xVal/nucFac;
    adjY=(yVal-pOrig(1))./(pOrig(2)*nucFac);
    pAdj=robustfit(adjX(fitNucs),adjY(fitNucs));
    
    adjLineX=medLineX./nucFac;
    adjLineY=(medLineY-pOrig(1))./(pOrig(2)*nucFac);
    pAdj=fliplr(polyfit(adjLineX,adjLineY,1));

        
    mycAdd2=adjY-adjX;
    adjColl{xId}=[adjX,adjY];
    
    
    % figure
    c=c+1;
    subplot(4,ceil(numel(smpID)/2),subPlot(c))
    hold off
    nonNan=all([xVal yVal]>0,2);
    dscatter(xVal(nonNan),yVal(nonNan));    
    hold on    
    scatter(medLineX,medLineY,'ro','filled')
    scatter(xVal(intNucs),yVal(intNucs),'ro','filled','DisplayName','qPCR nucleosomes')
    baseLine=plot(xlim(),xlim()*pOrig(2)+pOrig(1),'k-','LineWidth',2,'DisplayName',sprintf('BaseLine: %.2fx%+.2f',pOrig(2),pOrig(1)));
    %plot(xlim(),xlim()*pSimple(2)+pSimple(1),'r-','LineWidth',2,'DisplayName',sprintf('%.2fx%+.2f',pSimple(2),pSimple(1)))
    xlabel([strjoin(expTypes{xId,1:4}),' (relative)'])
    ylabel([strjoin(expTypes{yId,1:4}),' (relative)'])
    %title(expTypes.con{yId})
    title(sprintf('%s %s ChIPSeq cr:%.2f',expTypes.tag2{yId},expTypes.gt{yId},corr(xVal,yVal,'rows','pairwise')))
    set(legend(baseLine),'Location','best')
    if mod(c,ceil(numel(smpID)/2))==0
        ylabel(colorbar(),'Nucleosome density')
    end
    colormap(gca,Colo1(20:end,:))
    text(xVal(intNucs),yVal(intNucs),GP.gene_infoR64.nameNew(nucMeta.gene(intNucs)));
    subplot(4,ceil(numel(smpID)/2),subPlot(c)+ceil(numel(smpID)/2))
    hold off
    scatter(adjX,adjY,[],mycAdd2,'.');
    colormap(gca,cMapModel(50:end,:))
    hold on
    scatter(adjX(intNucs),adjY(intNucs),[],'ro','filled');
    text(adjX(intNucs),adjY(intNucs),num2str(mycAdd2(intNucs),2));
    baseLine=plot(xlim(),xlim()*pAdj(2)+pAdj(1),'k-','LineWidth',2,'DisplayName','Replication dependent integration')
    set(legend(baseLine),'Location','best')
    %xlabel('"abs" nucleosome Occ')
    %ylabel('Myc integration events')
    ylabel('Myc/cell cycle')
    xlabel('HA/median(HA)')

    caxis([-1 5])
    if mod(c,ceil(numel(smpID)/2))==0
        ylabel(colorbar(),'non S-phase myc cintegration events')
    end
    title(sprintf('%.1fk/cell cycle outside SPhase (%.0f/min)',round(sum(mycAdd2,'omitnan')/1000,1),round(sum(mycAdd2,'omitnan')/ccTime(xId))))
end
set(gcf,'Color',[1 1 1])
figureName='H3oldforMain';
%save_gf(gcf,figureName)
%%
%% turnover numbers for individual samples
smeta.mycAdd(:)=NaN;
clear mycAdd pKeep
figure
c=0;
%smeta.expId(smeta.expId==37,:)=34
%smeta.expId(smeta.expId==37+54,:)=34+54

for xId=find(strcmp(smeta.ab,'ha')&ismember(smeta.expId,smpID))'
    yId=find(strcmp(smeta.name,strrep(smeta.name{xId},'ha','myc')));
    fitNucs=fitNucsColl{2};
    if numel(yId)>0 
        if contains(smeta.tag2(xId),'h3') & contains(smeta.gt(xId),'wt')
            xVal=medNuc(:,34);%
        else
            xVal=medNuc(:,smeta.expId(xId));%nucData(:,xId);%
        end
        yVal=nucData(:,yId);%medNuc(:,smeta.expId(yId));%
        
        [pRobust,S]=robustfit(xVal(fitNucs),yVal(fitNucs));
        
        medNucs=abs(1-xVal./median(xVal,'omitnan'))<=0.05;
        medX=median(xVal(medNucs),'omitnan');
        medY=median(yVal(medNucs));        
        pSimple=[0 medY./medX];  
        %nucBin=sum(nucData(:,xId)>(haBins.*medX),2);
        nucBin=sum(xVal>(haBins),2);
        lineNucs=nucBin>0 & nucBin<numel(haBins);
        medLineX=accumarray(nucBin(lineNucs),xVal(lineNucs),[],@median);
        medLineY=accumarray(nucBin(lineNucs),yVal(lineNucs),[],@(x)quantile(x,0.475));
        pMid=fliplr(polyfit(medLineX,medLineY,1)); 
        pBlack=[0 median(yVal,'omitnan')./median(xVal,'omitnan')];
        pMix=[pMid(1) pMid(2)];
        
        pOrig=eval(fittingChoice);%fliplr(pMid)%;        
        adjX=xVal/median(xVal,'omitnan');
        adjY=(yVal-pOrig(1))./(pOrig(2)*median(xVal,'omitnan'));
        mycAdd(:,xId)=adjY-adjX;
        pKeep(xId,:)=pOrig;      
        if c==20
            c=0;
            figure;
        end
        c=c+1;        
        subplot(4,5,c) 
        nonNan=all([xVal yVal]>0,2);
        dscatter(xVal(nonNan),yVal(nonNan))
        hold on
        scatter(medLineX,medLineY,'ro','filled');
        plot([0 30],[0 30]*pOrig(2)+pOrig(1),'r-')
        plot([0 30],[0 30].*pExp(smeta.expId(xId),2)+pExp(smeta.expId(xId),1),'k-')        
        %pOrig=pExp(smeta.expId(xId),:)
        title(extractBetween(smeta.name(xId),'ha-','.mat'))
    else
        mycAdd(:,xId)=NaN;
    end
end
smeta.mycAdd(1:xId)=sum(mycAdd,'omitnan')';

% calcualte tor singal in individual samples
for xId=find(strcmp(smeta.ab,'ha')&ismember(smeta.expId,smpID))'
    yId=find(strcmp(smeta.name,strrep(smeta.name{xId},'ha','myc')));
    %abs(26-nucMeta.ToR)
end
%close all
figure
c=0
for e=smpID'
    c=c+1;
    subplot(2,2,1)
    selSmp=smeta.expId==e & smeta.bad==0 & smeta.mycAdd>0;
    scatter(c.*ones(sum(selSmp),1),smeta.mycAdd(selSmp),[],[.7 .7 .7],'o','filled')
    hold on
    scatter(c,sum(diff(adjColl{e},1,2),'omitnan'),'ko','filled')
    medianIntAll(e)=sum(diff(adjColl{e},1,2),'omitnan')
    medianInt(e)=median(smeta.mycAdd(selSmp));
    stdInt(e)=std(smeta.mycAdd(selSmp));
    title('cal. myc integration events')
    xlabel('Sample')
    set(gca,'XTick',1:numel(smpID),'XTickLabel',strcat(expTypes.tag2(smpID),'-',expTypes.gt(smpID)),'XTickLabelRotation',45);
    
    subplot(2,2,2)
    scatter(c.*ones(sum(selSmp),1),pKeep(selSmp,1),[],[.7 .7 .7],'o','filled')
    hold on
    scatter(c,pExp(e,1),'ko','filled')    
    title('"myc noise"')
    set(gca,'XTick',1:numel(smpID),'XTickLabel',strcat(expTypes.tag2(smpID),'-',expTypes.gt(smpID)),'XTickLabelRotation',45);
    
    
    subplot(2,2,3)
    scatter(c.*ones(sum(selSmp),1),pKeep(selSmp,2),[],[.7 .7 .7],'o','filled')
    hold on
    scatter(c,pExp(e,2),'ko','filled')    
    title('baseline slope')
    set(gca,'XTick',1:numel(smpID),'XTickLabel',strcat(expTypes.tag2(smpID),'-',expTypes.gt(smpID)),'XTickLabelRotation',45);
    
    
    subplot(2,2,4)    
    scatter(sum(diff(adjColl{e},1,2),'omitnan'),median(smeta.mycAdd(selSmp)),[],[.7 .7 .7],'o','filled')
    text(sum(diff(adjColl{e},1,2),'omitnan'),median(smeta.mycAdd(selSmp)),strjoin(expTypes{e,[2,4]}))

    hold on
    xlabel('median')
    ylabel('mean single sample')    
    %sumMycMedian(e)=sum(diff(adjColl{e},1,2),'omitnan');
    plot(xlim,xlim,'k--')
    title('calc myc events')
    sumMyc.(fittingChoice)(e)=sum(diff(adjColl{e},1,2),'omitnan');
end
suptitle(fittingChoice)
%buildtable
hisNucs=[2733,2734,12699,12700,54681,54682,54683];
hisProd=500;

ccFull=ccTime;
ccTime=0.75.*ccFull
c=0;
clear sumTable
for e=smpID'
    c=c+1;
    clear tableLine
    addMyc2=diff(adjColl{e},1,2);
    tableLine.mycInt=medianIntAll(e);
    tableLine.mycStd=stdInt(e);
    tableLine.mycMin=medianIntAll(e)/ccTime(e);
    tableLine.stdMin=stdInt(e)/ccTime(e);
    tableLine.mycFrac=(hisProd-tableLine.mycMin-tableLine.stdMin.*[-1 0 1])./hisProd;
    tableLine.totalMin=(tableLine.mycMin+tableLine.stdMin.*[-1 0 1])./tableLine.mycFrac;    
    
    hta1Events=mycAdd(hisNucs(2),smeta.expId==e & smeta.bad==0 & smeta.mycAdd>0)./tableLine.mycFrac(2)
    hta1Temp=ccTime(e)./(hta1Events)
    tableLine.hta1Time=mean(hta1Temp);
    tableLine.hta1Std=std(hta1Temp);    
    
    tdh3Events=median(mycAdd(nucMeta.gene==GP.gene_table.TDH3 & nucMeta.order>1,smeta.expId==e & smeta.bad==0 & smeta.mycAdd>0))./tableLine.mycFrac(2)
    tdh3Temp=ccTime(e)./(tdh3Events)
    tableLine.tdh3Time=mean(tdh3Temp);
    tableLine.tdh3Std=std(tdh3Temp); 
    
    tabelWStdLine=
    
    sumTable(c)=tableLine    
end

scatter(sumMyc.pRobust(smpID),sumMyc.pMid(smpID),'filled')
text(sumMyc.pRobust(smpID),sumMyc.pMid(smpID),strcat(expTypes.tag2(smpID),'-',expTypes.gt(smpID)))

title('median Integration in baseline fit')
xlabel('pRobust')
ylabel('pMid')
%% numbers
hisProduction=[0:10:1000];
mycEvents=[0:10:300]';
mycFraction=(hisProduction-mycEvents)./hisProduction;
subplot(1,2,2)
imagesc(mycFraction,'XData',[-0.5 1000.5],'YData',[-.5 300.5],[0 1])
xlabel('Histone (myc) production')
ylabel('Intergation events (that result in cleavage)')
ylabel(colorbar(),'myc Fraction')
hold on
plot([500 500],ylim(),'k-','DisplayName','20% hisPro of Sphase -Yoav');
scatter([500 500],[205 187 ],'ro','filled')
text([500 500],[205 187 ],{'HHT1','HHT2'})
title('steady state model of histone tunrover')

%% compare additional turnover between H2A and 

clear all
load('forFigure15')
load('./extData/addData.mat','h3To','h3AdjMid')
GP=load('group_imp.mat');
h3frac=0.60;

xId=find(contains(expTypes.ab,'ha')&contains(expTypes.gt,'wt')&contains(expTypes.con,'async')&strcmp(expTypes.tag2,'h2b'))
yId=find(contains(expTypes.ab,'myc')&contains(expTypes.gt,'wt')&contains(expTypes.con,'async')&strcmp(expTypes.tag2,'h2b'))
hisNucs=find(nucMeta.order<0 & ismember(nucMeta.gene,GP.groups{23}{2}{45}) & h3To>2);

fastTo=find(h3To>2 &nucMeta.order<=-1 &medNuc(:,yId)<=max(medNuc(hisNucs,yId)))
figure;
turnOverFromAsync(medNuc(:,16),medNuc(:,70),'axis1',subplot(2,2,1),'axis2',subplot(2,2,2),'haBins',[7.5:1:16.5])
subplot(2,2,1)
title('H2B wt async cr: 0.59')
xlabel('HA level ChipSeq')\
ylabel('myc level ChipSeq')
ylabel(colorbar(),'nucleosome density')

temp=load('forFigure3add.mat','smeta','nucData')
xId=find(contians(temp.smeta.name)

subplot(2,2,3)
hold off
nonNan=~isnan(medNuc(:,yId))&~isnan(h3AdjMid(:,2));
dscatter(h3AdjMid(nonNan,2)./h3frac,medNuc(nonNan,yId))
hold on
scatter(h3AdjMid(hisNucs,2)./h3frac,medNuc(hisNucs,yId),'ro','filled')
p=robustfit(h3AdjMid(hisNucs,2)./h3frac,medNuc(hisNucs,yId))
plot([0 max(h3AdjMid(hisNucs,2))]./h3frac,[0 max(h3AdjMid(hisNucs,2))]./h3frac.*p(2)+p(1),'k-','Linewidth',2)
axis tight
xlabel('H3 incorperation events')
ylabel('H2B myc Levels')
yLim=ylim()
yyaxis right
ylim((yLim-p(1))./p(2));
ylabel('estiamted H2B incorperation events')

subplot(2,2,4)
load('./extData/addData.mat','geneLvl')
selNucs=nucMeta.order>0 & nonNan& ismember(nucMeta.gene,find(~isnan(geneLvl)));
dscatter(geneLvl(nucMeta.gene(selNucs))-5.1,medNuc(selNucs,yId))
hold on
midNucs=false(size(selNucs))
midNucs(selNucs)=abs(geneLvl(nucMeta.gene(selNucs))-median(geneLvl,'omitnan'))<0.25
scatter(median(geneLvl(nucMeta.gene(midNucs))-5.1),median(medNuc(midNucs,yId)),'ro','filled')
transFac=2.^(median(geneLvl(nucMeta.gene(midNucs))-5.1))./median(medNuc(midNucs,yId));

axis tight
yLim=ylim()
ylim(yLim)
yyaxis right
ylim(yLim.*transFac);
ylabel('est.H2B incoperation events based on transcription')

nonNan=~isnan(adjColl{xId}(:,2))&~isnan(medNuc(:,yId))
dscatter(adjColl{xId}(nonNan,2),medNuc(nonNan,yId))
%scatter(adjColl{xId}(nonNan,2),medNuc(nonNan,yId),[],xTo(nonNan),'.')
hold on
scatter(adjColl{xId}(fastTo,2),medNuc(fastTo,yId),[],'o','filled','DisplayName','histone promoter')
scatter(adjColl{xId}(hisNucs,2),medNuc(hisNucs,yId),[],'o','filled','DisplayName','histone promoter')
hold on
plot(xlim,xlim.*p(2)+p(1),'k-','LineWidth',2)
xlabel('H3 integration events outside sphase')
ylabel(strjoin(expTypes{yId,1:4}))
subplot(1,2,2)
adjY=(medNuc(:,yId)-p(1))./p(2)
scatter(adjColl{xId}(nonNan,2),adjY(nonNan),'.')
title(sum(adjY,'omitnan'))

nInt=sum(medNuc(:,yId),'omitnan')./p(1)/100;
title(sprintf('Model based: %d int/min, ; myc/ha: %.2f%%',round(nInt),600/(600+nInt)*100))
h2bfrac=600/(600+nInt);

% subplot(1,2,2)
% seqH2=1.06/318;
% seqH3=1.46/316;
% scatter(medNuc(:,xId+54)*seqH3,medNuc(:,yId)*seqH2,'.')
% hold on
% xlabel(sprintf('abs. %s',strjoin(expTypes{xId+54,1:4})))
% ylabel(sprintf('abs. %s',strjoin(expTypes{yId,1:4})))
% title('absolute myc/ha from Spike-In')
% scatter(medNuc(hisNucs,xId+54)*seqH3,medNuc(hisNucs,yId)*seqH2,'ro','filled')
% plot(xlim,xlim*h2bfrac/h3frac,'k-')

%% dont know
load('qPCRtargets.mat','qPCRtargets')
figure
intHa=find(contains(expTypes.ab,'ha')&contains(expTypes.gt,'wt')&contains(expTypes.con,'async')&ismember(expTypes.tag2,{'h3clean','h2b','h3x2'}))
lowNucs={'KOG1','AUS1'}
for i=1:numel(intHa)
    myci=find(ismember(expTypes(:,2:4),expTypes(intHa(i),2:4)) & strcmp(expTypes.ab,'myc'))
    tempTo=log2(medNuc(:,myci)+1)-log2(medNuc(:,intHa(i))+1);
    subplot(1,3,i)
    hold off
    [y,bins]=histcounts(tempTo,30)
    plot(movmean(bins,2,'omitnan','Endpoints','discard'),y,'Linewidth',2)
    xlabel(sprintf('log2(%s+1)-log2(%s+1)',strjoin(expTypes{intHa(i),1:4}),strjoin(expTypes{myci,1:4})))
    hold on
    scatter(tempTo(qPCRtargets.nuc(qPCRtargets.Var13==1),:),ones(sum(qPCRtargets.Var13==1),1).*max(ylim)*0.2,'filled')
    text(tempTo(qPCRtargets.nuc(qPCRtargets.Var13==1),:),ones(sum(qPCRtargets.Var13==1),1).*max(ylim)*0.2,qPCRtargets.gene(qPCRtargets.Var13==1),'Rotation',45)
    axis tight
    title(strjoin(expTypes{intHa(i),2:4}))
    repTo=mean(tempTo(qPCRtargets.nuc(ismember(qPCRtargets.gene,lowNucs))));
    plot(repTo.*[1 1]+1,[0 ])
    frac=mean((tempTo(~isnan(tempTo))-repTo)>1)
    title(sprintf('%s - %.2f',strjoin(expTypes{intHa(i),2:4}),frac))
end


%% compare replication singals

expTypes.tp=max(str2double(regexp(expTypes.con,'(?<=T)\d+','match','once')),-1)
h3myc=contains(expTypes.ab,'myc') & contains(expTypes.tag2,'h3');
h2amyc=contains(expTypes.ab,'myc') & contains(expTypes.tag2,'h2a');
corrRep(h2amyc)=corr(nucMeta.ToR,medNuc(:,h2amyc),'rows','pairwise');
corrRep(h3myc)=corr(abs(26-nucMeta.ToR),medNuc(:,h3myc),'rows','pairwise');
[tcs,~,expTypes.tcId]=unique(expTypes(:,[1,4]))
figure
for t=[1,2]
    selSmp=find(expTypes.tcId==t);
    [tp,idx]=sort(expTypes.tp(selSmp))
    plot(expTypes.tp(selSmp(idx)),pExp(selSmp(idx)),'-','LineWidth',2,'DisplayName',sprintf('baseline %s',tcs.tag2{t}))
    hold on
end

for t=[3,4]
    selSmp=find(expTypes.tcId==t);
    [tp,idx]=sort(expTypes.tp(selSmp))
    plot(expTypes.tp(selSmp(idx)),-corrRep(selSmp(idx)),'-','LineWidth',2,'DisplayName',sprintf('Tor Signal %s',tcs.tag2{t}))
    hold on
end
axis tight
xlabel('time after release')
ylabel('corr/slope')

%% aanlyse replication arrested cells in timecourses
clear all
load('forFigure15')
load('./extData/addData.mat','h3To','h3AdjMid')
highTo=h3To>log(2);
close all
figure('Color',[1 1 1])
c=0;
c=c+1;
xVal=h3AdjMid(:,1);
yVal=h3AdjMid(:,2);
subplot(2,3,c)
scatter(xVal,yVal,[],diff(h3AdjMid,1,2),'o','filled')
caxis([0.5 5])
    ylabel('est. number of Myc \newline integration events p. cell cycle')
xlabel('HA./median(HA)')
title('Async Cells')
haBins=[7.5:1:16.5]/median(medNuc(:,34),'omitnan');
nucBin=sum(xVal>(haBins),2);
lineNucs=nucBin>0 & nucBin<numel(haBins);
medLineX=accumarray(nucBin(lineNucs),xVal(lineNucs),[],@median);
medLineY=accumarray(nucBin(lineNucs),yVal(lineNucs),[],@(x)quantile(x,0.475));
pMid=fliplr(polyfit(medLineX,medLineY,1));
hold on
baseLine=plot(xlim(),xlim()*pMid(2)+pMid(1),'k-','LineWidth',2,'DisplayName',sprintf('BaseLine: %.2fx%+.2f',pMid(2),pMid(1)));
set(gca,'FontSize',15)

c=c+1;
alphaHA=30;
alphaMyc=30+54;
nonNan=~isnan(medNuc(:,alphaMyc))&~isnan(diff(h3AdjMid,1,2));
p=linortfit2(medNuc(highTo&nonNan,alphaMyc),diff(h3AdjMid(highTo&nonNan,:),1,2));
normFac=p(1)

xVal=medNuc(:,alphaHA)./median(medNuc(:,alphaHA),'omitnan');
yVal=medNuc(:,alphaMyc).*normFac;
subplot(2,3,c)
scatter(xVal,yVal,[],diff(h3AdjMid,1,2),'o','filled')
caxis([0.5 5])
    ylabel('est. number of Myc \newline integration events')
xlabel('HA./median(HA)')
title('Alpha Factor')
haBins=[7.5:1:16.5]/median(medNuc(:,alphaHA),'omitnan');
nucBin=sum(xVal>(haBins),2);
lineNucs=nucBin>0 & nucBin<numel(haBins);
medLineX=accumarray(nucBin(lineNucs),xVal(lineNucs),[],@median);
medLineY=accumarray(nucBin(lineNucs),yVal(lineNucs),[],@(x)quantile(x,0.475));
pMid=fliplr(polyfit(medLineX,medLineY,1));
hold on
baseLine=plot(xlim(),xlim()*pMid(2)+pMid(1),'k-','LineWidth',2,'DisplayName',sprintf('BaseLine: %.2fx%+.2f',pMid(2),pMid(1)));
set(gca,'FontSize',15)
    
    

intExps={'forFigure3add.mat','forFigure14.mat','forFigure17.mat'}
intSampleNames={'T38','fast1a-wt-nacl33-x17.mat','wt-37t50'}
titleString={'Oxidative Stress (H2O2)','NaCl','HeatShock (37C)'}
highTo=h3To>log(2);

%% stresss-depednent arrest

for i=1:3
    c=c+1;
    subplot(2,3,c)
    temp=load(intExps{i},'smeta','nucData')
    tempSmpHa=find(contains(temp.smeta.name,intSampleNames{i}) & strcmp(temp.smeta.ab,'ha') & contains(temp.smeta.tag2,{'h3','h3clean'}));
    tempSmpMyc=find(contains(temp.smeta.name,intSampleNames{i}) & strcmp(temp.smeta.ab,'myc') & contains(temp.smeta.tag2,{'h3','h3clean'}));    
    nonNan=~isnan(temp.nucData(:,tempSmpMyc))&~isnan(diff(h3AdjMid,1,2));
    p=linortfit2(temp.nucData(highTo&nonNan,tempSmpMyc),diff(h3AdjMid(highTo&nonNan,:),1,2));
    normFac=p(1)  ;  
    xVal=temp.nucData(:,tempSmpHa)./median(temp.nucData(:,tempSmpHa),'omitnan');
    yVal=temp.nucData(:,tempSmpMyc).*normFac;    
    scatter(xVal,yVal,[],diff(h3AdjMid,1,2),'o','filled')
    caxis([0.5 5])
    ylabel('est. number of Myc \newline integration events')
    xlabel('HA./median(HA)')
    title(titleString{i})
    haBins=[7.5:1:16.5]/median(temp.nucData(:,tempSmpHa),'omitnan');
    nucBin=sum(xVal>(haBins),2);
    lineNucs=nucBin>0 & nucBin<numel(haBins);
    medLineX=accumarray(nucBin(lineNucs),xVal(lineNucs),[],@median);
    medLineY=accumarray(nucBin(lineNucs),yVal(lineNucs),[],@(x)quantile(x,0.475));
    pMid=fliplr(polyfit(medLineX,medLineY,1));
    hold on
    baseLine=plot(xlim(),xlim()*pMid(2)+pMid(1),'k-','LineWidth',2,'DisplayName',sprintf('BaseLine: %.2fx%+.2f',pMid(2),pMid(1)));
    set(gca,'FontSize',15)
end
ylabel(colorbar(),'estimated non Sphase exch. events')
save_gf(gcf,'SN1')


%% replciation arrest H2A
figure('Position',[1852 374 1892 420],'Renderer','painters','Color',[1 1 1]);
subplot(1,3,1)
hold off
goodNucs=all(~isnan(medNuc(:,[15,69])),2);
turnOverFromAsync(medNuc(:,15),medNuc(:,69),'axis1',subplot(1,3,1),'axis2',subplot(1,3,2),'haBins',[5.5:1:14.5])
subplot(1,3,1)
%title(sprintf('H2A wt async (cr:%.2f)',corr(medNuc(goodNucs,15),medNuc(goodNucs,69))))
title('H2A Async Cells')
xlabel('HA level ChipSeq')
ylabel('myc level ChipSeq')
ylabel(colorbar(),'nucleosome density','FontSize',10)
xlim([0 40])
ylim([0 40])
temp=load('forFigure3add.mat','smeta','nucData')
xId=find(strcmp(temp.smeta.name,'ha-h2a1-wt-h202T38-x19.mat'));
yId=find(strcmp(temp.smeta.name,'myc-h2a1-wt-h202T38-x19.mat'));
subplot(1,3,2)
hold off
goodNucs=all(~isnan(temp.nucData(:,[xId,yId])),2);
turnOverFromAsync(temp.nucData(:,xId),temp.nucData(:,yId),'axis1',subplot(1,3,2),'axis2',subplot(1,3,3),'haBins',[5.5:1:14.5])
subplot(1,3,2)
xlim([0 40])
ylim([0 40])
%title(sprintf('H2A wt H2O2 (cr:%.2f)',corr(temp.nucData(goodNucs,xId),temp.nucData(goodNucs,yId))))
title('H2A Oxidative stress (H2O2)')

xlabel('HA level ChipSeq')
ylabel('myc level ChipSeq')
ylabel(colorbar(),'nucleosome density','FontSize',10)
subplot(1,3,3)
hold off
imagesc(corr([medNuc(:,[15,69]),temp.nucData(:,[xId,yId])],'rows','complete'))
ylabel(colorbar(),'Pearson Correlation','FontSize',10)
caxis([0 1])
set(gca,'xtick',1:4,'ytick',1:4,'XTickLabel',{'Async HA','Async myc','H2O2 HA','H2O2 myc'},'YTickLabel',{'Async HA','Async myc','H2O2 HA','H2O2 myc'},'XTickLabelRotation',0)
title('H2A all nucleoosmes')


