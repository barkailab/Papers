%%Figure 2A
clear all
smFish=dir('/home/labs/barkailab/felixj/Documents/MATLAB/projects/shio/*/*smFISH_38*f6.mat')
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/shio/extData/dorsalGradient.mat')
close all
figure('Units','normalized','Position',[0.0073 0.4019 1 0.4139],'Color',[1 1 1])
c=0;
intSmp=[2]
paperData{1,1}=smFish(intSmp,:);
paperData{1,2}='Fig2A';
for i=intSmp
    clear xLimLine
    c=c+1;
    matFile=[smFish(i).folder,'/' smFish(i).name];
    out=load(matFile);
    midPoint=mean(out.midLine);
    rotAngle=180-atand((out.midLine(1,1)-out.midLine(end,1))./(out.midLine(1,2)-out.midLine(end,2)));    
    midLineRot=nan(size(out.midLine));
    midLineRot(:,1)=[(out.midLine(:,1)-midPoint(1))*cosd(rotAngle) + (out.midLine(:,2)-midPoint(2))*sind(rotAngle)];
    midLineRot(:,2)=[-(out.midLine(:,1)-midPoint(1))*sind(rotAngle) + (out.midLine(:,2)-midPoint(2))*cosd(rotAngle)];      
    nucDis=sqrt(polyarea(out.roiPts.X,out.roiPts.Y)/min(400,sum(inpolygon(out.nuclei.Centroid(:,1),out.nuclei.Centroid(:,2),out.roiPts.X,out.roiPts.Y))));
    disBins=[-20.5:20.5]*nucDis;
    for j=find(contains(out.probes,'twist'))
        subplot(numel(intSmp),3,1)
        hold off
        locs=out.ts{j};
        if isfield(out,'roiPts')
            intLocs=inpolygon(locs.Centroid(:,1),locs.Centroid(:,2),out.roiPts.X,out.roiPts.Y);
        else
            intLocs=true(size(locs,1),1);
        end        
        locs.Xnew=[(locs.Centroid(:,1)-midPoint(1))*cosd(rotAngle) + (locs.Centroid(:,2)-midPoint(2))*sind(rotAngle)];
        locs.Ynew=[-(locs.Centroid(:,1)-midPoint(1))*sind(rotAngle) + (locs.Centroid(:,2)-midPoint(2))*cosd(rotAngle)];
        locs.binIdx=sum(locs.midDis>disBins,2)

        scatter(locs.Xnew(intLocs)/nucDis,locs.Ynew(intLocs)/nucDis,[],locs.int(intLocs),'filled','LineWidth',0.5,'MarkerEdgeColor',[.4 .4 .4])
        axis tight
        caxis(quantile(locs.int(intLocs),[0.1 0.9]))
        if contains(out.probes{j},'twist')
            tssColor=brewermap(64,'Greens');
            tssColor=tssColor(1:50,:)
        else
            tssColor=brewermap(64,'OrRd');
            tssColor=tssColor(1:50,:)
        end
        colormap(gca,tssColor)
        yLim=ylim();
        xLim=xlim();
        hold on
        plot(midLineRot(:,1)/nucDis,midLineRot(:,2)/nucDis,'k-','LineWidth',2)
        xlim(xLim+[- 1 1]);
        ylim(yLim);
        title(sprintf('%s,%s,%d,%.2f',strrep(extractBefore(smFish(i).name,'f6.mat'),'_',' '),out.probes{j},sum(intLocs),nucDis))
        %xlabel('Position relative to midline')
        yticks([])
        ylabel(colorbar(),'TS spot intensity')
        set(gca,'BoxStyle','full')
        xlim(ceil(max(abs(locs.midDis(locs.binIdx>0 & locs.binIdx<42 & intLocs)))/nucDis+1).*[-1 1])
        
        subplot(numel(intSmp),3,2)
        hold off
        yyaxis left
        plot(movmean(disBins,2,'Endpoints','discard')/nucDis,accumarray(locs.binIdx(locs.binIdx>0 & locs.binIdx<42 & intLocs),1,[41 1]),'DisplayName',out.probes{j},'LineWidth',2,'Color',tssColor(33,:))
        set(gca,'YColor',[0 0 0])
        yyaxis right
        plot([-14:15],doralLvl,'k--','DisplayName','dorsal gradient')
        set(gca,'YColor',[0 0 0])
        %xlabel('Position relative to midline')
        xlim(ceil(max(abs(locs.midDis(locs.binIdx>0 & locs.binIdx<42 & intLocs)))/nucDis+1).*[-1 1])
        subplot(numel(intSmp),3,3)
        hold off
        plot(movmean(disBins,2,'Endpoints','discard')/nucDis,accumarray(locs.binIdx(locs.binIdx>0 & locs.binIdx<42 & intLocs),locs.int(locs.binIdx>0 & locs.binIdx<42 & intLocs),[41 1],@median),'DisplayName',out.probes{j},'LineWidth',2,'Color',tssColor(33,:))
        hold on
        xLimLine{j}=quantile(locs.midDis(intLocs),[0 1])/nucDis   
        yyaxis right
        plot([-14:15],doralLvl,'k--','DisplayName','dorsal gradient')
        set(gca,'YColor',[0  0 0])
        %xlabel('Position relative to midline')        
        xlim(ceil(max(abs(locs.midDis(locs.binIdx>0 & locs.binIdx<42 & intLocs)))/nucDis+1).*[-1 1])
    end
end
save_gf(gcf,'Figure2A_T48Twist','type',{'svg','jpg'},'paper','benny21')



%% Figure 2B
clearvars -except paperData
smFish=dir('/home/labs/barkailab/felixj/Documents/MATLAB/projects/shio/*/*smFISH_38*f6.mat')
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/shio/extData/dorsalGradient.mat')
close all
figure('Units','normalized','Position',[[0,0,1,1]],'Color',[1 1 1])
c=0;
intSmp=[4,3,5]
paperData{2,1}=smFish(intSmp,:)
paperData{2,2}='Fig2B'
for i=intSmp
    clear xLimLine
    c=c+1;
    matFile=[smFish(i).folder,'/' smFish(i).name];
    out=load(matFile);
    midPoint=mean(out.midLine);
    rotAngle=180-atand((out.midLine(1,1)-out.midLine(end,1))./(out.midLine(1,2)-out.midLine(end,2)));    
    midLineRot=nan(size(out.midLine));
    midLineRot(:,1)=[(out.midLine(:,1)-midPoint(1))*cosd(rotAngle) + (out.midLine(:,2)-midPoint(2))*sind(rotAngle)];
    midLineRot(:,2)=[-(out.midLine(:,1)-midPoint(1))*sind(rotAngle) + (out.midLine(:,2)-midPoint(2))*cosd(rotAngle)];      
    nucDis=sqrt(polyarea(out.roiPts.X,out.roiPts.Y)/sum(inpolygon(out.nuclei.Centroid(:,1),out.nuclei.Centroid(:,2),out.roiPts.X,out.roiPts.Y)));
    disBins=[-20.5:20.5]*nucDis;
    for j=find(~cellfun('isempty',out.ts))
        subplot(3,numel(intSmp),(j-1)*numel(intSmp)+c)
        hold off
        locs=out.ts{j};
        if isfield(out,'roiPts')
            intLocs=inpolygon(locs.Centroid(:,1),locs.Centroid(:,2),out.roiPts.X,out.roiPts.Y);
        else
            intLocs=true(size(locs,1),1);
        end        
        locs.Xnew=[(locs.Centroid(:,1)-midPoint(1))*cosd(rotAngle) + (locs.Centroid(:,2)-midPoint(2))*sind(rotAngle)];
        locs.Ynew=[-(locs.Centroid(:,1)-midPoint(1))*sind(rotAngle) + (locs.Centroid(:,2)-midPoint(2))*cosd(rotAngle)];
        scatter(locs.Xnew(intLocs)/nucDis,locs.Ynew(intLocs)/nucDis,[],locs.int(intLocs),'filled','LineWidth',0.5,'MarkerEdgeColor',[.4 .4 .4])
        axis tight
        caxis(quantile(locs.int(intLocs),[0.1 0.9]))
        if contains(out.probes{j},'twist')
            tssColor=brewermap(64,'Greens');
            tssColor=tssColor(1:50,:)
        else
            tssColor=brewermap(64,'OrRd');
            tssColor=tssColor(1:50,:)
        end
        colormap(gca,tssColor)
        yLim=ylim();
        xLim=xlim();
        hold on
        plot(midLineRot(:,1)/nucDis,midLineRot(:,2)/nucDis,'k-','LineWidth',2)
        xlim([-14 14]);
        ylim(yLim);
        title(sprintf('%s,%s,%d,%.2f',strrep(extractBefore(smFish(i).name,'f6.mat'),'_',' '),out.probes{j},sum(intLocs),nucDis))
        %xlabel('Position relative to midline')
        yticks([])
        ylabel(colorbar(),'TS spot intensity')
        set(gca,'BoxStyle','full')
        
        subplot(3,numel(intSmp),2*numel(intSmp)+c)
        locs.binIdx=sum(locs.midDis>disBins,2)
        plot(movmean(disBins,2,'Endpoints','discard')/nucDis,accumarray(locs.binIdx(locs.binIdx>0 & locs.binIdx<42 & intLocs),1,[41 1]),'DisplayName',out.probes{j},'LineWidth',2,'Color',tssColor(33,:))
        hold on
        xLimLine{j}=quantile(locs.midDis(intLocs),[0 1])/nucDis
    end
    subplot(3,numel(intSmp),2*numel(intSmp)+c)
    plot([-14:15],rescale(doralLvl,min(ylim()),max(ylim)),'k--','DisplayName','dorsal gradient')
    out.nuclei.midDis=disToLines(out.nuclei.Centroid(:,1),out.nuclei.Centroid(:,2),out.midLine(:,1),out.midLine(:,2));
    out.nuclei.binIdx=sum(out.nuclei.midDis>disBins,2);
    intNucs=inpolygon(out.nuclei.Centroid(:,1),out.nuclei.Centroid(:,2),out.roiPts.X,out.roiPts.Y)
    plot(movmean(disBins,2,'Endpoints','discard')/nucDis,accumarray(out.nuclei.binIdx(out.nuclei.binIdx>0 & out.nuclei.binIdx<42 & intNucs),1,[41 1]),'k-','DisplayName','nuclei','LineWidth',2)
    xLimLine=cat(1,xLimLine{:})
    %xlim([min(xLimLine(:,1)) max(xLimLine(:,2))]+[-1 1])
    xlim([-14 14])
    %title('TS number')
end
save_gf(gcf,'Figure2B_T48Twist','type',{'svg','jpg'},'paper','benny21')
%% figure 3
clearvars -except paperData

clear all
load('dorsalGradient.mat','doralLvl');
countData=readtable('shariCounting.xlsx','ReadVariableNames',false)
%remove empty lines 
countData(cellfun('isempty',countData.Var1)&cellfun('isempty',countData.Var2),:)=[];
embryoRow=[find(~cellfun('isempty',countData.Var1));size(countData,1)]
intCols=false(1,size(countData,2));
for i=1:numel(intCols)
    intCols(i)=isa(countData{:,i},'double')
end
close all
probeColor=       brighten([brewermap(64,'Reds');brewermap(1,'Greens')],-0.2);
probeColor=probeColor([40,end],:)
figure('Units','normalized','Position',[[0,0,1,1]],'Color',[1 1 1])
for i=1:numel(embryoRow)-1
    tempRows=ismember([1:size(countData,1)]',embryoRow(i):embryoRow(i+1)-1) & ~cellfun('isempty',countData.Var2);
    tempTable=countData(tempRows,:)
    probeId=regexp(tempTable.Var2,'(.*)\s([0,1,2])','tokens','once');
    tempTable(cellfun('isempty',probeId),:)=[]
    probeId=cat(1,probeId{:});
    [uProbes,~,uid]=unique(probeId(:,1));
    for u=1:numel(uProbes)
        subplot(3,numel(embryoRow)-1,i)
        primedFraction=sum(tempTable{uid==u,intCols}.*str2double(probeId(uid==u,2)))./sum(tempTable{uid==u,intCols})/2;
        primedFraction=primedFraction(~isnan(primedFraction))
        nCells=(sum(~isnan(primedFraction))-1)/2;
        plot(-nCells:nCells,movmean(primedFraction(~isnan(primedFraction)),1),'LineWidth',2,'DisplayName',uProbes{u},'Color',probeColor(u,:))
        hold on
        if ~strcmp(uProbes{u},'twist')
            subplot(3,numel(embryoRow)-1,i+numel(embryoRow)-1)
            locDorsal=interp1(-14:15,doralLvl,-nCells:nCells)./max(doralLvl);
            scatter(locDorsal,primedFraction(~isnan(primedFraction)),[],probeColor(u,:),'filled','DisplayName',uProbes{u})
            pLine=robustfit(locDorsal,primedFraction(~isnan(primedFraction)));            
            xlabel('Dorsal level [AU]')
            ylabel('primed fraction')
            hold on
            plot(quantile(locDorsal,[0 1]),pLine(2)*quantile(locDorsal,[0 1])+pLine(1),'-','LineWidth',1.5,'Color',brighten(probeColor(u,:),-0.4))
            %text(mean(quantile(locDorsal,[0 1])),pLine(2)*mean(quantile(locDorsal,[0 1]))+pLine(1),sprintf('cr: %.2f',corr(locDorsal',primedFraction(~isnan(primedFraction))','rows','pairwise')),'Color',brighten(probeColor(u,:),-0.4))
            subplot(3,numel(embryoRow)-1,i+2*(numel(embryoRow)-1))
            p=-log(1-primedFraction)/15;
            scatter(locDorsal,p(~isnan(primedFraction)),[],probeColor(u,:),'filled','DisplayName',uProbes{u})
            ylabel('priming prob/min')
            xlabel('Dorsal level [AU]')
            hold on
            pLine=robustfit(locDorsal,p(~isnan(primedFraction))); 
            plot(quantile(locDorsal,[0 1]),pLine(2)*quantile(locDorsal,[0 1])+pLine(1),'-','LineWidth',1.5,'Color',brighten(probeColor(u,:),-0.4))
           % text(mean(quantile(locDorsal,[0 1])),pLine(2)*mean(quantile(locDorsal,[0 1]))+pLine(1),sprintf('cr: %.2f',corr(locDorsal',p(~isnan(primedFraction))','rows','pairwise')),'Color',brighten(probeColor(u,:),-0.4))           
        end
    end
    subplot(3,numel(embryoRow)-1,i)
    plot(-14:15,rescale(doralLvl,0,1),'k--','LineWidth',2,'DisplayName','dorsal gradient')
    xlabel('Distance from midLine')
    ylabel('Fraction primed promoters')
    title(countData.Var1{embryoRow(i)})
    xlim(nCells.*[-1 1])
    %legend()
end
save_gf(gcf,'Figure3_T48primingProb','type',{'svg','jpg'},'paper','benny21')

%% figure 4
clearvars -except paperData
smFish=dir('/home/labs/barkailab/felixj/Documents/MATLAB/projects/shio/*/*smFISH_35*f6.mat')
selNames={'with 3''','w/o 3'''}
load('dorsalGradient.mat','doralLvl');
close all
figure('Units','normalized','Position',[[0,0,1,1]],'Color',[1 1 1]);
cls=lines(3);
cls=[[0 .75 0];brighten(brewermap(1,'Reds'),-0.6)];
paperData{3,1}=smFish(1,:)
paperData{3,2}='Fig4'

for i=1:numel(smFish)
    out=load([smFish(i).folder '/' smFish(i).name])
    nucDis=sqrt(polyarea(out.roiPts.X,out.roiPts.Y)/sum(inpolygon(out.nuclei.Centroid(:,1),out.nuclei.Centroid(:,2),out.roiPts.X,out.roiPts.Y)));
    fiveP=find(contains(out.probes,'-5'));
    threeP=find(contains(out.probes,'-3'));
    if isfield(out,'roiPts')
        sel3P=inpolygon(out.ts{threeP}.Centroid(:,1),out.ts{threeP}.Centroid(:,2),out.roiPts.X,out.roiPts.Y);
        locs3P=out.ts{threeP}(sel3P,:);
        sel5P=inpolygon(out.ts{fiveP}.Centroid(:,1),out.ts{fiveP}.Centroid(:,2),out.roiPts.X,out.roiPts.Y);
        locs5P=out.ts{fiveP}(sel5P,:);
        midWidth=1000;
        clear sel3P sel5P
    else
        locs5P=out.ts{fiveP}(all(out.ts{fiveP}.Centroid>20,2),:);
        locs3P=out.ts{threeP}(all(out.ts{threeP}.Centroid>20,2),:);
        midWidth=500
    end
    disMat=pdist2(locs5P.Centroid,locs3P.Centroid);
    [val,idx]=min(disMat,[],1);
    idx(val>5)=0;
    for j=find(accumarray(idx(idx>0)',1)>1)'
        sel3P=idx==j;
        idx(sel3P)=idx(sel3P).*(val(sel3P)==min(val(sel3P)));
    end
    locs5P.has3P(idx(idx>0))=true;
    locs5P.int3P(idx(idx>0))=locs3P.int(idx>0);
    clear disMat idx    

    disBins=nucDis.*[-15.5:1:15.5];
    binMark=nucDis.*[-16:1:16];
    binIdx=sum(locs5P.midDis>=disBins,2)+1;
    subplot(3,numel(smFish),i)
    hold off
    subplot(3,numel(smFish),i+numel(smFish))
    hold off
    c=0;
    for j=[true,false]
        c=c+1;
        sel5P=locs5P.has3P==j;
        binMedian=accumarray(binIdx(sel5P),locs5P.int(sel5P),[numel(disBins)+1,1],@median,NaN);
        binStd=accumarray(binIdx(sel5P),locs5P.int(sel5P),[numel(disBins)+1,1],@std,NaN);
        binN=accumarray(binIdx(sel5P),1,[numel(disBins)+1,1],@sum);        
        subplot(3,numel(smFish),i)
        plot(binMark/nucDis,smoothdata(binN,'gaussian',3),'Color',cls(c,:),'LineWidth',2,'DisplayName',sprintf('%s',selNames{c}))
        hold on
        if j==true
            crInd=corr(locs5P.int(sel5P),interp1((-14:15).*nucDis,doralLvl,locs5P.midDis(sel5P)),'rows','pairwise');
            crMed=corr(binMedian(3:end-1),doralLvl,'rows','pairwise');
            subplot(3,numel(smFish),i+numel(smFish))
            lineColl(c)=plot(binMark/nucDis,smoothdata(binMedian,'gaussian',3),'Color',cls(2,:),'LineWidth',2,'DisplayName',sprintf('%s,%.2f / %.2f',selNames{c},crInd,crMed));  
            hold on
            %plot(binMark/nucDis,smoothdata(binMedian,'gaussian',3)-0.02.*range(ylim),'--','Color',[0 .75 0],'LineWidth',2,'DisplayName',sprintf('%s,%.2f / %.2f',selNames{c},crInd,crMed));
            binDorsal=interp1(-14:15,doralLvl,binMark/nucDis);
            subplot(3,numel(smFish),i+numel(smFish)*2)
            scatter(binDorsal,binMedian,[],cls(2,:),'o','filled')
            hold on
            p=robustfit(binDorsal,binMedian);
            
            plot(xlim,xlim*p(2)+p(1),'-','Color',cls(2,:))
            %scatter(interp1((-14:15).*nucDis,doralLvl,locs5P.midDis(sel5P)),locs5P.int(sel5P),[],cls(2,:),'.')
            %p=robustfit(interp1((-14:15).*nucDis,doralLvl,locs5P.midDis(sel5P)),locs5P.int(sel5P))
            %plot(xlim,xlim*p(2)+p(1),'--','Color',cls(2,:))

        end
    end
    subplot(3,numel(smFish),i)
    xlim(binMark(quantile(binIdx,[0 1])+[-1 1])/nucDis)
    xlabel('Distance from MidLine')
    ylabel('Number of TSS')
    title(strrep(regexp(smFish(i).name,'\d+_\d+_\d+','match'),'_',' '))
    plot(-14:15,rescale(doralLvl,0,max(ylim)),'k--')
    %legend()
    subplot(3,numel(smFish),i+numel(smFish))
    xlabel('Distance from MidLine')
    xlim(binMark(quantile(binIdx(locs5P.has3P),[0 1])+[-1 1])/nucDis)
    plot(-14:15,rescale(doralLvl,min(ylim),max(ylim)),'k--')
    %legend()
    ylabel('Median Intensity')
end
set(gcf,'color','w')
save_gf(gcf,'Figure4A_T48loadingDashedLine','type',{'svg','jpg'},'paper','benny21')

%% Figure 5
clearvars -except paperData
close all
figure('Units','normalized','Position',[[0,0,1,1]],'Color',[1 1 1]);
smFish=dir('/home/labs/barkailab/felixj/Documents/MATLAB/projects/shio/*/*smFISH_35*f6.mat')
load('dorsalGradient.mat','doralLvl');
lineCol=brewermap(256,'greens')
lineCol=lineCol(156,:)
paperData{4,1}=smFish(2,:)
paperData{4,2}='Fig5'

for i=1:numel(smFish)
    out=load([smFish(i).folder '/' smFish(i).name])
    threeP=find(contains(out.probes,'3'));
    if ~isfield(out,'roiPts')
        out.roiPts=table([0;0;size(out.nucleiMask,2);size(out.nucleiMask,2)],[0;size(out.nucleiMask,1);size(out.nucleiMask,1);0],'VariableNames',{'X','Y'})
    elseif numel(out.roiPts)==0
        out.roiPts=table([0;0;size(out.nucleiMask,2);size(out.nucleiMask,2)],[0;size(out.nucleiMask,1);size(out.nucleiMask,1);0],'VariableNames',{'X','Y'})
    else
        nucDis=sqrt(polyarea(out.roiPts.X,out.roiPts.Y)/sum(inpolygon(out.nuclei.Centroid(:,1),out.nuclei.Centroid(:,2),out.roiPts.X,out.roiPts.Y)));
    end
    pixelWidth=size(out.nucleiMask,2);
    pixelheight=size(out.nucleiMask,1);
    tempPic=out.pic(:,:,threeP);
    subplot(3,numel(smFish),i)
    selPixel=inpolygon(repmat([1:pixelWidth],pixelheight,1),repmat([1:pixelheight]',1,pixelWidth),out.roiPts.X,out.roiPts.Y);%&out.nucleiMask==0;
    mRNALvl=accumarray(round(out.disMat(selPixel)-min(out.disMat(selPixel))+1),tempPic(selPixel),[],@median);
    hold on
    %   crMrna=corr(mRNALvl,interp1((-14:15).*nucDis,doralLvl,min(out.disMat(selPixel))+[0.5:numel(mRNALvl)-.5])','rows','pairwise');
    plot((min(out.disMat(selPixel))+[0.5:numel(mRNALvl)-.5])/nucDis,movmean(mRNALvl,50),'-','LineWidth',2,'DisplayName',sprintf('%s: %.2f mRNA',out.probes{threeP},0),'Color',lineCol)
    hold on
    xLim=[floor(min(out.disMat(selPixel))/nucDis), ceil(max(out.disMat(selPixel))/nucDis)];
    %selPixel=true(size(out.disMat));
    %mRNALvl=accumarray(round(out.disMat(selPixel)-min(out.disMat(selPixel))+1),tempPic(selPixel),[],@median);
    %plot((min(out.disMat(selPixel))+[0.5:numel(mRNALvl)-.5])/nucDis,movmean(mRNALvl,50),'--','LineWidth',2,'DisplayName',sprintf('%s: %.2f mRNA',out.probes{threeP},0),'Color',lineCol)
    axis tight
    xlim(xLim)
    plot((-14:15),rescale(doralLvl,min(ylim),max(ylim)),'k--')
    xlabel('Distance from midline')
    ylabel('median probe intensity')
    title(strrep(smFish(i).name,'_',' '))
    yLim=ylim()
    subplot(3,numel(smFish),i+numel(smFish))
    p=polyfit(out.midLine(:,1),out.midLine(:,2),1);
    rotDeg=90-atand(p(1));
    
    imageMat=imrotate(wiener2(tempPic,[100,100]),-rotDeg);
    hold off
    imagesc(imageMat,'XData',-quantile(out.disMat(:),[1 0])/nucDis,'Ydata',(1:size(imageMat,1))/nucDis,quantile(imageMat(:),[0,0.999]))
    hold on
    colormap(gca,brewermap(256,'greens'))
    colorbar()
    %plot(midLineRot(:,1),midLineRot(:,2),'k-')
    %     nPixel=sum(imageMat>0);
    %     [~,mPos]=max(nPixel);
    %     for c=1:min(mPos-1,numel(nPixel)-mPos)
    %         pixelArea(c)=(2*c+1).*min(nPixel(mPos+[-c:c]));
    %     end
    %     [~,bestC]=max(pixelArea)
    
    subplot(3,numel(smFish),i+2*numel(smFish))
    %     imagesc(imrotate(mRNALvl(round(out.disMat-min(out.disMat,[],'all'))+1),-rotDeg),'XData',-quantile(out.disMat(:),[1 0])/nucDis,'Ydata',(1:size(imageMat,1))/nucDis)
    %     colormap(gca,brewermap(256,'greens'))
    %     colorbar()
    crMrna=corr(interp1(-14:15,doralLvl,(min(out.disMat(selPixel))+[0.5:numel(mRNALvl)-.5])/nucDis)', mRNALvl,'rows','pairwise');
    scatter(interp1(-14:15,doralLvl,(min(out.disMat(selPixel))+[0.5:numel(mRNALvl)-.5])/nucDis), mRNALvl,[],lineCol,'.','DisplayName',num2str(crMrna,'%.2f'))  
    axis tight
    ylim([yLim(1),max(ylim)])
    xlim([min(xlim) 1])
    legend()
end
save_gf(gcf,'Figure5A_mRNAGreen','type',{'svg','jpg'},'paper','benny21')

close all
clear all
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/shio/extData/dorsalGradient.mat')
p = 0.13
dorsalGradient=movmean(doralLvl./max(doralLvl),3);

pGradient=p*doralLvl/max(doralLvl);
rnaLevel=@(p,T,a)(a*(T+(exp(-T.*p)-1)./p))
openFraction=@(p,T)((1-exp(-T.*p)))
plot(0:1:100,rnaLevel(0.1,0:1:100,1))
%dorsal=normpdf([-15:15]',0,5.5)/normpdf(0,0,5.5);
lineColor=parula(20);
close all
figure('Units','normalized','Position', [0.1755 0.4472 0.7385 0.3231],'Color',[1 1 1])
c=0;
pMid=0.13
for t=0:1:15
    subplot(1,3,1)
    lineObj2(t+1)=plot([-14.5:14.5]',openFraction(max(dorsalGradient-0.5,0)*pMid*2,t),'Color',lineColor(t+1,:),'LineWidth',2,'DisplayName',sprintf('mrna at %d',t))
    hold on
    xlim([-10 10])
    title('fraction of primed promoters')
    xlabel('relative position to midline')
    
    subplot(1,3,2)
    lineObj(t+1)=plot([-14.5:14.5]',max(0,rnaLevel(max(dorsalGradient-0.5,0)*pMid*2,t,1)),'Color',lineColor(t+1,:),'LineWidth',2,'DisplayName',sprintf('mrna at %d',t))
    hold on
    xlim([-10 10])
    title('mRNA w/o gradient')
    xlabel('relative position to midline')

    subplot(1,3,3)
    lineObj(t+1)=plot([-14.5:14.5]',max(0,rnaLevel(max(dorsalGradient-0.5,0)*pMid*2,t,1)).*dorsalGradient,'Color',lineColor(t+1,:),'LineWidth',2,'DisplayName',sprintf('mrna at %d',t))
    hold on
    xlim([-10 10])
    title('mRNA with int Gradient')  
    colormap(parula(16))
    ylabel(colorbar(),'minutes in cycle 14');    
    xlabel('relative position to midline')

    caxis([-.5 15.5])    
end
save_gf(gcf,'Figure5B_model','type',{'svg','jpg'},'paper','benny21')



%% Figure 6
clearvars -except paperData
smFish=dir('/home/labs/barkailab/felixj/Documents/MATLAB/projects/shio/*/*smFISH_50*f6.mat')
smFish=smFish([10,11,12])
load('dorsalGradient.mat','doralLvl');
close all
c=0;
figure('Units','normalized','Position',[[0,0,1,1]],'Color',[1 1 1])
lineCol=[ 0.8500    0.3250    0.0980];
paperData{5,1}=smFish([1,2],:)
paperData{5,2}='Fig6'

for i=find(contains({smFish.name},'mist'))
    c=c+1;
    out=load([smFish(i).folder '/' smFish(i).name])
    if ~isfield(out,'roiPts')
        out.roiPts=table([0;0;size(out.nucleiMask,2);size(out.nucleiMask,2)],[0;size(out.nucleiMask,1);size(out.nucleiMask,1);0],'VariableNames',{'X','Y'})
        nucDis=40;
    elseif numel(out.roiPts)==0
        out.roiPts=table([0;0;size(out.nucleiMask,2);size(out.nucleiMask,2)],[0;size(out.nucleiMask,1);size(out.nucleiMask,1);0],'VariableNames',{'X','Y'})
        nucDis=40;
    else
        nucDis = sqrt(polyarea(out.roiPts.X,out.roiPts.Y)/sum(inpolygon(out.nuclei.Centroid(:,1),out.nuclei.Centroid(:,2),out.roiPts.X,out.roiPts.Y)));
    end
    bins=-16:16
    subplot(3,sum(contains({smFish.name},'mist')),c)
    hold off
    subplot(3,sum(contains({smFish.name},'mist')),c+sum(contains({smFish.name},'mist')))
    hold off
    for j=find(contains(out.probes,'mist'))
        if ~ismember('midDis',out.ts{j}.Properties.VariableNames)
            out.ts{j}.midDis =disToLines(out.ts{j}.Centroid(:,1),out.ts{j}.Centroid(:,2),out.midLine(:,1),out.midLine(:,2));
        end
        locs=out.ts{j};
        [~,binIdx]=ismember(round(locs.midDis./nucDis),bins);        
        intLocs=inpolygon(locs.Centroid(:,1),locs.Centroid(:,2),out.roiPts.X,out.roiPts.Y) & binIdx>0;     
        crInd=corr(locs.int(intLocs),interp1((-14:15).*nucDis,doralLvl,locs.midDis(intLocs)),'rows','pairwise');
        nLoc=accumarray(binIdx(intLocs),1,[numel(bins),1],@sum,0);
        crN=corr(nLoc(3:end-1),doralLvl,'rows','pairwise');
        intLoc=accumarray(binIdx(intLocs),locs.int(intLocs),[numel(bins),1],@median,NaN);
        subplot(3,sum(contains({smFish.name},'mist')),c)
        plot(bins,smoothdata(nLoc,'gaussian',3),'LineWidth',2,'DisplayName',sprintf('%s: %.2f',out.probes{j},crN),'Color',lineCol)
        hold on
        subplot(3,sum(contains({smFish.name},'mist')),c+sum(contains({smFish.name},'mist')))
        plot(bins,smoothdata(intLoc,'gaussian',3),'LineWidth',2,'DisplayName',sprintf('%s: %.2f',out.probes{j},crInd),'Color',lineCol)
        axis tight
        %text(min(xlim),max(ylim),sprintf('%.2f',crInd),'Color',lineCol,'VerticalAlignment','top','HorizontalAlignment','left')
        hold on
        
        subplot(3,sum(contains({smFish.name},'mist')),c+2*sum(contains({smFish.name},'mist')))
        tempPic=out.pic(:,:,j);
        selPixel=inpolygon(repmat([1:1024],1024,1),repmat([1:1024]',1,1024),out.roiPts.X,out.roiPts.Y);%&out.nucleiMask==0;
        mRNALvl=accumarray(round(out.disMat(selPixel)-min(out.disMat(selPixel),[],'all')+1),tempPic(selPixel),[],@median);
        hold on        
        crMrna=corr(mRNALvl,interp1((-14:15).*nucDis,doralLvl,min(out.disMat(selPixel))+[0.5:numel(mRNALvl)-.5])','rows','pairwise');
        plot((min(out.disMat(selPixel))+[0.5:numel(mRNALvl)-.5])/nucDis,movmean(mRNALvl,50),'LineWidth',2,'DisplayName',sprintf('%s: %.2f mRNA',out.probes{j},crMrna),'Color',lineCol)
        axis tight
        %text(min(xlim),max(ylim),sprintf('%.2f',crMrna),'Color',lineCol,'VerticalAlignment','top','HorizontalAlignment','left')     
        
    end
    subplot(3,sum(contains({smFish.name},'mist')),c)
    xlim(bins(max(1,min(quantile(binIdx,[0 1])+[-1 1],numel(bins)))))
    plot((-14:15),rescale(doralLvl,min(ylim),max(ylim)),'k--','DisplayName','DorsalLvl') 
   % legend()
    xlabel('Distance from midline')
    ylabel('# of TSS')
    title(sprintf('%s %.2f',strrep(regexp(smFish(i).name,'\d+_.*?_\d+','match','once'),'_',' ' ),nucDis))
    %axis tight
    xlim([-15 15])

    
    subplot(3,sum(contains({smFish.name},'mist')),c+sum(contains({smFish.name},'mist')))
    xlim(bins(max(1,min(quantile(binIdx,[0 1])+[-1 1],numel(bins)))))
    plot((-14:15),rescale(doralLvl,min(ylim),max(ylim)),'k--','DisplayName','DorsalLvl')      
    %legend()    
    xlabel('Distance from midline')
    ylabel('median Intensity')
    %axis tight
    xlim([-15 15])
    
    subplot(3,sum(contains({smFish.name},'mist')),c+2*sum(contains({smFish.name},'mist')))
    plot((-14:15),rescale(doralLvl,min(ylim),max(ylim)),'k--','DisplayName','DorsalLvl')      
    %legend()    
    xlabel('Distance from midline')
    ylabel('median Cytoplamisc intensity (mRNA)')
    xlim(bins(max(1,min(quantile(binIdx,[0 1])+[-1 1],numel(bins)))))
    clearvars -except doralLvl smFish c lineCol paperData
    %axis tight
    xlim([-15 15])

end
save_gf(gcf,'Figure6Mist','type',{'svg','jpg'},'paper','benny21')

%% Figure 7A
clearvars -except paperData
smFish=dir('/home/labs/barkailab/felixj/Documents/MATLAB/projects/shio/*/*smFISH_38*f6.mat')
smFish=smFish([2])
load('dorsalGradient.mat','doralLvl');
close all
figure('Units','normalized','Position',[0.2547 0.0787 0.2234 0.7917],'Color',[1 1 1])
lineCol=brighten([brewermap(1,'Greens');brewermap(1,'Reds')],-0.5);;
paperData{6,1}=smFish
paperData{6,2}='Fig7'

for i=1:numel(smFish)
    out=load([smFish(i).folder '/' smFish(i).name])
    if ~isfield(out,'roiPts')
        out.roiPts=table([0;0;size(out.nucleiMask,2);size(out.nucleiMask,2)],[0;size(out.nucleiMask,1);size(out.nucleiMask,1);0],'VariableNames',{'X','Y'})
        nucDis=40;
    elseif numel(out.roiPts)==0
        out.roiPts=table([0;0;size(out.nucleiMask,2);size(out.nucleiMask,2)],[0;size(out.nucleiMask,1);size(out.nucleiMask,1);0],'VariableNames',{'X','Y'})
        nucDis=40;
    else
        nucDis = sqrt(polyarea(out.roiPts.X,out.roiPts.Y)/min(400,sum(inpolygon(out.nuclei.Centroid(:,1),out.nuclei.Centroid(:,2),out.roiPts.X,out.roiPts.Y))));
    end
    bins=-16:16
    subplot(3,numel(smFish),i)
    hold off
    subplot(3,numel(smFish),i+numel(smFish))
    hold off
    subplot(3,numel(smFish),i+2*numel(smFish))
    hold off
    for j=find(~cellfun('isempty',out.probes))
        if ~ismember('midDis',out.ts{j}.Properties.VariableNames)
            out.ts{j}.midDis =disToLines(out.ts{j}.Centroid(:,1),out.ts{j}.Centroid(:,2),out.midLine(:,1),out.midLine(:,2));
        end
        locs=out.ts{j};
        [~,binIdx]=ismember(round(locs.midDis./nucDis),bins);        
        intLocs=inpolygon(locs.Centroid(:,1),locs.Centroid(:,2),out.roiPts.X,out.roiPts.Y) & binIdx>0;     
        crInd=corr(locs.int(intLocs),interp1((-14:15).*nucDis,doralLvl,locs.midDis(intLocs)),'rows','pairwise');
        nLoc=accumarray(binIdx(intLocs),1,[numel(bins),1],@sum,0);
        crN=corr(nLoc(3:end-1),doralLvl,'rows','pairwise');
        intLoc=accumarray(binIdx(intLocs),locs.int(intLocs),[numel(bins),1],@median,0);
        subplot(3,numel(smFish),i)
        plot(bins,smoothdata(nLoc,'gaussian',3),'LineWidth',2,'DisplayName',sprintf('%s: %.2f',out.probes{j},crN),'Color',lineCol(j,:))
        hold on
        subplot(3,numel(smFish),i+numel(smFish))
        plot(bins,smoothdata(intLoc,'gaussian',3),'LineWidth',2,'DisplayName',sprintf('%s: %.2f',out.probes{j},crInd),'Color',lineCol(j,:))
        hold on
        axis tight
        subplot(3,numel(smFish),2*i+numel(smFish))
        tempPic=out.pic(:,:,j);
        picSize=size(tempPic);
        selPixel=inpolygon(repmat([1:picSize(2)],picSize(1),1),repmat([1:picSize(1)]',1,picSize(2)),out.roiPts.X,out.roiPts.Y);%&out.nucleiMask==0;
        mRNALvl=accumarray(round(out.disMat(selPixel)-min(out.disMat(selPixel),[],'all')+1),tempPic(selPixel),[],@median);
        hold on        
        crMrna=corr(mRNALvl,interp1((-14:15).*nucDis,doralLvl,min(out.disMat(selPixel))+[0.5:numel(mRNALvl)-.5])','rows','pairwise');
        plot((min(out.disMat(selPixel))+[0.5:numel(mRNALvl)-.5])/nucDis,movmean(mRNALvl,50),'LineWidth',2,'DisplayName',sprintf('%s: %.2f mRNA',out.probes{j},crMrna),'Color',lineCol(j,:))
        axis tight
    end
    subplot(3,numel(smFish),i)
    xlim(bins(max(1,min(quantile(binIdx,[0 1])+[-1 1],numel(bins)))))
    plot((-14:15),rescale(doralLvl,min(ylim),max(ylim)),'k--','DisplayName','DorsalLvl') 
   % legend()
    xlabel('Distance from midline')
    ylabel('# of TSS')
    title(sprintf('%s %.2f',strrep(regexp(smFish(i).name,'\d+_.*?_\d+','match','once'),'_',' ' ),nucDis))
    axis tight
    
    subplot(3,numel(smFish),i+1*numel(smFish))
    xlim(bins(max(1,min(quantile(binIdx,[0 1])+[-1 1],numel(bins)))))
    plot((-14:15),rescale(doralLvl,min(ylim),max(ylim)),'k--','DisplayName','DorsalLvl')      
    %legend()    
    xlabel('Distance from midline')
    ylabel('median Intensity')
    %axis tight
    
    subplot(3,numel(smFish),i+2*numel(smFish))
    plot((-14:15),rescale(doralLvl,min(ylim),max(ylim)),'k--','DisplayName','DorsalLvl')      
    %legend()    
    xlabel('Distance from midline')
    ylabel('median Cytoplamisc intensity (mRNA)')
    clearvars -except doralLvl smFish c lineCol
    axis tight
end
save_gf(gcf,'Figure7A_T48-twsit','type',{'svg','jpg'},'paper','benny21')


%% figure 7B
clearvars -except paperData
smFish=dir('/home/labs/barkailab/felixj/Documents/MATLAB/projects/shio/*/*smFISH_35*f6.mat')
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/shio/extData/dorsalGradient.mat')
close all
figure('Units','normalized','Position',[0.2734 0.0861 0.4245 0.7157],'Color',[1 1 1])
c=0;
intSmp=[3,1]
paperData{7,1}=smFish(intSmp,:);
paperData{7,2}='Fig7B';

for i=intSmp
    clear xLimLine
    c=c+1;
    matFile=[smFish(i).folder,'/' smFish(i).name];
    out=load(matFile);
    midPoint=mean(out.midLine);
    rotAngle=180-atand((out.midLine(1,1)-out.midLine(end,1))./(out.midLine(1,2)-out.midLine(end,2)));    
    midLineRot=nan(size(out.midLine));
    midLineRot(:,1)=[(out.midLine(:,1)-midPoint(1))*cosd(rotAngle) + (out.midLine(:,2)-midPoint(2))*sind(rotAngle)];
    midLineRot(:,2)=[-(out.midLine(:,1)-midPoint(1))*sind(rotAngle) + (out.midLine(:,2)-midPoint(2))*cosd(rotAngle)];      
    nucDis=sqrt(polyarea(out.roiPts.X,out.roiPts.Y)/sum(inpolygon(out.nuclei.Centroid(:,1),out.nuclei.Centroid(:,2),out.roiPts.X,out.roiPts.Y)));
    disBins=[-20.5:20.5]*nucDis;
    cc=0;
    for j=[3,1]
        cc=cc+1;
        subplot(3,numel(intSmp),(cc-1)*numel(intSmp)+c)
        hold off
        locs=out.ts{j};
        if isfield(out,'roiPts')
            intLocs=inpolygon(locs.Centroid(:,1),locs.Centroid(:,2),out.roiPts.X,out.roiPts.Y);
        else
            intLocs=true(size(locs,1),1);
        end        
        locs.Xnew=[(locs.Centroid(:,1)-midPoint(1))*cosd(rotAngle) + (locs.Centroid(:,2)-midPoint(2))*sind(rotAngle)];
        locs.Ynew=[-(locs.Centroid(:,1)-midPoint(1))*sind(rotAngle) + (locs.Centroid(:,2)-midPoint(2))*cosd(rotAngle)];
        scatter(locs.Xnew(intLocs)/nucDis,locs.Ynew(intLocs)/nucDis,[],locs.int(intLocs),'filled','LineWidth',0.5,'MarkerEdgeColor',[.4 .4 .4])
        axis tight
        caxis(quantile(locs.int(intLocs),[0.1 0.9]))
        if contains(out.probes{j},'3')
            %load('yellowish.mat','yellowish')
            %tssColor=yellowish;
            %lineCol=brighten(tssColor(end,:),-0.4);
            tssColor=brewermap(64,'Greens');
            tssColor=tssColor(1:50,:);
            lineCol=tssColor(35,:);
        else
            tssColor=brewermap(64,'OrRd');
            tssColor=tssColor(1:50,:)
            lineCol=tssColor(50,:);
        end
        colormap(gca,tssColor)
        yLim=ylim();
        xLim=xlim();
        hold on
        plot(midLineRot(:,1)/nucDis,midLineRot(:,2)/nucDis,'k-','LineWidth',2)
        xlim([-10 10]);
        ylim(yLim);
        title(sprintf('%s,%s,%d,%.2f',strrep(extractBefore(smFish(i).name,'f6.mat'),'_',' '),out.probes{j},sum(intLocs),nucDis))
        xlabel('Position relative to midline')
        yticks([])
        ylabel(colorbar(),'TS spot intensity')
        set(gca,'BoxStyle','full')
        
        subplot(3,numel(intSmp),2*numel(intSmp)+c)
        locs.binIdx=sum(locs.midDis>disBins,2);
        plot(movmean(disBins,2,'Endpoints','discard')/nucDis,accumarray(locs.binIdx(locs.binIdx>0 & locs.binIdx<42 & intLocs),1,[41 1]),'DisplayName',out.probes{j},'LineWidth',2,'Color',lineCol)
        hold on
        xLimLine{j}=quantile(locs.midDis(intLocs),[0 1])/nucDis
    end
    subplot(3,numel(intSmp),2*numel(intSmp)+c)
    plot([-14:15],rescale(doralLvl,min(ylim()),max(ylim)),'k--','DisplayName','dorsal gradient')
    out.nuclei.midDis=disToLines(out.nuclei.Centroid(:,1),out.nuclei.Centroid(:,2),out.midLine(:,1),out.midLine(:,2));
    out.nuclei.binIdx=sum(out.nuclei.midDis>disBins,2);
    intNucs=inpolygon(out.nuclei.Centroid(:,1),out.nuclei.Centroid(:,2),out.roiPts.X,out.roiPts.Y)
    %plot(movmean(disBins,2,'Endpoints','discard')/nucDis,accumarray(out.nuclei.binIdx(out.nuclei.binIdx>0 & out.nuclei.binIdx<42 & intNucs),1,[41 1]),'k-','DisplayName','nuclei','LineWidth',2)
    xLimLine=cat(1,xLimLine{:})
    xlim([min(xLimLine(:,1)) max(xLimLine(:,2))]+[-1 1])
    %title('TS number')
end
save_gf(gcf,'Figure7B_T48-53','type',{'svg','jpg'},'paper','benny21')


%% calculate POL2 numbers
clearvars -except paperData
smFish=[dir('./*/*.nd2');dir('./*mist-ind/*.czi')]
close all
paperData{8,1}=smFish([1,9],:);
paperData{8,2}='Pol2'

for i=[1,9]
    try
    out=load([smFish(i).folder '/' smFish(i).name(1:end-4) 'f6.mat'])
    figure('Units','normalized','Position',[[0,0,1,1]],'Color',[1 1 1]);
    c=0;
    for j=find(~cellfun('isempty',out.ts))
        c=c+1
        locs=out.ts{j};
        subplot(sum(~cellfun('isempty',out.ts)),2,c)
        if ismember('MaxIntensity',locs.Properties.VariableNames)
            if iscell(locs.MaxIntensity)
                locs.MaxIntensity=cat(1,locs.MaxIntensity{:});
            end
            dscatter(log(locs.Volume(locs.int>0)),log(locs.MaxIntensity(locs.int>0)))
            hold on
            xlabel('Spot intensity')
            ylabel('Spot maximum')
            p=robustfit(log(locs.int(locs.int>0)),log(locs.MaxIntensity(locs.int>0)));
            %plot(xlim,xlim*p(2)+p(1),'k-')
            plot(xlim,[1;1].*(log(median(locs.MaxIntensity))+1),'k-')
            tss=locs.int>0& log(locs.MaxIntensity)>=(log(median(locs.MaxIntensity))+1);
            mrna=locs.int>0& log(locs.MaxIntensity)<(log(median(locs.MaxIntensity))+1);
        else
            scatter(locs.Volume,locs.int,[],locs.nucDis,'o','filled')
            caxis([0 5])
            xlabel('spot volume')
            ylabel('spot intensity - background')
            tss=locs.int>0& locs.Volume>=300;
            mrna=locs.int>0& locs.Volume<300;
        end
        %title(sprintf('%s-%s',out.channels{j},out.probes{j}))
        c=c+1;
        subplot(sum(~cellfun('isempty',out.ts)),2,c)
        [yAll,bins]=histcounts(log(locs.int(tss|mrna)),50,'Normalization','probability');
        [yTss,~]=histcounts(log(locs.int(tss)),bins,'Normalization','probability');
        [yMrna,~]=histcounts(log(locs.int(mrna)),bins,'Normalization','probability');
        hold off
        plot(movmean(bins,2,'omitnan','Endpoints','discard'),yAll,'-','LineWidth',2,'DisplayName','All spots')
        hold on
        plot(movmean(bins,2,'omitnan','Endpoints','discard'),yTss,'-','LineWidth',2,'DisplayName',sprintf('Big %.2fK',median(locs.int(tss))/1000))
        plot(movmean(bins,2,'omitnan','Endpoints','discard'),yMrna,'-','LineWidth',2,'DisplayName',sprintf('small %.2fK',median(locs.int(mrna))/1000))
        legend()
        xlabel('log(spot intensity')
        ylabel('probability')
    end
    end
end
save_gf(gcf,'Pol2Calc_Mist','type',{'svg','jpg'},'paper','benny21')
%% save data for export
clearvars -except paperData
for i=1:size(paperData,1)
    paperData{i,2}=repmat(paperData(i,2),numel(paperData{i,1}),1)
end
allFiles=struct2table(cat(1,paperData{:,1}))
allFiles.Figure=cat(1,paperData{:,2});
allFiles=removevars(allFiles,{'bytes','date','datenum','isdir'});
for i=1:size(allFiles,1)
    if contains(allFiles.name{i},'.mat')
        allFiles.matPath{i}=[allFiles.folder{i} '/' allFiles.name{i}];
        picFile=dir([allFiles.folder{i} '/../*/' extractBefore(allFiles.name{i},'f6.'),'*.czi']);
        allFiles.picPath{i}=[picFile.folder '/' picFile.name];
    else
        allFiles.picPath{i}=[allFiles.folder{i} '/' allFiles.name{i}]
        matFile=dir([allFiles.folder{i} '/../*/' extractBefore(allFiles.name{i},'.'),'*f6.mat']);
        allFiles.matPath{i}=[matFile.folder '/' matFile.name];
    end
end
%% cp files to folder
copyCmd='cp %s %s'
uploadFolder='/home/labs/barkailab/felixj/Documents/MATLAB/projects/shio/forUpload/';
for i=1:size(allFiles,1)
    system(sprintf(copyCmd,allFiles.matPath{i},uploadFolder))    
    system(sprintf(copyCmd,allFiles.picPath{i},uploadFolder))    
end
allFiles.name=regexp(allFiles.picPath,'([A-Za-z0-9_]+.czi)|([A-Za-z0-9_]+.nd2)','match','once')
writetable(allFiles,[uploadFolder '/figures.xls'])
probeTable=readtable('probes.xlsx');
[~,idx]=ismember(regexp(allFiles.picPath,'([A-Za-z0-9_]+.czi)|([A-Za-z0-9_]+.nd2)','match','once'),regexp(probeTable.Var1,'([A-Za-z0-9_]+.czi)|([A-Za-z0-9_]+.nd2)','match','once'));
probeTable=probeTable(unique(idx),:)
writetable(probeTable,[uploadFolder '/probes.xls'])
