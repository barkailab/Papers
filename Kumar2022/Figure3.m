load('redData.mat')
%% add truncations specific data
sTable.bad(abs(sTable.end-tfs.length(sTable.tfid))>1 &~ismember(sTable.start,[0,2])&sTable.start>=0,:)=1; % ignore internal deletions
termName={'Nterm','Cterm'}
close all
c=0;
figure('Color',[1 1 1],'Position', [1 41 1920 962],'Renderer','painters')
sTable.addEffect(:)=nan;
sTable.addCut=repmat({[]},height(sTable),1);
for i=find(accumarray(sTable.tfid,1)>1)'
    for Cterm=[0,1]
        selSmp=find((sTable.tfid==i&sTable.C==Cterm&sTable.bad==0&sTable.lostBind==0)|(sTable.tfid==i&sTable.rmv==0&sTable.bad==0&sTable.lostBind==0));
        selSmp=selSmp(~ismember(sTable.type(selSmp),[3,5]));
        xPos=sTable.rmv(selSmp);
        xPos(xPos<0)=numel(tfs.disOrder{i})+xPos(xPos<0);    
        if numel(selSmp)>1
            if c==15
                saveas(gcf,['/home/labs/barkailab/felixj/Documents/MATLAB/projects/Divya/Fig3/IDR' num2str(i) '.jpg'])
                close gcf
                c=0;
                figure('Color',[1 1 1],'Position', [1 41 1920 962],'Renderer','painters')
            end
            c=c+1;
            pl=floor((c-1)/5)*10+mod((c-1),5)+1;
            
            crMat=corr(sumProm(promType<3,selSmp),'rows','pairwise');
            subplot(6,5,pl)
            xLim=[0 max(xPos)+50];
            plot(xPos,crMat(:,1),'Linewidth',2)
            xlim(xLim)
            ylabel('crWT')
            xticks([])
            title(sprintf('%s-%s',tfs.name{i},termName{Cterm+1}))
            
            yyaxis right
            bar(xPos,-[0;diff(crMat(:,1))],'FaceAlpha',.5)
            xlim(xLim)
            xticks([])
            ylabel('Effect')
            
            subplot(6,5,pl+5)
            cutFrag=zeros(numel(selSmp),numel(tfs.disOrder{i}));
            for s=1:numel(selSmp)
                if sTable.rmv(selSmp(s))>=0
                    selRes=sTable.start(selSmp(s)):sTable.end(selSmp(s));
                    selRes=selRes(selRes>0);
                    cutFrag(s,selRes)=1;
                else
                    selRes=abs(sTable.start(selSmp(s))):abs(sTable.end(selSmp(s)));
                    selRes=selRes(selRes>0);
                    cutFrag(s,:)=1;
                    cutFrag(s,selRes)=0;
                end
            end
            addCut=[zeros(1,numel(tfs.disOrder{i}));diff(cutFrag,1,1)];
            cutIdr=cell2mat(arrayfun(@(x)accumarray(sum(tfs.disOrder{i}(addCut(x,:)==1)'>=[0,0.25,0.5,0.75],2),1,[4,1]),1:numel(selSmp),'UniformOutput',false))
            bar(xPos,cutIdr','stacked')
            xlim(xLim)
            %selSmp=find((sTable.tfid==i&sTable.C==Cterm&sTable.bad==0)|(sTable.tfid==i&sTable.rmv==0&sTable.bad==0));
            set(gca,'ColorOrder',flipud(bone(4)))
            ylabel('# residues')
            xlabel('AA removed')
            sTable.addEffect(selSmp)=[0;-diff(crMat(:,1))];
            sTable.addCut(selSmp)=mat2cell(addCut,ones(height(addCut),1),width(addCut));
        end
    end
end
saveas(gcf,['/home/labs/barkailab/felixj/Documents/MATLAB/projects/Divya/Fig3/IDR' num2str(i) '.jpg'])
close gcf

sTable.meanIdr(~isnan(sTable.addEffect))=arrayfun(@(x)mean(tfs.disOrder{sTable.tfid(x)}(sTable.addCut{x}==1)),find(~isnan(sTable.addEffect)));
sTable.idrFrac(~isnan(sTable.addEffect))=arrayfun(@(x)mean(tfs.disOrder{sTable.tfid(x)}(sTable.addCut{x}==1)>0.5),find(~isnan(sTable.addEffect)));
sTable.addRmv(~cellfun('isempty',sTable.addCut))=cellfun(@(x)sum(x),sTable.addCut(~cellfun('isempty',sTable.addCut)))

%% Figure 5B
close all
[~,Ts]=ismember({'MOT3','SKN7','SWI5'},tfs.name)
cMapGen=flipud(brighten(brewermap(128,'RdYlBu'),0.5))
cMapBG=brewermap(128,'Greys')
cMapBG=brighten(cMapBG(5:69,:),.45);
cMap=brewermap(128,'BuPu')
cMap=cMap(10:end,:)

figure('Color',[1 1 1],'Position', [317 202 543 695],'Renderer','painters');
c=0
for t=Ts
    c=c+1;
    subplot(3,2,c)
    hold off
    selSmp=find(sTable.tfid==t&sTable.bad==0&sTable.lostBind==0&ismember(sTable.type,[1,2,4]))
    intSmp=selSmp(ismember(sTable.rmv(selSmp),[0 max(sTable.rmv(selSmp))]))
    selSmp=selSmp(sTable.C(selSmp)==sTable.C(intSmp(2)) | sTable.rmv(selSmp)==0);
    selTargets=find(any(((sumProm(:,intSmp)-mean(sumProm(promType<3,intSmp),'omitnan'))./std(sumProm(promType<3,intSmp)))>2.5,2) & promType<3);%promType<3%
    [~,idx]=sort(diff(log2(sumProm(selTargets,intSmp)),1,2));
    selTargets=selTargets(idx)
    imagesc(log2(sumProm(selTargets,selSmp)+700)-mean(log2(sumProm(selTargets,selSmp)+700),2),[-2 2])
    %xlabel('AA removed')
    xticks(1:numel(selSmp))
    xticklabels(sTable.rmv(selSmp))
    ylabel(sprintf('WT + DBD targets (%d)',numel(selTargets)))
    yticks([])
    %ylabel(colorbar(),'FC binding')
    title(tfs.name{t})
    colormap(gca,cMapGen)
         if c==5
ylabel(colorbar('Location','north'),'binding FC-log2')
     end
    c=c+1
    subplot(3,2,c)
    hold off
    ps = corr(sumProm(selTargets,selSmp),sumProm(selTargets,intSmp));
    xFit=fit(sTable.rmv(selSmp),ps(:,1),'linearinterp');
    yFit=fit(sTable.rmv(selSmp),ps(:,2),'linearinterp');
    missingPoint=unique(min(0:50:max(sTable.rmv(selSmp)+50),max(sTable.rmv(selSmp))));
    scatter(yFit(missingPoint),xFit(missingPoint),10,[.7 .7 .7],'filled')
    hold on
    scatter(ps(:,2),ps(:,1),100,sTable.rmv(selSmp)/tfs.nDLength(t)*100,'filled')
    scatter(ps(ismember(selSmp,intSmp),2),ps(ismember(selSmp,intSmp),1),100,'ko','Linewidth',2)
    %scatter(ps(sTable.C(selSmp)==1,1),ps(sTable.C(selSmp)==1,2),100,[.3 .3 .3],'Linewidth',2)
    caxis([0 100])
    %ylabel(colorbar(),'% nonNDBD truncated')
    xlabel('cr with DBD')
    ylabel('cr with WT')
    %title('Correlation over targets')
             bgIdx=round(rescale(-diag(corr(sumProm(promType<3,intSmp)),1),0.5,height(cMapBG)+0.4,'InputMin',-1,'InputMax',-0.2))
    bgIdx=round(rescale(diag(corr(sumProm(promType<3,intSmp)),1),0.5,height(cMapBG)+0.4,'InputMin',.5,'InputMax',1))

     %set(gca,'Color',cMapBG(bgIdx,:))
     xLim=xlim();
     yLim=ylim();
     %rectangle('Position',[xLim(1),yLim(1),diff(xLim),diff(yLim)],'FaceColor',cMapBG(bgIdx,:),'linestyle','none')
     % to add a background color accoridng to the wt-dbd correlation

     colormap(gca,cMap)
     if c==6
         ylabel(colorbar('Location','north'),'% nonDBD deleted')
     end
end
saveas(gcf,'/home/labs/barkailab/felixj/Documents/MATLAB/projects/Divya/svg/Fig3B-rev.svg')

%% calculate big jump for all terminal DBDs - processing
selTruncs=sTable.bad==0&sTable.lostBind==0&ismember(sTable.type,[2,4,5]);
tfs.C=accumarray(sTable.tfid(selTruncs),sTable.C(selTruncs),[height(tfs),1],[],0);
tfs.N=accumarray(sTable.tfid(selTruncs),sTable.C(selTruncs)==0,[height(tfs),1],[],0);
tfs.NC=accumarray(sTable.tfid(selTruncs),sTable.rmv(selTruncs)<0,[height(tfs),1],[],0);




Ts=find(tfs.nTrunc>1 &  tfs.dIdx>0)'
close all
figure('Color',[1 1 1],'Position', [1 41 1920 962],'Renderer','painters');
c=0;
bigJump=nan(height(tfs),1);
nTrunc=nan(height(tfs),1);
sTable.wtTargetCr(:)=NaN;
sTable.dbdTargetCr(:)=NaN;
for t=Ts
    c=c+1;
    subplot(6,10,c)
    selSmp=find(sTable.tfid==t&sTable.bad==0&sTable.lostBind==0&(sTable.type==1 | (abs(sTable.rmv)>0 & sTable.C == sTable.C(tfs.dIdx(t)))))
    rmvTemp=sTable.rmv(selSmp);
    rmvTemp(rmvTemp<0)=numel(tfs.disOrder{t})+rmvTemp(rmvTemp<0);
    intSmp=tfs{t,{'wtIdx','dIdx'}}
    if ~ismember(intSmp(2),selSmp)
        intSmp(2)=selSmp(end);
    end
    selTargets=find(any(((sumProm(:,intSmp)-mean(sumProm(promType<3,intSmp),'omitnan'))./std(sumProm(promType<3,intSmp)))>2.5,2) & promType<3);%promType<3%
    %     [pc,ps,pl,~,pe]=pca(log2(sumProm(selTargets,selSmp)+700)','NumComponents',2)
    ps = corr(sumProm(selTargets,selSmp),sumProm(selTargets,intSmp));
    
    sTable.wtTargetCr(selSmp)=ps(:,1);
    sTable.dbdTargetCr(selSmp)=ps(:,2);
    
    xFit=fit(rmvTemp,ps(:,1),'linearinterp');
    yFit=fit(rmvTemp,ps(:,2),'linearinterp');
    missingPoint=unique(min(0:50:max(rmvTemp+50),max(rmvTemp)));
    scatter(yFit(missingPoint),xFit(missingPoint),10,[.7 .7 .7],'filled')
    hold on
    scatter(ps(:,2),ps(:,1),100,rmvTemp/tfs.nDLength(t),'filled')
    scatter(ps(ismember(selSmp,intSmp),2),ps(ismember(selSmp,intSmp),1),100,'ko','Linewidth',2)
    %scatter(ps(sTable.C(selSmp)==1,1),ps(sTable.C(selSmp)==1,2),100,[.3 .3 .3],'Linewidth',2)
    caxis([0 1])
    %xlabel(sprintf('cr with %s',sTable.label{intSmp(1)}))
    %ylabel(sprintf('cr with %s',sTable.label{intSmp(2)}))
    %title(sprintf('%s (n:%d)',tfs.name{t},numel(selTargets)))
    %text(ps(:,1),ps(:,2),sTable.label(selSmp))
    xticks([])
    yticks([])
    stepDist=diag(squareform(pdist([xFit(missingPoint),yFit(missingPoint)])),1)
    bigJump(t)=max((stepDist./diff(missingPoint'))./median(stepDist./diff(missingPoint')))
    nTrunc(t)=numel(selSmp)-1;
    % define new dbd based on truncations 
    crDBD(t)=diag(corr(sumProm(promType<3,intSmp)),1);
    dbdIdx(t)=intSmp(2);
    
    bgIdx=round(rescale(-crDBD(t),0.5,height(cMapBG)+0.4,'InputMin',-1,'InputMax',-0.2))
    set(gca,'Color',cMapBG(bgIdx,:))
    colormap(gca,cMap)
    text(mean(xlim),max(ylim),sprintf('%s: %.2f',tfs.name{t},bigJump(t)),'HorizontalAlignment','center','VerticalAlignment','bottom','Fontweight','bold')
end
bigJump(nTrunc==1)=nan;

%% summary graph 3
figure('Color',[1 1 1],'Position', [636 243 610 423],'Renderer','painters');
truncTh=2;
clear xVal yVal
xVal=log2(accumarray(sTable.tfid(sTable.bad==0&sTable.lostBind==0),sTable.rmv(sTable.bad==0&sTable.lostBind==0),[],@max)./nTrunc);
yVal=(-log2(bigJump)/4)+1%rescale(-log2(bigJump),0,1,'InputMax',0,'InputMin',-4)
hold off
scatter(xVal,yVal,rescale(1-crDBD,30,300),nTrunc','filled')
xlabel('Average length of truncation (log2)')
ylabel('graduality (~median vs. max (log2))')
text(xVal(nTrunc>truncTh),yVal(nTrunc>truncTh),tfs.name(nTrunc>truncTh))
text(xVal(nTrunc<=truncTh),yVal(nTrunc<=truncTh),tfs.name(nTrunc<=truncTh),'COlor',[.5 .5 .5])
ylim([0 1])
colormap([.7 .7 .7;brewermap(128,'reds')])
caxis([truncTh 10])
ylabel(colorbar(),'# truncs')
saveas(gcf,'/home/labs/barkailab/felixj/Documents/MATLAB/projects/Divya/svg/Fig3F-summary.svg')
%% plot graphs - FIgure 3C

simpleTs=find(nTrunc>truncTh&min(tfs{:,{'N','C'}},[],2)<2 & ~ismember(tfs.name,{'MOT3','SKN7','SWI5','LYS14'}))'
[~,idx]=sortrows([nTrunc(simpleTs)<5 | crDBD(simpleTs)'>=0.6,bigJump(simpleTs)])
simpleTs=simpleTs(idx);
close all
figure('Color',[1 1 1],'Position', [1216 33 713 954],'Renderer','painters');
c=0;
for t=simpleTs
    c=c+1;
    subplot(6,4,c)
    selSmp=find(sTable.tfid==t&sTable.bad==0&sTable.lostBind==0&(sTable.type==1 | (sTable.rmv>0 & sTable.C == sTable.C(tfs.dIdx(t)))))   
    intSmp=tfs{t,{'wtIdx','dIdx'}}   
    if ~ismember(intSmp(2),selSmp)
        intSmp(2)=selSmp(end);
    end
    selTargets=find(any(((sumProm(:,intSmp)-mean(sumProm(promType<3,intSmp),'omitnan'))./std(sumProm(promType<3,intSmp)))>2.5,2) & promType<3);%promType<3%
    %     [pc,ps,pl,~,pe]=pca(log2(sumProm(selTargets,selSmp)+700)','NumComponents',2)
    ps = corr(sumProm(selTargets,selSmp),sumProm(selTargets,intSmp));
    xFit=fit(sTable.rmv(selSmp),ps(:,1),'linearinterp');
    yFit=fit(sTable.rmv(selSmp),ps(:,2),'linearinterp');
    missingPoint=unique(min(0:50:max(sTable.rmv(selSmp)+50),max(sTable.rmv(selSmp))));
    scatter(yFit(missingPoint),xFit(missingPoint),10,0.5*[1 1 1],'filled')
    hold on
    scatter(ps(:,2),ps(:,1),100,sTable.rmv(selSmp)/tfs.nDLength(t),'filled')
    scatter(ps(ismember(selSmp,intSmp),2),ps(ismember(selSmp,intSmp),1),100,'ko','Linewidth',2)
    %scatter(ps(sTable.C(selSmp)==1,1),ps(sTable.C(selSmp)==1,2),100,[.3 .3 .3],'Linewidth',2)
    caxis([0 1])
    %xlabel(sprintf('cr with %s',sTable.label{intSmp(1)}))
    %ylabel(sprintf('cr with %s',sTable.label{intSmp(2)}))
    %title(sprintf('%s (n:%d)',tfs.name{t},numel(selTargets)))
    %text(ps(:,1),ps(:,2),sTable.label(selSmp))
    xticks([])
    yticks([])
    bgIdx=round(rescale(-diag(corr(sumProm(promType<3,intSmp)),1),0.5,height(cMapBG)+0.4,'InputMin',-1,'InputMax',-0.2))
    bgIdx=round(rescale(diag(corr(sumProm(promType<3,intSmp)),1),0.5,height(cMapBG)+0.4,'InputMin',.5,'InputMax',1))

    set(gca,'Color',cMapBG(bgIdx,:))
    
     xLim=xlim();
     yLim=ylim();
     %rectangle('Position',[xLim(1),yLim(1),diff(xLim),diff(yLim)],'FaceColor',cMapBG(bgIdx,:),'linestyle','none')
    colormap(gca,cMap)
    %text(mean(xlim),max(ylim),sprintf('%s: %.2f',tfs.name{t},rescale(-log2(bigJump(t)),0,1,'InputMax',0,'InputMin',-4)),'HorizontalAlignment','center','VerticalAlignment','bottom','Fontweight','bold')
    if nTrunc(t)<5 | crDBD(t)'>=0.6
        text(mean(xlim),max(ylim),sprintf('%s',tfs.name{t}),'HorizontalAlignment','center','VerticalAlignment','bottom','Fontweight','bold','Color',.6*[1 1 1])
    else
        text(mean(xlim),max(ylim),sprintf('%s',tfs.name{t}),'HorizontalAlignment','center','VerticalAlignment','bottom','Fontweight','bold')
    end
end
set(gcf, 'InvertHardCopy', 'off');

saveas(gcf,'/home/labs/barkailab/felixj/Documents/MATLAB/projects/Divya/svg/Fig3D-smallOneSided-rev.svg')

%% double sided E
compTs=find(min(tfs{:,{'N','C'}},[],2)>=2)
[~,idx]=sortrows([sTable.rmv(tfs.dIdx(compTs))<0,bigJump(compTs)])
compTs=compTs(idx)
figure('Renderer','painters','Color',[1 1 1], 'Position', [649 638 1263 138],'InvertHardCopy', 'off');
c=0;
for t=compTs'
    c=c+1;
    selSmp=find(sTable.tfid==t&sTable.bad==0&sTable.lostBind==0&sTable.type~=3)
    intSmp=tfs{t,{'wtIdx','dIdx'}}
    
    rmvTemp=sTable.rmv(selSmp);
    rmvTemp(rmvTemp<0)=numel(tfs.disOrder{t})+rmvTemp(rmvTemp<0);
    selTargets=find(any(((sumProm(:,intSmp)-mean(sumProm(promType<3,intSmp),'omitnan'))./std(sumProm(promType<3,intSmp)))>2.5,2) & promType<3);%promType<3%
    
    %     [pc,ps,pl,~,pe]=pca(log2(sumProm(selTargets,selSmp)+700)','NumComponents',2)
    ps = corr(sumProm(selTargets,selSmp),sumProm(selTargets,intSmp));
    subplot(1,6,c)
    hold off
    if 1==1
        fitSmp=sTable.C(selSmp)==sTable.C(intSmp(2)) | sTable.type(selSmp)==1;
        xFit=fit(rmvTemp(fitSmp),ps(fitSmp,1),'linearinterp');
        yFit=fit(rmvTemp(fitSmp),ps(fitSmp,2),'linearinterp');
        missingPoint=unique(min(0:50:max(rmvTemp(fitSmp))+50,max(rmvTemp(fitSmp))));
        scatter(yFit(missingPoint),xFit(missingPoint),10,[.7 .7 .7],'filled')
        hold on
        scatter(ps(fitSmp,2),ps(fitSmp,1),100,'o','MarkerEdgeColor',[1 .3 .3],'Linewidth',1.5)
    elseif sTable.rmv(intSmp(2))<0
        scatter(ps(sTable.rmv(selSmp)<0,2),ps(sTable.rmv(selSmp)<0,1),100,'o','MarkerEdgeColor',[1 .3 .3],'Linewidth',1.5)
        hold on
        scatter(ps(sTable.C(selSmp)==1,2),ps(sTable.C(selSmp)==1,1),100,'o','MarkerEdgeColor',[.3 .3 .3],'Linewidth',1.5)
    end
    cVal=rmvTemp/tfs.nDLength(t);
    %cVal(cVal<0)=1-cVal(cVal<0);
    scatter(ps(:,2),ps(:,1),100,cVal,'filled')
    %text(mean(xlim),max(ylim),sprintf('%s: %.2f',tfs.name{t},rescale(-log2(bigJump(t)),0,1,'InputMax',0,'InputMin',-4)),'HorizontalAlignment','center','VerticalAlignment','bottom','Fontweight','bold')
    text(mean(xlim),max(ylim),sprintf('%s',tfs.name{t}),'HorizontalAlignment','center','VerticalAlignment','bottom','Fontweight','bold')

    %bgIdx=round(rescale(-diag(corr(sumProm(promType<3,intSmp)),1),0.5,height(cMapBG)+0.4,'InputMin',-1,'InputMax',-0.2))
    bgIdx=round(rescale(diag(corr(sumProm(promType<3,intSmp)),1),0.5,height(cMapBG)+0.4,'InputMin',.5,'InputMax',1))

    set(gca,'Color',cMapBG(bgIdx,:))
         xLim=xlim();
     yLim=ylim();
     %rectangle('Position',[xLim(1),yLim(1),diff(xLim),diff(yLim)],'FaceColor',cMapBG(bgIdx,:),'linestyle','none')
     
    colormap(gca,cMap)
    xlim(round(xlim,1))
    xticks(xlim)
       ylim(round(ylim,1))
    yticks(ylim)
    caxis([0 1])
end
saveas(gcf,'/home/labs/barkailab/felixj/Documents/MATLAB/projects/Divya/svg/Fig3G-doubleside.svg')
figure('Renderer','painters','Color',[1 1 1]);
colormap(gca,cMapBG)
caxis([.5 1])
ylabel(colorbar('southoutside'),'Correlation with TF')
saveas(gcf,'/home/labs/barkailab/felixj/Documents/MATLAB/projects/Divya/svg/Fig3-BGColor.svg')

%% hsitogram IDR vs. nonIDR effect H
figure('Renderer','painters','Color',[1 1 1])
idr=(sTable.meanIdr>=0.5)&(sTable.lostBind==0)&(sTable.bad==0);
notIdr=(sTable.meanIdr<0.5)&(sTable.lostBind==0)&(sTable.bad==0);

xVal=accumarray(sTable.tfid(idr),max(0,sTable.addEffect(idr)));
yVal=accumarray(sTable.tfid(notIdr),max(0,sTable.addEffect(notIdr)));
% scatter(xVal,yVal,100,cellfun(@(x)sum(x),tfs.nonDbd),'filled')
% ylabel('effect of structured regions')
% xlabel('effect of IDRs')
% ylabel(colorbar(),'nonDBD length')
% hold on
% scatter(xVal(zcs),yVal(zcs),100,'ko','linewidth',2)
% plot(xlim,xlim,'k--')
figure('Renderer','painters','Color',[1 1 1])
zcs=ismember(tfs.name,descs.TF(ismember(descs.FAMILY,'Zinc cluster')))
clear a
sel=(xVal>0|yVal>0)&~zcs;
[a(:,1),b]=histcounts(xVal(sel),-0.025:0.05:1.025,'Normalization','cumcount')
[a(:,2),b]=histcounts(yVal(sel),-0.025:0.05:1.025,'Normalization','cumcount')
hold off
truncName={'idr','structured'}
for i=1:width(a)
    plot(movmean(b,2,'Endpoints','discard'),sum(sel)-a(:,i),'Linewidth',2,'Displayname',truncName{i})
    hold on
end
xlabel('truncaation effect')
ylabel('#tfs')
legend()
saveas(gcf,'/home/labs/barkailab/felixj/Documents/MATLAB/projects/Divya/svg/Fig3H-idrEffect.svg')
%% Figure 3G
sTable.bad(abs(sTable.end-tfs.length(sTable.tfid))>1 &~ismember(sTable.start,[0,2])&sTable.start>=0,:)=0
sTable.type(abs(sTable.end-tfs.length(sTable.tfid))>1 &~ismember(sTable.start,[0,2])&sTable.start>=0,:)=5
selTruncs=sTable.bad==0&sTable.lostBind==0&ismember(sTable.type,[2,4,5]);
tfs.C=accumarray(sTable.tfid(selTruncs),sTable.C(selTruncs),[height(tfs),1],[],0);
tfs.N=accumarray(sTable.tfid(selTruncs),sTable.C(selTruncs)==0,[height(tfs),1],[],0);
tfs.NC=accumarray(sTable.tfid(selTruncs),sTable.rmv(selTruncs)<0,[height(tfs),1],[],0);

for i=1:height(sTable.type)
    fullSeq=ones(size(tfs.disOrder{sTable.tfid(i)}));
    fullSeq(max(sTable.start(i),1):sTable.end(i))=0
    sTable.seq{i}=fullSeq;
end

% calculate big jump for all terminal DBDs
Ts=find(tfs.nTrunc>1 &  tfs.dIdx>0)'
close all
figure('Color',[1 1 1],'Position', [1 41 1920 962],'Renderer','painters');

bigJump=nan(height(tfs),1);
nTrunc=nan(height(tfs),1);
sTable.wtTargetCr(:)=NaN;
sTable.dbdTargetCr(:)=NaN;
subplot(1,1,1)
c=0;
c=0;
truncEffect=nan(height(sTable),2);
delEffect=nan(height(sTable),2);

linEffect=nan(height(sTable),2);
relDBD=nan(height(sTable),1)
for s=find(sTable.type==5)'
    t=sTable.tfid(s)
    truncSmp=[find(sTable.tfid==t&sTable.bad==0&sTable.lostBind==0&ismember(sTable.type,[1,2,4]))]
    truncMat=[zeros(size(tfs.disOrder{t}));diff(1-cell2mat(sTable.seq(truncSmp)))];
    intSmp=tfs{t,{'wtIdx','dIdx'}}
    selTargets=find(any(((sumProm(:,intSmp)-mean(sumProm(promType<3,intSmp),'omitnan'))./std(sumProm(promType<3,intSmp)))>2.5,2) & promType<3);%promType<3%
    crs=corr(sumProm(selTargets,truncSmp),sumProm(selTargets,intSmp))
    truncEffect(s,:)=sum((sum(truncMat & (sTable.seq{s}==0),2)./sum(truncMat,2)).*[0,0;diff(crs)],'omitnan');
    delEffect(s,:)=diff(corr(sumProm(selTargets,[intSmp(1),s]),sumProm(selTargets,intSmp)));
    temp=corr(sumProm(selTargets,[intSmp(1),s]),sumProm(selTargets,intSmp(2)));
    relDBD(s)=diff(temp)./(1-temp(1));
    linSeq=ones(size(sTable.seq{s}));
    if mode(sTable.start(truncSmp))==2
        linSeq(2:sTable.rmv(s)+1)=0;
    else
        linSeq(end-sTable.rmv(s):end)=0;
    end    
    linEffect(s,:)=sum((sum(truncMat & (linSeq==0),2)./sum(truncMat,2)).*[0,0;diff(crs)],'omitnan');
end

set(gcf,'Color',[1 1 1],'Renderer','painters','Position',[350 210 1118 642])
subplot(1,1,1)
subplot(1,2,1)
hold off
%scatter(1+truncEffect(:,1),1+delEffect(:,1),200,-delEffect(:,1)+linEffect(:,1),'filled')
scatter(1+truncEffect(:,1),1+delEffect(:,1),200,sTable.tfid,'filled')
hold on
plot(xlim,xlim,'k--')
%caxis([0 .5])
text(1+truncEffect(:,1),1+delEffect(:,1),sTable.label)

xlabel('expected (trunc)')
ylabel('observed (del)')
title('wt corr')
axis square

%text(-truncEffect(:,1),-delEffect(:,1),sTable.label)
%colormap(gca,brighten(redblue(),0.5))
saveas(gcf,'Figure3Final.svg')
