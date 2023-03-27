%% 
load('redData.mat')

%% Figure 4C scatter DBD nonDBD
zSP=(sumProm-mean(sumProm(promType<3,:)))./std(sumProm(promType<3,:))
figure('Renderer','painters','Color',[1 1 1])
c=0;
for t=find(contains(tfs.name,{'UME6','DOT6'}))'
    c=c+1;
    subplot(2,2,c)
    selWt=find(sTable.type==1&sTable.tfid==t);
    selnD=find(sTable.type==3&sTable.tfid==t);
    selD=find(sTable.type==2&sTable.tfid==t);
    
    scatter(zSP(promType<3,selWt),zSP(promType<3,selnD),[],zSP(promType<3,selD),'filled','MarkerEdgeColor',[1 1 1].*.6)
    %colormap([.7 .7 .7; .3 .3 .3])
    caxis([0 8])
    xlabel(sTable.label{selWt})
    ylabel(sTable.label{selnD})
   
    colormap(gca,brewermap(128,'oranges'))
    xlim([-1 max(xlim)])
    ylim([-1 max(ylim)])
     text(xlim*[.95;.05],ylim*[.05;.95],num2str(corr(sumProm(promType<3,selWt),sumProm(promType<3,selnD)),'%.2f'))
     ylabel(colorbar(),'DBD signal')
end
figure('Renderer','painters','Color',[1 1 1])
c=0;
cMap=flipud(bone())
for t=find(contains(tfs.name,{'UME6','DOT6'}))'
    c=c+1;
    subplot(2,2,c)
    selWt=find(sTable.type==1&sTable.tfid==t);
    selnD=find(sTable.type==3&sTable.tfid==t);
    selD=find(sTable.type==2&sTable.tfid==t);
    crSP=corr(sumProm(promType<3,selWt),sumProm(promType<3,selnD));
    scatter(zSP(promType<3,selWt),zSP(promType<3,selnD),20,crSP.*ones(sum(promType<3),1),'filled')
    %colormap([.7 .7 .7; .3 .3 .3])
    caxis([0 1])
    xlabel(sTable.label{selWt})
    ylabel(sTable.label{selnD})
       xlim([-1.5 max(xlim)])
    ylim([-1.5 max(ylim)])
     text(xlim*[.95;.05],ylim*[.05;.95],num2str(crSP,'%.2f'))
     colormap(gca,cMap)

     c=c+1;
     subplot(2,2,c)
    crSP=corr(sumProm(promType<3,selWt),sumProm(promType<3,selD));   
     scatter(zSP(promType<3,selWt),zSP(promType<3,selD),20,crSP.*ones(sum(promType<3),1),'filled')
    %colormap([.7 .7 .7; .3 .3 .3])
    caxis([0 1])
    xlabel(sTable.label{selWt})
    ylabel(sTable.label{selD})
    xlim([-1.5 max(xlim)])
    ylim([-1.5 max(ylim)])
     text(xlim*[.95;.05],ylim*[.05;.95],num2str(crSP,'%.2f'))     
     colormap(gca,cMap)
end
saveas(gcf,'/home/labs/barkailab/felixj/Documents/MATLAB/projects/Divya/svg/Fig4C-Example.svg')

%% DBD vs. nonDBD correlations
selTfs=find(tfs.ndOk==1&tfs.dOk==1&~ismember(tfs.name,'FHL1'))
figure('Renderer','painters','Color',[1 1 1],'Position',[12 502 1701 349]);

subplot(1,1,1)
scatter(tfs.nDCr(selTfs),tfs.dCr(selTfs),100,sTable.repTelo(tfs.ndIdx(selTfs)),'filled','Linewidth',1,'MarkerEdgeColor',[0 0 0])
text(tfs.nDCr(selTfs),tfs.dCr(selTfs),tfs.name(selTfs))
xlabel('wt-nonDBD')
ylabel('wt-DBD')
ylabel(colorbar(),'nonDBD rep. correlation')
caxis([0.8 1])
colormap(gca,brewermap(128,'BuPU'))
xlim([0 1])
ylim([0 1])
axis square
saveas(gcf,'/home/labs/barkailab/felixj/Documents/MATLAB/projects/Divya/svg/Fig4C-Dbd-nonDBDCr.svg')

%% 4F,G correlation maps
%load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Divya/PughFeatures.mat','PughFeatures')
%stmProms=find(promType<3&~all(sumProm==0,2))%find(contains(PughFeatures.FeatureClassLevel1,'STM')&promType<3);%
tfTypes={'FL','DBD','nonDBD'}

figure('Renderer','painters','Color',[1 1 1])
clear promMat promCr cluN sizeClu
c=0
for t=[2,1,3]
    c=c+1;
    selSmp=find(sTable.type==t&sTable.lostBind==0&sTable.bad==0&~contains(sTable.full,'FHL1_nonDBD'))
    Y = pdist(sumProm(promType<3,selSmp)','correlation');
    Z = linkage(Y,'average');
    smpOrder = optimalleaforder(Z,Y);
    selSmp=selSmp(smpOrder);
    subplot(1,40,c*10+[-9 -1])
    %figure('Renderer','painters','Color',[1 1 1])
    crTemp{t}=corr(sumProm(promType<3,selSmp))
    imagesc(crTemp{t},[0 1])
    set(gca,'Ytick',1:numel(selSmp),'Yticklabel',sTable.TF(selSmp),'FontSize',5)
    title(sprintf('%s -%d',tfTypes{t},numel(selSmp)))
    colormap(gca,brewermap(128,'blues'))
    if t==3
        ylabel(colorbar('Location','west'),'correlation')
    end
    axis square
end

%figure('Renderer','painters','Color',[1 1 1])
subplot(2,4,4)
hold off
cMap=flipud(brewermap(5,'blues'));
for t=1:3
    plot(sort(1-squareform(1-crTemp{t})),2*[prod(size(crTemp{t})+[-1 0])/2:-1:1]./prod(size(crTemp{t})+[-1 0]),'Linewidth',2,'DisplayName',tfTypes{t},'Color',cMap(t+1,:))
    hold on
    xlim([0 1])
    xlabel('pairwise correlation')
    ylabel('cumultative fraction')
end
legend('Location','best')
%% Figure 4H
selND=find(ismember(sTable.TF,{'SKO1','CRZ1'})&sTable.type==3);
selTargets=find(any(((sumProm(:,selND)-mean(sumProm(promType<3,selND),'omitnan'))./std(sumProm(promType<3,selND)))>2.5,2) & promType<3);
[~,idx]=sort(diff(log2(sumProm(selTargets,selND)+700),1,2));
selTargets=selTargets(idx)
%tfOrder={'SOK2','SKO1','GIS1','LYS14','BAS1','MSN2','CRZ1','FKH2','SWI5','SWI4','TEC1'};
tfOrder={'SOK2','SKO1','GIS1','LYS14','BAS1','MSN2','CRZ1','SWI5','CHA4','SUM1','TEC1'};
selSmp=find(sTable.type==3&ismember(sTable.TF,tfOrder))
[~,idx]=ismember(tfOrder,sTable.TF(selSmp));
selSmp=selSmp(idx)
[~,top40]=maxk(sumProm(:,selSmp),40)
zSP=(sumProm-mean(sumProm(promType<3,:)))./std(sumProm(promType<3,:))
figure('Renderer','painters','Color',[1 1 1])
imagesc(zSP(top40(:),selSmp)',[0 10])
colormap(gca,brewermap(128,'Oranges'))
set(gca,'Ytick',1:numel(selSmp),'Yticklabel',sTable.TF(selSmp))
hold on
plot(xlim()',[1.5:1:numel(selSmp)].*[1;1],'-','Linewidth',1.5,'Color',0.7.*[1 1 1])
plot([40.5:40:numel(top40)].*[1;1],ylim','-','Linewidth',1.5,'Color',0.7.*[1 1 1])
ylabel(colorbar(),'cc')
saveas(gcf,'/home/labs/barkailab/felixj/Documents/MATLAB/projects/Divya/svg/Fig4F-nonDBDtargets.svg')
%% Figure 4I
intTfs={'SOK2','CRZ1'}
figure('Renderer','painters','Color',[1 1 1],'Position',[661 148 966 746]);
c=0;
for i=1
    for j=2
        c=c+1;
        subplot(2,2,c)
        hold off
        selNDx=find(ismember(sTable.TF,intTfs(i))&sTable.type==3)        
        selNDy=find(ismember(sTable.TF,intTfs(j))&sTable.type==3)
        
        selTFx=find(ismember(sTable.TF,intTfs(i))&sTable.type==1)        
        selTFy=find(ismember(sTable.TF,intTfs(j))&sTable.type==1)
        scatter(zSP(:,selNDx),zSP(:,selNDy),[],diff(log2(sumProm(:,[selTFx selTFy])+700),1,2),'filled','MarkerEdgeColor',[1 1 1].*.3)
        xlabel(sTable.TF{selNDx})
        ylabel(sTable.TF{selNDy})
        caxis([-3 3])
        ylabel(colorbar(),'wt ratios')
        colormap(gca,flipud(brighten(brewermap(128,'RdYlBu'),0.5)))
        hold on
        plot(xlim,xlim,'k--')
        axis square
        xlim([-1 25])
        ylim([-1 25])
    end
end
saveas(gcf,'/home/labs/barkailab/felixj/Documents/MATLAB/projects/Divya/svg/Fig4F-XCheck.svg')