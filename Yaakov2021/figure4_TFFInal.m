function figure3_TF()
clear all
load('/home/labs/barkailab/LAB/checTS20190916.mat','medianNorm')
%intChecProfiles={'Msn2','Yap1','Yap1_H2O2','Yap1_dtail'};
intChecProfiles={'Msn2','Yap1_H2O2','Msn2_dZF','Yap1_dDBD_30','Msn2_Zfonly','Yap1_dtail'}%,'Yap2','Crz1','Nrg1','SKO1','Mig1','Mig2','Mig3','Mig4'};%,'Msn2_Zfonly','Msn2_dZF',,'Yap1_dtail'};
%intChecProfiles={'Rpn4','Reb1','Rgm1','Swi5','ACA1','CST6','Com2'}
sumProm=nan(6701,numel(intChecProfiles));
%intChecProfiles={'Msn2','Yap1_H2O2','Msn2_dZF','Yap1_dDBD_30'}
for i=1:numel(intChecProfiles)
    checProfiles(:,i)=medianNorm.(intChecProfiles{i});
end
GP=load('group_imp.mat');
checProfiles=[cell2mat(checProfiles') nan(size(checProfiles,2),GP.chr_len(17))]';
checProfiles=movmean(checProfiles,25,1,'omitnan');
clearvars -except checProfiles intChecProfiles;
intChecFactors={'Abf1_GZ','Reb1_GZ','Rap1_GZ','Cbf1_DK','Hsf1_DK'}
for i=1:numel(intChecFactors)
    temp=load('ChecProfile.mat',intChecFactors{i});
    checProfiles(:,end+1)=temp.(intChecFactors{i});
    intChecProfiles(end+1)=intChecFactors(i);
end
clearvars -except checProfiles intChecProfiles;
[data,smeta]=getNaamaData('figure',3);
%[data,smeta]=getNaamaData('figure',14);

smeta.tp=str2double(regexp(smeta.con,'\d+$','match','once'));
smeta.tp(isnan(smeta.tp)&contains(smeta.con,'async'),:)=-1;
smeta.tp(isnan(smeta.tp)&contains(smeta.con,'alpha'),:)=-5;
%[tcs,~,smeta.tcId]=unique(smeta(:,[12 7 11 14]),'stable')
%[tcs,~,smeta.tcId]=unique(smeta(:,[7 14]),'stable')
[tcs,~,smeta.tcId]=unique(smeta(:,[7 9 14]),'stable')

[~,idx]=sortrows(smeta,{'tag2','ab','tcId','tp'});
smeta=smeta(idx,:);
data=data(:,idx);
% hisNucs=[486002,486186;486281,486448;2275035,2275190;2275211,2275351];
% hisPos=zeros(size(data,1),1);
% hisMed=nan(size(hisNucs,1),size(data,2))
% for i=1:size(hisNucs,1)
%     hisPos(hisNucs(i,1):hisNucs(i,2),1)=i;
%     hisMed(i,:)=mean(data(hisPos==i,:),1)
% end
% hisMed=mean(log(hisMed),1)
% data(:,smeta.tcId==2)=data(:,smeta.tcId==2).*exp(-hisMed(smeta.tcId==2));
load('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/H2O2/ExpH2O2.mat')
load('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/H2O2/ExpH2O2F.mat')
ExpData(2)=H2O2_YPD;
ExpData(3)=H2O2_SD;
ExpData(1).Name='gilad'
ExpData(2).Name='YPD'
ExpData(3).Name='SC'
clearvars H2O2_SD H2O2_YPD
for i=1:numel(ExpData)
    selTp=ExpData(i).time>=0 & ExpData(i).time<=3;
    dynExp{i}=ExpData(i).expression-median(ExpData(i).expression(:,selTp),2,'omitnan');
end
maxChg=max(movmean(dynExp{1},3,2,'omitnan'),[],2);
%% show expression of target genes
GP=load('group_imp.mat');
metaProfile=meta_profile(checProfiles,'promoter',700,'useORF',false,'afterTSS',100,'scaled',false,'disCDS',false);
promLen=load('promoterLengthsGP.mat','promoterLengthsGP')
promLen=promLen.promoterLengthsGP;
sumProm2=nan(6701,size(checProfiles,2));
for i=find(~isnan(promLen)&~isnan(GP.gene_infoR64.stein_tss(:,1)))'
    promBases=[GP.gene_infoR64.stein_tss(i,2)-GP.gene_infoR64.dir(i)*promLen(i):GP.gene_infoR64.dir(i):GP.gene_infoR64.stein_tss(i,2)]+GP.chrIdx(GP.gene_infoR64.position(i,1));
    sumProm2(i,:)=mean(checProfiles(promBases,:),1)*700;
end
sumProm=permute(sum(metaProfile.cube(:,1:700,:),2),[1,3,2]);
targetTh=15000;
sumProm=sumProm;
minNumberOfTargets=100;
for i=1:size(sumProm,2)
    targetGenes{i}=find(sumProm(:,i)>targetTh);  
    if numel(targetGenes{i})<minNumberOfTargets
        [~,targetGenes{i}]=maxk(sumProm(:,i),minNumberOfTargets);
    end
    medTfExp(i,:)=median(dynExp{1}(targetGenes{i},:),1,'omitnan');
    [~,idx]=sort(sumProm(targetGenes{i},i),'descend');
    targetGenes{i}=targetGenes{i}(idx);
end

figure
for i=1:size(medTfExp,1)
    subplot(ceil(size(medTfExp,1)/2),2,i)
    plot(ExpData(1).time,medTfExp(i,:),'-','LineWidth',2)
    hold on
    title(intChecProfiles{i})
   axis('tight')
end

% figure
% subplot(2,3,1)
% imagesc(medTfExp(selTf,:)')
% xticks(1:numel(selTf))
% xticklabels(intChecProfiles(selTf))
% title('target gene expression')
% subplot(1,3,2)
% imagesc(dynExp{1}(cat(1,targetGenes{selTf}),:)',[-1 1]*1)
% subplot(1,3,3)
% plot(medTfExp(selTf,:)','LineWidth',2)

ws=1000;
clear peakPos peakGenes
for i=1:size(sumProm,2)      
    peakPos{i}=zeros(numel(targetGenes{i}),1)
    for g=1:numel(targetGenes{i})
        tss=GP.gene_infoR64.felixTss(targetGenes{i}(g),:);
        %tss=GP.gene_infoR64.stein_tss(targetGenes{i}(g),:);
        tss=GP.chrIdx(tss(1))+tss(2:3);
        if GP.gene_infoR64.dir(targetGenes{i}(g)) ==1
            [~,peakPos{i}(g)]=findpeaks(movmean(checProfiles(tss(1)-ws:tss(1),i),5),tss(1)-ws:tss(1),'Npeaks',1,'SortStr','descend');
        else
            [~,peakPos{i}(g)]=findpeaks(movmean(checProfiles(tss(1):tss(1)+ws,i),5),tss(1):tss(1)+ws,'Npeaks',1,'SortStr','descend'); 
        end
    end
    [peakPos{i},idx]=sort(peakPos{i});
    targetGenes{i}=targetGenes{i}(idx);
end
goodGene=cell(size(sumProm,2),1);
tssIdx=nan(6701,2);
tssIdx(~isnan(GP.gene_infoR64.felixTss(:,1)),:)=GP.chrIdx(GP.gene_infoR64.felixTss(~isnan(GP.gene_infoR64.felixTss(:,1)),1))+GP.gene_infoR64.felixTss(~isnan(GP.gene_infoR64.felixTss(:,1)),2:3)

%% determine motif direction
scGenome=load('SC_genome.mat')
scGenome=cat(2,scGenome.SC_genome.Sequence);
motifs={'AGGGG'}
for i=1
    figure;
    for p=1:numel(peakPos{i})
        motifOcc{1}=regexp(scGenome(peakPos{i}(p)+[-100:100]),sprintf('%s',motifs{i}));        
        motifOcc{2}=regexp(scGenome(peakPos{i}(p)+[-100:100]),sprintf('%s',seqrcomplement(motifs{i})));
        subplot(7,8,p)
        plot(checProfiles(peakPos{i}(p)+[-100:100],i))
        hold on
        scatter(motifOcc{1},zeros(size(motifOcc{1})),'filled')
        scatter(motifOcc{2},zeros(size(motifOcc{2})),'filled')
        if numel(motifOcc{1})>numel(motifOcc{2})
            motifDir{i}(p)=1;
        elseif numel(motifOcc{1})<numel(motifOcc{2})
            motifDir{i}(p)=-1;
        else
            motifDir{i}(p)=0;
        end
    end
end
for i=1:size(sumProm,2)      
    [uniPeaks,~,peakId]=unique(peakPos{i},'stable')
    goodGene{i}=zeros(numel(peakId),2)
    for p=1:max(peakId)
        selPeaks=find(peakId==p);
        nonIndGene{i}(selPeaks,1)=[repmat(numel(selPeaks)==1,numel(selPeaks),1) & maxChg(targetGenes{i}(selPeaks))<0.5 & corr(ExpData(1).expression(targetGenes{i}(selPeaks),:)',medTfExp(i,:)')<.3];
        crSel=corr(dynExp{1}(targetGenes{i}(selPeaks),:)',medTfExp(i,:)','rows','pairwise');
        if  numel(selPeaks)==1
             goodGene{i}(selPeaks,1)=true
        else
             goodGene{i}(selPeaks,1)=false
        end
        %goodGene{i}(selPeaks,1)= crSel==max(crSel);
        %goodGene{i}(selPeaks,2)=crSel;
    end
end
filtPeaks=cell(size(peakPos));
filtGenes=cell(size(peakPos));
for i=1:numel(peakPos)
%     if numel(motifDir)>=i
%         selGene=goodGene{i}(:,1)==1 & abs(motifDir{i})'>0;
%         filtDir{i}=motifDir{i}(selGene)
%     else
        selGene=goodGene{i}(:,1)==1
%    end
    filtPeaks{i}=peakPos{i}(selGene);
    filtGenes{i}=targetGenes{i}(selGene);
    
%     filtPeaks{i}=peakPos{i}(nonIndGene{i});
%     filtGenes{i}=targetGenes{i}(nonIndGene{i});
end

subplot(1,2,1)
hold off
tfX=1;
tfY=11;
hold off
scatter(sumProm(:,1),sumProm(:,2),[],relChg,'.','DisplayName','All genes')
caxis([0 1])
hold on;plot(xlim,[1 1].*targetTh,'k--','DisplayName','Singal Th');
plot([1 1].*targetTh,ylim(),'k--','DisplayName','Singal Th');
selTf=[tfX,tfY]
styleType={'o','s'}
for i=selTf
    scatter(sumProm(filtGenes{i},1),sumProm(filtGenes{i},2),[],relChg(filtGenes{i}),styleType{i},'filled','DisplayName',sprintf('%s-targets',intChecProfiles{i}))
    scatter(sumProm(targetGenes{i}(nonIndGene{i}),1),sumProm(targetGenes{i}(nonIndGene{i}),2),[],[1 0 0],styleType{i},'DisplayName',sprintf('%s-non Ind',intChecProfiles{tfY}))
end

targetGenes{i}(nonIndGene{i})
ylabel(colorbar(),'log2 expression change in H2O2')
xlabel(intChecProfiles{tfX})
ylabel(intChecProfiles{tfY})
tfEnv=500;
selPlot=[3,4,7,8];
tfColor=brewermap(2,'accent');
c=0;

for i=1:numel(selTf)    
    posi=acol(filtPeaks{selTf(i)}'+[-tfEnv:tfEnv]');
    for j=1:numel(selTf);
        c=c+1;        
        surTf=mean(reshape(checProfiles(posi,selTf(j)),2*tfEnv+1,[]),2);
        subplot(2,4,selPlot(c))
        hold off
        plot([-tfEnv:tfEnv],surTf,'Color',tfColor(j,:),'Linewidth',2)  
        xlim([-tfEnv tfEnv])
        ylim([0 250])
        xlabel(sprintf('Rel to %s BS',intChecProfiles{selTf(i)}))
        ylabel(sprintf('%s signal',intChecProfiles{selTf(j)}))
    end
end
save_gf(gcf,'SF3_TftargetSel')


[tcs,~,smeta.tcId]=unique(smeta(:,{'ab','tag2','gt'}));
selPlot=[11,12,5,6,9,10,3,4];
%tcCLim=[0,9;0,12;0 7;0 25];
%tcCLim=[0,9;0,9;0,9;0,9;0,9;0,9;0 30;0 30;0 30;0 30;0 30;0 30];
selTcs=[1,2,3,4];
cAxis{1}=[0 8]
cAxis{2}=[0 8]
cAxis{3}=[0 25]
cAxis{4}=[0 25]
selTf=[1,2,8,9]
intChecProfiles=strrep(intChecProfiles,'_',' ')
figure
c=0;
for t=selTcs%1:size(tcs,1)
    selSmp=find(smeta.tcId==t &~smeta.bad);
    for i=selTf
        peakSur=nan(numel(filtGenes{i}),2*tfEnv+1,numel(selSmp));
        for g=1:numel(filtGenes{i})
            if GP.gene_infoR64.dir(filtGenes{i}(g))==1
                posi=filtPeaks{i}(g)+[-tfEnv:tfEnv];
            else
                posi=filtPeaks{i}(g)+[tfEnv:-1:-tfEnv];
            end
            peakSur(g,:,:)=data(posi,selSmp);
        end
        medSur=squeeze(median(peakSur,1,'omitnan'));
        c=c+1;
        subplot(numel(selTcs),numel(selTf),c)
        %subplot(numel(selTf),size(tcs,1),(i-1)*size(tcs,1)+t)
        hold off
        imagesc(medSur')
        xticks([501]);        
        xticklabels(sprintf('%s (%d)',intChecProfiles{i},numel(filtPeaks{i})))
        yticks([])
        if i==selTf(1)
            ylabel(sprintf('\\fontsize{15} %s',[tcs.tag2{t} '-' tcs.ab{t} '-' tcs.gt{t}]))
            yticks(1:numel(selSmp))
            yticklabels(smeta.con(selSmp))
        elseif i==selTf(end)
            colorbar('east')
        end
        caxis(cAxis{t})
        if t==selTcs(1)
            title(intChecProfiles{i})
        end
    end
end

%%
hisEnv=[-700,300];
figure
for t=1:size(tcs,1)
    imageMat=nan(sum(smeta.tcId==t),range(hisEnv)+1,numel(GP.groups{23}{2}{45}));
    for g=1:numel(GP.groups{23}{2}{45})
        selPos=tssIdx(GP.groups{23}{2}{45}(g),1)+(hisEnv(1):hisEnv(2)).*GP.gene_infoR64.dir(GP.groups{23}{2}{45}(g));
        imageMat(:,:,g)=data(selPos,smeta.tcId==t)';        
    end
    subplot(2,2,t)
    imagesc(mean(imageMat,3,'omitnan'),'XData',[-700 300])
    colorbar()
    yticks(1:size(imageMat,1))
    yticklabels(smeta.tp(smeta.tcId==t))
    title(strjoin(tcs{t,:}))
    xlabel('hsitone loci')
    xticks([-300,0,300])
    ylabel('time after HS')
end

%% 
close all
for i=1:size(peakSur,1)
    if mod(i-1,12)==0
        figure
    end    
    subplot(3,4,mod(i-1,12)+1)
    imagesc(squeeze(peakSur(i,:,:))')
end

%%

expSel=1;
filtTfExp=nan(size(checProfiles,2),size(dynExp{expSel},2));
for i=1:3
    filtTfExp(i,:)=median(dynExp{expSel}(filtGenes{i},:),'omitnan');
end
subplot(3,3,1)
hold off
for i=1:numel(selTf)
    plot(ExpData(expSel).time,filtTfExp(selTf(i),:),'-','Color',tfColor(i,:),'LineWidth',2,'DisplayName',sprintf('%s targets',intChecProfiles{selTf(i)}))
    hold on
    axis tight
end
ylabel('median log(ExpChange)')
xlabel('time in H2O2')
% for t=1:size(tcs,1)
%     selSmp=find(smeta.tcId==t &~smeta.bad)    
%     for i=selTf        
%         clear peakSur
%         for g=1:numel(peakGenes{i})
%             if GP.gene_infoR64.dir(peakGenes{i}(g))==1
%                 posi=peakPos{i}(g)+[-300:300];
%             else
%                 posi=peakPos{i}(g)+[300:-1:-300];
%             end
%             peakSur(g,:,:)=data(posi,selSmp);
%         end
%         medSur{i}=squeeze(median(peakSur,1,'omitnan'))
%     end
%     subplot(2,2,t)
%     imagesc(cat(1,medSur{selTf})')
%     colorbar()
%     borders=cumsum(cellfun(@(x)size(x,1),medSur(selTf)));
%     hold
%     plot(borders.*[1;1]+.5,ylim'.*ones(1,numel(borders)),'k-');
%     xticks(movmean([0 borders],2,'Endpoints','discard'));
%     xticklabels(intChecProfiles(selTf))
%     title([tcs.tag2{t} '-' tcs.ab{t}])
%     yticks(1:numel(selSmp))
%     yticklabels(smeta.tp(selSmp))    
% end



suptitle(intChecProfiles{selTf})

[tcs,~,smeta.tcId]=unique(smeta(:,[14 7]),'stable');
selTf=3;
figure
for t=1:size(tcs,1)
    selSmp=find(smeta.tcId==t &~smeta.bad)
    clear peakSur
    for g=1:numel(peakGenes{selTf})
        if GP.gene_infoR64.dir(peakGenes{selTf}(g)) ==1
            posi=peakPos{selTf}(g)+[-300:300];
        else
            posi=peakPos{selTf}(g)+[300:-1:-300];
        end
        peakSur(g,:,:)=data(posi,selSmp);
    end
    subplot(2,2,t)
    imagesc(squeeze(median(peakSur,1,'omitnan'))','Xdata',[-300 300])
    colorbar()
    title([tcs.tag2{t} '-' tcs.ab{t}])
    yticks(1:numel(selSmp))
    yticklabels(smeta.tp(selSmp))
%     for s=1:numel(selSmp)
%         subplot(3,5,s)
%         imagesc(peakSur(:,:,s),'Xdata',[-300 300],[0 15])
%         title(smeta.name(selSmp(s)))
%     end
    
end
suptitle(intChecProfiles{selTf})
end