function figure3_Ver3()
clear all
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/nucMeta015AddFried.mat')
[data,smeta]=getNaamaData('figure',14);
%[data,smeta]=getNaamaData('figure',7);

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
nucData(mean(nucData(:,strcmp(smeta.ab,'ha')),2)<0.01,:)=NaN;
nucData(mean(nucData(:,strcmp(smeta.ab,'ha')),2)>30,:)=NaN;
clear data

smeta.tp=str2double(regexp(smeta.con,'\d+$','match','once'));
smeta.tp(isnan(smeta.tp),:)=-1;
[~,idx]=sortrows(smeta,{'tag2','ab','tp'});
smeta=smeta(idx,:);
nucData=nucData(:,idx);
save('forFigure3add.mat')

clear all
load('forFigure3add.mat')
[tcs,~,smeta.tcId]=unique(smeta(:,[14 7]),'stable');
async=load('asyncMed.mat','medNuc','expTypes');
load('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/H2O2/ExpH2O2.mat')
load('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/H2O2/ExpH2O2F.mat')
ExpData(2,1)=H2O2_SD;
ExpData(3,1)=H2O2_YPD;;
ExpData(1).name='Exp';
ExpData(2).name='SD';
ExpData(3).name='YPD';
clearvars H2O2_SD H2O2_YPD

for i=1:size(ExpData,1)
    dynExp{i}=ExpData(i).expression-mean(ExpData(i).expression(:,:),2,'omitnan');
end
%% calculate dynamic Nucs

regev=load ('/home/labs/barkailab/LAB/data/DataExternal/EXPRESSION/regev_pnas_2008_rnaseq/data');
geneLvl=median(regev.regev_rnaseq,2,'omitnan');
clearvars regev i

%% find & cluster variable genes
scatter(mean(ExpData(1).expression(:,:),2,'omitnan'),std(movmean(ExpData(1).expression(:,:),2,2,'omitnan','Endpoints','discard'),[],2,'omitnan'),'.')
dynTh=getline();
%dynGenes=inpolygon(mean(ExpData(1).expression,2,'omitnan'),std(movmean(ExpData(1).expression,2,2,'omitnan','Endpoints','discard'),[],2,'omitnan'),dynTh(:,1),dynTh(:,2))
%hold on
%scatter(mean(ExpData(1).expression(dynGenes,:),2,'omitnan'),std(movmean(ExpData(1).expression(dynGenes,:),2,2,'omitnan','Endpoints','discard'),[],2,'omitnan'),'.')
expClu=nan(6701,1);
nClu=2;
expClu(dynGenes)=kmeans(movmean(dynExp{1}(dynGenes,:),2,2,'omitnan','Endpoints','discard'),nClu,'MaxIter',500,'Replicates',100,'Options',statset('UseParallel',true));
[~,expOrder]=sort(accumarray(expClu(expClu>0),median(dynExp{1}(expClu>0,11:14),2,'omitnan'),[],@(x)median(x,'omitnan')),'descend');
expClu(expClu>0)=changem(expClu(expClu>0),1:max(expClu),expOrder); %% induced genes should be expClu==1 now.
%% fitting of H2A myc decline for TEV
GP=load('group_imp.mat');
nClu=2
subRed=nan(6701,1);
subRed(expClu==2)=kmeans(movmean(dynExp{1}(expClu==2,:),2,2,'omitnan','Endpoints','discard'),nClu,'MaxIter',500,'Replicates',100,'Options',statset('UseParallel',true));;
groupColor=lines(4);
geneGroups={find(expClu==2 & subRed==1 & relChg<-1),find(expClu==2 & subRed==2 & relChg<-1),find(expClu==2 & ismember([1:6701]',GP.groups{7}{2}{65}) & relChg<-1),find(expClu==2 & ismember([1:6701]',GP.groups{7}{2}{66}) & relChg<-1)};
groupNames={'red1','red2','ribo-red','ribi-red'}
figure
for g=1:numel(geneGroups)
    medExp(g,:)=median(dynExp{1}(geneGroups{g},:),'omitnan');
    plot(ExpData(1).time,medExp(g,:),'Linewidth',2,'Color',groupColor(g,:),'DisplayName',groupNames{g})
    hold on
end
axis tight
xlabel('time')
ylabel('log2 mRNA change')
figure;
selSmp=smeta.tcId==2%;&smeta.tp>=-1 &smeta.tp<=28;
decayFunc=@(t,M0,M1,k) M1+((M0-M1).*2.^(-t./k));
fitEndPt=5;
for g=1:numel(geneGroups)
    selNucs=nucMeta.order>1 & ismember(nucMeta.gene,geneGroups{g});
    subplot(2,4,g)
    hold off
    yVal=nucData(selNucs,selSmp); 
    haVal=nucData(selNucs,smeta.tcId==1);
    errorbar(smeta.tp(selSmp),mean(yVal,'omitnan'),std(yVal,'omitnan')./sqrt(sum(~isnan(yVal))),'Color',groupColor(g,:),'DisplayName','Myc')    
    hold on
    fitFunc=@(p)sum((decayFunc([0:1:(fitEndPt-2)].*5,mean(mean(yVal(:,1:2))),0.95*mean(yVal(:,fitEndPt)),p(1))-mean(yVal(:,2:fitEndPt),'omitnan')).^2);
    p=fminunc(fitFunc,[5]);
    plot([0:1:(fitEndPt-2)].*5+3,decayFunc([0:1:(fitEndPt-2)].*5,mean(mean(yVal(:,1:2))),0.95.*mean(yVal(:,5)),p(1)),'--','LineWidth',2,'DisplayName','Fit')    
    plot(smeta.tp(smeta.tcId==1),mean(haVal),'k-','DisplayName','HA');
    ylabel('abs. myc')
    xlabel('time in H2O2')
    title(sprintf('%s-t_{1/2}:%.2f',groupNames{g},p))
    subplot(2,4,g+4)
    hold off
    yVal=nucData(selNucs,selSmp)./mean(nucData(selNucs,selSmp),2,'omitnan');
    haVal=nucData(selNucs,smeta.tcId==1)./mean(nucData(selNucs,smeta.tcId==1),2);
    errorbar(smeta.tp(selSmp),mean(yVal,'omitnan'),std(yVal,'omitnan')./sqrt(sum(~isnan(yVal))),'Color',groupColor(g,:),'DisplayName','myc')    
    hold on

    fitFunc=@(p)sum((decayFunc([0:1:(fitEndPt-2)].*5,mean(mean(yVal(:,1:2))),0.95*mean(yVal(:,fitEndPt)),p(1))-mean(yVal(:,2:fitEndPt),'omitnan')).^2);
    p=fminunc(fitFunc,[5]);
    plot([0:1:(fitEndPt-2)].*5+3,decayFunc([0:1:(fitEndPt-2)].*5,mean(mean(yVal(:,1:2))),0.95.*mean(yVal(:,5)),p(1)),'--','LineWidth',2,'DisplayName','Fit')    
    
    plot(smeta.tp(smeta.tcId==1),mean(haVal),'k-','DisplayName','HA');
    ylabel('normalized myc')
    xlabel('time in H2O2')
    title(sprintf('%s-t_{1/2}:%.2f',groupNames{g},p))
end
plot(mean(nucData(selNucs,smeta.tcId==2),'omitnan'))
save_gf(gcf,'Fig3_fit')
%% calculate nucleosome dynamics, i.e. change
dynNuc=nan(size(nucData));
for i=1:max(smeta.tcId)
    dynNuc(:,smeta.tcId==i&smeta.bad==0)=log(nucData(:,smeta.tcId==i&smeta.bad==0)+1)-log(median(nucData(:,smeta.tcId==i&smeta.bad==0)+1,2,'omitnan'));%log(medNuc(:,asynCon(i))+1);;    
    %dynNuc(:,smeta.tcId==i&smeta.bad==0)=(nucData(:,smeta.tcId==i&smeta.bad==0)+1)-(median(nucData(:,smeta.tcId==i&smeta.bad==0&smeta.tp<=0)+1,2,'omitnan'));
    tempMat=nucData(:,smeta.tcId==i&smeta.bad==0);
    rsMat=nan(size(tempMat));
    for j=find(~all(isnan(tempMat),2))'
        scaleI=[mean(maxk(tempMat(j,1:10),2)),mean(mink(tempMat(j,1:10),2))];
        rsMat(j,:)=rescale(tempMat(j,:),0,1,'InputMax',scaleI(1),'InputMin',scaleI(2));
        %rsMat(j,rsMat(j,:)<=0)=NaN;
    end
    rsNuc(:,smeta.tcId==i&smeta.bad==0)=rsMat;
end
%% fit H2A decline with rescaling (not necessary)
fitSmp= smeta.tcId==4 & smeta.tp<45  & smeta.tp>23;
crVal=corr(dynNuc(:,fitSmp)',mean(dynNuc(selNucs,fitSmp))');
selNucs=find(ismember(nucMeta.gene,find(expClu==2 & relChg<-1)) & nucMeta.order>0 & crVal>0.5 & range(nucData(:,fitSmp),2)>3);
for i=1:numel(selNucs)
    tpI=smeta.tp(fitSmp);
    mycI=nucData(selNucs(i),fitSmp);
    limitsI=[mycI(1) mean(mink(nucData(selNucs(i),smeta.tcId==4 & smeta.tp>23),2))];    
    scoreFunc=@(p) sum((sigfunc(tpI',p(1),28,limitsI) -mycI).^2)
    p(i)=fminunc(scoreFunc,[35]);
end
hist(p(1:361),100)
sigfunc(smeta.tp(fitSmp)',10,0,[11,7])
tempMat=nucData(selNucs,smeta.tcId==2);
rsMat=nan(size(tempMat));
clear f1 f2 f s1 s2 sAll fAll
sigfunc=@(x,turnPt,startPt,limits)limits(1)+diff(limits)./(1+exp(-(x-turnPt).*4./(turnPt-startPt)));
for i=1:sum(selNucs)
    scaleI=[mean(maxk(tempMat(i,1:10),1)),mean(mink(tempMat(i,1:10),1))];
    rsMat(i,:)=rescale(tempMat(i,:),0,1,'InputMax',scaleI(1),'InputMin',scaleI(2));
    rsMat(i,rsMat(i,:)<=0)=NaN;
    if sum(isnan(rsMat(i,:)))<2
        [fAll(i,:),sR]=robustfit([0,5,10,15],log2(rsMat(i,2:5)));
        sAll(i).normr=sR.robust_s;
    else
        sAll(i).normr=NaN;     
        fAll(i,:)=[NaN NaN];
    end
    if all(~isnan(rsMat(i,2:4)))
        [f1(i,:),s1(i)]=polyfit([0,5,10],log2(rsMat(i,2:4)),1); 
    else
        si.normr=NaN;
        si.df=1;
        si.R=[NaN,NaN;NaN,NaN];
        f1(i,:)=[NaN NaN];
        s1(i)=si;
    end
    if all(~isnan(rsMat(i,3:5)))
        [f2(i,:),s2(i)]=polyfit([0,5,10],log2(rsMat(i,3:5)),1); 
    else        
        si.normr=NaN;
        si.df=1;
        si.R=[NaN,NaN;NaN,NaN];
        f2(i,:)=[NaN NaN];
        s2(i)=si;
    end
end
s1=struct2table(s1);
s2=struct2table(s2);
sAll=struct2table(sAll);
fAll=fAll(:,[2,1]);
figure;
subplot(2,2,1)
imagesc(rsMat)
ylabel('all ORF nucleosme of Rep. genes')
xlabel('time')
xticks(1:sum(smeta.tcId==2))
xticklabels(smeta.tp(smeta.tcId==2))
ylabel(colorbar(),'rescaled H2A myc')
subplot(2,2,3)
plot(smeta.tp(smeta.tcId==2),mean(rsMat,1,'omitnan'),'Linewidth',2)
axis('tight')
ylabel('mean rescaled')
xlabel('time')
subplot(2,2,2)
hold off
varNames={'1','2','All'}
histBins=[-.5:0.02:.5];
for i=1:numel(varNames)
    sData=eval(sprintf('s%s',varNames{i}));
    fData=eval(sprintf('f%s',varNames{i}));
    y=histcounts(fData(:,1),histBins,'Normalization','pdf');
    plot(movmean(histBins,2,'Endpoints','discard'),y,'-','DisplayName',sprintf('%s - all Nucs:%.2f',varNames{i},mean(fData(:,1),'omitnan')),'Linewidth',2)
    hold on
    goodNucs= range(nucData(selNucs,find(smeta.tcId==2,6)),2)>5;
    y=histcounts(fData(goodNucs,1),histBins,'Normalization','pdf');
    plot(movmean(histBins,2,'Endpoints','discard'),y,'-','DisplayName',sprintf('%s - sel Nucs:%.2f',varNames{i},mean(fData(goodNucs,1),'omitnan')),'Linewidth',2)
end
xlabel('fit 1/T_{0.5}')

for i=1:sum(selNucs)
    hold off
    scatter(smeta.tp(smeta.tcId==2),log2(rsMat(i,:)),'o','filled')
    p=robustfit([0,5,10,15],log2(rsMat(i,2:5)))
    hold on
    plot([0,5,10,15]+3,(3+[0,5,10,15])*p(2)+p(1),'--')
    pause
end
%% all ind. nucleosomes
intNuc=2;
relChg=mean(ExpData(1).expression(:,9:11),2,'omitnan')-mean(ExpData(1).expression(:,2:3),2,'omitnan');
geneGroups={find(expClu==1 & relChg>1 &ismember([1:6701]',nucMeta.gene(nucMeta.order==intNuc))),find(expClu==2 & relChg<-1 & ismember([1:6701]',nucMeta.gene(nucMeta.order==intNuc)))}
groupNames={'ind','rep'};
for g=1:numel(geneGroups)
    [~,idx]=sort(relChg(geneGroups{g}),'descend');
    geneGroups{g}=geneGroups{g}(idx);
end
figure
imagesc(dynExp{1}(cat(1,geneGroups{:}),:)',[-2 2])


title('RNA')
colorbar()
yticks(1:2:size(dynExp{1},2))
yticklabels(ExpData(1).time(1:2:size(dynExp{1},2)))
hold on
linePos=[0 cumsum(cellfun('prodofsize',geneGroups))]+.5;
plot(linePos.*[1;1],ylim'.*[1 1 1],'k--','Linewidth',2)
xticks(movmean(linePos,2,'Endpoints','discard'))
xticklabels(groupNames)
figure;
selNucs=cell(numel(geneGroups),1);
for g=1:numel(geneGroups)
    temp=find(ismember(nucMeta.gene,geneGroups{g}) & nucMeta.order==intNuc);
    [~,idx]=sort(relChg(nucMeta.gene(temp)),'descend');
    selNucs{g}=temp(idx);
end
linePos=[0 cumsum(cellfun('prodofsize',selNucs))']+.5;
intTc=[2,1,4,3];
for t=1:numel(intTc)
    selSmp=smeta.tcId==intTc(t) & smeta.bad==0;
    subplot(1,numel(intTc),t)
    hold off
    imagesc(dynNuc(cat(1,selNucs{:}),selSmp)',[-.25 .25])
    colorbar()
    title(sprintf('%+d nuc, %s',intNuc,strjoin(tcs{intTc(t),:})))
    yticks(1:sum(selSmp));
    yticklabels(smeta.tp(selSmp))
    hold on
    plot(linePos.*[1;1],ylim'.*[1 1 1],'k--','Linewidth',2)
    xticks(movmean(linePos,2,'Endpoints','discard'))
    xticklabels(groupNames)
end
save_gf(gcf,sprintf('Fig3_%dnucIndivGenes',intNuc))

%% summary curves recreate 
intNuc=[1,2]
geneGroups={
            find(expClu==1 & relChg>1 &ismember([1:6701]',nucMeta.gene(ismember(nucMeta.order,intNuc)))),
            find(expClu==2 & relChg<-1 & geneLvl<10 &ismember([1:6701]',nucMeta.gene(ismember(nucMeta.order,intNuc)))) 
            find(expClu==2 & relChg<-1 & geneLvl>10 &ismember([1:6701]',nucMeta.gene(ismember(nucMeta.order,intNuc))))
            };        
groupNames={'ind','rep','high rep'}

intNuc=[-2,-1,1:3]
geneGroups={
            find(expClu==1 & relChg>1 &ismember([1:6701]',nucMeta.gene(ismember(nucMeta.order,intNuc)))),
            find(expClu==2 & relChg<-1 &ismember([1:6701]',nucMeta.gene(ismember(nucMeta.order,intNuc)))) 
            };        
groupNames={'ind','rep','high rep'}


for i=1:size(ExpData,1)
    dynExp{i}=ExpData(i).expression-mean(ExpData(i).expression(:,2:3),2,'omitnan');
end
dynNuc=nan(size(nucData));
for i=1:max(smeta.tcId)
    dynNuc(:,smeta.tcId==i&smeta.bad==0)=log(nucData(:,smeta.tcId==i&smeta.bad==0)+1)-log(median(nucData(:,smeta.tcId==i&smeta.bad==0 & smeta.tp<=3)+1,2,'omitnan'));%log(medNuc(:,asynCon(i))+1);;    
    %dynNuc(:,smeta.tcId==i&smeta.bad==0)=log(nucData(:,smeta.tcId==i&smeta.bad==0)+1)-log(median(nucData(:,smeta.tcId==i&smeta.bad==0&smeta.tp<=0)+1,2,'omitnan'));
end
intTc=[1,2;3,4]
close all
lineColors=flipud(brewermap(numel(geneGroups),'Set1'));
for n=1:numel(intNuc)
    %figure
    selNucs=cell(numel(geneGroups),1);
    for g=1:numel(geneGroups)
        selNucs{g}=find(nucMeta.order==intNuc(n) & ismember(nucMeta.gene,geneGroups{g}));
        medExp(g,:)=median(dynExp{1}(nucMeta.gene(selNucs{g}),:),'omitnan');
    end
    for t=1:size(intTc,1)
%         subplot(2,size(intTc,1),t)
%         for g=1:numel(geneGroups)
%             plot(ExpData(1).time,medExp(g,:),'Linewidth',2,'Color',lineColors(g,:))
%             hold on
%         end
%         axis tight
%         ylabel('RNA')
%         hold on
%         plot(xlim,[0 0],'k--')
%       
%         hold off
        subplot(size(intTc,1),numel(intNuc),(t-1).*numel(intNuc)+n)  
        hold off
        for tt=1:size(intTc,2)
            selSmp=smeta.tcId==intTc(t,tt) & smeta.bad==0;
            for g=1:numel(geneGroups)
                lineName=sprintf('%s %s %s',tcs.tag2{intTc(t,tt)},tcs.ab{intTc(t,tt)},groupNames{g});
                %subplot(numel(intNuc),size(intTc,1),(n-1).*size(intTc,1)+t)
                plot(smeta.tp(selSmp),mean(dynNuc(selNucs{g},selSmp),'omitnan')','Color',brighten(lineColors(g,:),-1.5+tt),'Linewidth',2,'DisplayName',lineName)
                hold on
            end
        end        
        %subplot(2,size(intTc,1),t+size(intTc,1))
        axis tight
        legend()
        xlabel('time')
        ylabel('log change')
        plot(xlim,[0 0],'k--')
        title(sprintf('%+d nuc',intNuc(n)))
        
        %subplot(2,size(intTc,1),t)
        %xlim(quantile(smeta.tp(selSmp),[0 1]))        
    end
    %save_gf(gcf,sprintf('Fig3_%dnucSumDynamics',intNuc(n)),1)
end
save_gf(gcf,sprintf('FigS3_AllnucAllTCSumDynamics',intNuc(n)))
%% correlation summary
relChg2=mean(ExpData(1).expression(:,9:11),2,'omitnan')-median(ExpData(1).expression(:,:),2,'omitnan');
close all
intNuc=[-2,-1,1,2]
for n=1:numel(intNuc)
    selNucs=ismember(nucMeta.gene,find(dynGenes & abs(relChg)>1)) & nucMeta.order==intNuc(n);
    crMat=corr(dynNuc(selNucs,:),relChg(nucMeta.gene(selNucs)),'rows','pairwise');
    subplot(2,2,n)
    for t=1:size(tcs)
        selSmp=smeta.tcId==t & smeta.bad==0;
        plot(movmean(smeta.tp(selSmp),2),movmean(crMat(selSmp),2),'-','LineWidth',2,'DisplayName',strjoin(tcs{t,:}))
        hold on
    end
    axis('tight')
    title(intNuc(n))
    plot(xlim(),[0 0 ],'k--')
    xlabel('time')
    ylabel(sprintf('correlation with %s','relChg'))
end
save_gf(gcf,'Fig3_corrSumSmooth');
%% lines curves
dynNuc=nan(size(nucData));
for i=1:max(smeta.tcId)
    dynNuc(:,smeta.tcId==i&smeta.bad==0)=log(nucData(:,smeta.tcId==i&smeta.bad==0)+1)-log(median(nucData(:,smeta.tcId==i&smeta.bad==0 & smeta.tp<=3)+1,2,'omitnan'));%log(medNuc(:,asynCon(i))+1);;    
    %dynNuc(:,smeta.tcId==i&smeta.bad==0)=log(nucData(:,smeta.tcId==i&smeta.bad==0)+1)-log(median(nucData(:,smeta.tcId==i&smeta.bad==0&smeta.tp<=0)+1,2,'omitnan'));
end
absDir=geneLvl.*sign(relChg)-sign(relChg)*3;
absDir(~dynGenes|abs(relChg)<1)=NaN;

selGenes=~ismember(1:6701,GP.groups{23}{2}{45})'%expClu==1 & relChg>1;
%intPar={'relChg','absDir'}
%movPar=[0.1;.5]
%parLim=[-2 2;-9 9]
intPar={'relChg'}
intNuc={[2 3]}
close all
figure
smpCol={[find(ismember(smeta.tcId,[1,2]) & ismember(smeta.tp,[13]))'];
    [find(ismember(smeta.tcId,[1,2]) & ismember(smeta.tp,[18]))'];
    [find(ismember(smeta.tcId,[1,2]) & ismember(smeta.tp,[23]))'];
    [find(ismember(smeta.tcId,[3,4]) & ismember(smeta.tp,[18]))'];
    [find(ismember(smeta.tcId,[3,4]) & ismember(smeta.tp,[23]))']}
%intSmp=[find(ismember(smeta.tcId,[1,2]) & smeta.tp==23)']
%lineColor=lines(numel(intSmp))
movPar=[0.2,1]
parLim=[-2.25 2.25;5 12]
stepSize=5;
c=0
for co=1:numel(smpCol)
    intSmp=smpCol{co}
    if contains(smeta.tag2(intSmp(1)),{'2'})
        lineColor=[100,100,100;127,63,152]/256;
    else
        lineColor=[100,100,100;0,174,239]/256;
    end
    for p=1:numel(intPar)
        parData=eval(intPar{p});
        for n=1:numel(intNuc)
            c=c+1
            selNucs=find(ismember(nucMeta.order,intNuc{n}) & all(~isnan(dynNuc(:,intSmp)),2) &ismember(nucMeta.gene,find(~isnan(parData) &selGenes)));
            [xVal,idx]=sort(parData(nucMeta.gene(selNucs)));
            selNucs=selNucs(idx);
            %subplot(2,numel(intPar)*numel(intNuc),(p-1)*numel(intNuc)+n)
            %subplot(2,numel(intPar)*numel(intNuc),p+numel(intPar)*(n-1))
            subplot(2,3,c)
            hold off
            clear lineObj
            for s=1:numel(intSmp)
                lineName=sprintf('%s-%d',strjoin(smeta{intSmp(s),[7,8,9]}),smeta.tp(intSmp(s)))
                corr(xVal,dynNuc(selNucs,intSmp(s)))
                yVal=movmean(dynNuc(selNucs,intSmp(s)),movPar(p),'omitnan','SamplePoints',makeUnique(xVal,4));
                yValE=movstd(dynNuc(selNucs,intSmp(s)),movPar(p),'omitnan','SamplePoints',makeUnique(xVal,4))./sqrt(movsum(~isnan(dynNuc(selNucs,intSmp(s))),movPar(p),'omitnan','SamplePoints',makeUnique(xVal,4)));
                lineObj(s)=plot(xVal(1:stepSize:end),yVal(1:stepSize:end),'-','Color',lineColor(s,:),'DisplayName',lineName)
                hold on
                fill([xVal(1:stepSize:end);flipud(xVal(1:stepSize:end))],[yVal(1:stepSize:end)-yValE(1:stepSize:end);flipud(yVal(1:stepSize:end)+yValE(1:stepSize:end))],brighten(lineColor(s,:),0),'LineStyle','none','FaceAlpha',.3)
            end
            axis tight
            xlim(parLim(p,:))
            xlabel(intPar{p})
            ylabel(sprintf('movavg %.2f',movPar(p)))
            legend(lineObj)
            title(intNuc(n))
        end
    end
end
save_gf(gcf,sprintf('Fig3Lines_DynRel_geneBody',smeta.tag2{intSmp(1)},smeta.tp(intSmp(1))))
close all

%% summary matrix

intGenes={find(relChg>1 & expClu==1);
    find(relChg<-1 & expClu==2);
    GP.groups{6}{2}{1};
    GP.groups{7}{2}{66}}
intNuc={[-2,-1],[2]}
sumMat=cell(1,numel(intGenes))
for g=1:numel(intGenes)
    sumMat{g}=nan(numel(intNuc),size(dynNuc,2));
    for n=1:numel(intNuc)
        selNucs=ismember(nucMeta.gene,intGenes{g}) & ismember(nucMeta.order,intNuc{n}) & ~ismember(nucMeta.gene,GP.groups{23}{2}{45});
        sumMat{g}(n,:)=smoothdata(median(dynNuc(selNucs,:),'omitnan'),'gaussian',3);
    end    
end
intTp=[-1;unique(smeta.tp(smeta.tp>0))]
selTc=[2,4,1,3]
nucNames={'promoter (-2,-1)','geneBody (2,3)'}
geneNames={'induced','repressed','Stress genes','ribi'}
for g=1:numel(intGenes)
    for n=1:numel(intNuc)
        subplot(2,2,n+2*(g-1))
        tcMat=nan(numel(intTp),numel(selTc));
        for t=1:numel(selTc)
            subplot(4,1,t)
            selSmp=find(smeta.tcId==selTc(t) &ismember(smeta.tp,intTp));
            [~,idx]=ismember(smeta.tp(selSmp),intTp);
            tcMat(idx,t)=sumMat{g}(n,selSmp)
        end
        imagesc(tcMat)
        xticks(1:numel(selTc))
        xticklabels(strcat(tcs.ab(selTc),'-',tcs.tag2(selTc)))
        title(sprintf('%s-%s',geneNames{g},nucNames{n}))
        colorbar()
        caxis(max(abs(caxis()))*[-1 1])
        yticks(1:numel(intTp))
        yticklabels(intTp)
    end    
end
save_gf(gcf,sprintf('Fig3Matrix_allNUcleosomeTypes',smeta.tag2{intSmp(1)},smeta.tp(intSmp(1))))
close all
g=1
n=2
figure
tcMat=nan(numel(intTp),numel(selTc));
load('greytopurple.mat');
load('greytoblue.mat')

cLim={[-0.13 0.13],[-0.25 0.25],0.13*[-1 1],0.4*[-1 1]}
figure
for g=1:2
    subplot(5,2,g)
    plot(ExpData(1).time,median(dynExp{1}(intGenes{g},:),'omitnan')-median(dynExp{1}(intGenes{g},2),'omitnan'),'-')
    xlim([-1.5 69.5])
    ylabel('RNA')
    colorbar()
    title(geneNames{g})
    for t=1:numel(selTc)
        subplot(5,2,2*(t+1)+g-2)
        selSmp=find(smeta.tcId==selTc(t) &ismember(smeta.tp,intTp));
        [~,idx]=ismember(smeta.tp(selSmp),intTp);
        tcMat(idx,t)=sumMat{g}(n,selSmp)
        imagesc(tcMat(:,t)')
        ylabel(strcat(tcs.ab(selTc(t)),'-',tcs.tag2(selTc(t))))
        colorbar()
        if contains(tcs.tag2(selTc(t)),'h2')
            %colormap(gca,Colo(29+[-28:28],:))
            colormap(gca,yellowH2)
        else
           colormap(gca,yellowH3b)
           %colormap(gca,Colo1(40+[-24:24],:))
        end
        xticks([])
        if t==4
            xticks(1:numel(intTp))
            xticklabels(intTp)
        end
        caxis(cLim{t})
    end
end
set(gcf,'Color',[1 1 1])
save_gf(gcf,sprintf('Fig3Matrix_%s_%s',geneNames{g},extractBefore(nucNames{n},' ')))
end




colNucs=
%% correlation 
for n=1:numel(intNuc)
    figure
    for g=1:numel(geneGroups)
        selNucs=nucMeta.order==intNuc(n) & ismember(nucMeta.gene,find(expClu==g)); %geneGroups{g}
         subplot(1,2,g)
        imagesc(corr(dynNuc(selNucs,:),[geneLvl(nucMeta.gene(selNucs)),relChg(nucMeta.gene(selNucs))],'rows','pairwise'),[-.3 .3])
        yticks(1:size(smeta,1))
        yticklabels(strcat(smeta.tag2,smeta.ab,sprintfc('%d',smeta.tp)))
        title(sprintf('%s %d',groupNames{g},sum(selNucs)))
        xticks(1:2)
        xticklabels({'abs','rel'})
    end
    colorbar()
    suptitle(num2str(intNuc(n)))
end
%% correlation all tcs all nucs
%async=load('forFigure2add.mat');
dynNuc=nan(size(nucData));
asyncCon=[3,14,5,16];
for i=1:max(smeta.tcId)
   dynNuc(:,smeta.tcId==i&smeta.bad==0)=log(nucData(:,smeta.tcId==i&smeta.bad==0)+1)-log(median(nucData(:,smeta.tcId==i&smeta.bad==0 & smeta.tp<=3)+1,2,'omitnan'));%log(medNuc(:,asynCon(i))+1);;
   %dynNuc(:,smeta.tcId==i&smeta.bad==0)=log(nucData(:,smeta.tcId==i&smeta.bad==0)+1)-log(async.medNuc(:,asyncCon(i))+1);
end
intNuc={-2,-1,1,2:100};
blackList=ismember(nucMeta.gene,GP.groups{23}{2}{45});
geneGroups={find(expClu==2 & relChg<-1 ),find(expClu==1 & relChg>1)};
groupNames={'Rep','Ind'}
crMat=nan(size(dynNuc,2),numel(intNuc).*2,numel(geneGroups));
for n=1:numel(intNuc)
    for g=1:numel(geneGroups)
        selNucs=find(ismember(nucMeta.gene,geneGroups{g}) & ismember(nucMeta.order,intNuc{n}) &~blackList); 
        crMat(:,n+[0,numel(intNuc)],g)=corr(dynNuc(selNucs,:),[relChg(nucMeta.gene(selNucs)) geneLvl(nucMeta.gene(selNucs))],'rows','pairwise');
    end
end
close all
figure;
border=cumsum([.5;accumarray(smeta.tcId(smeta.bad==0),1)])
parCols={1:4,5:8}
parName={'relChg','geneLvl'}
c=0;
for g=1:numel(geneGroups)
    for p=1:numel(parCols)
        c=c+1;
        subplot(1,numel(parCols)*numel(geneGroups),c)
        imagesc(crMat(smeta.bad==0,parCols{p},g),[-.35 .35])
        title([groupNames{g} ' vs. ' parName{p}])
        xticks(1:4)
        xticklabels(repmat(intNuc,1,1))
        yticks(1:size(smeta,1))
        yticklabels(strcat(smeta.tag2(smeta.bad==0),smeta.ab(smeta.bad==0),sprintfc(' %d',smeta.tp(smeta.bad==0))))    
        hold on 
        plot(xlim'*ones(1,numel(border)),border'.*[1;1],'k--')
        colorbar()
    end
end

%% only show H3 in ORF for 
selNuc=[4]
selTc=[3,4]
for g=1:numel(geneGroups)
    for i=1:numel(selTc)
        selSmp=ismember(smeta.tcId,selTc(i)) & smeta.bad==0;
        subplot(numel(selTc),numel(geneGroups),(i-1).*numel(geneGroups)+g)
        imagesc(crMat(selSmp,selNuc+[0 numel(intNuc)],g)')
        title(sprintf('%s-%s',strjoin(tcs{selTc(i),:}),groupNames{g}))
        xticks(1:numel(selSmp))
        xticklabels(smeta.tp(selSmp))
        yticks(1:2)
        yticklabels({'relChg','geneLvl'})
    end
end
save_gf(gcf,'Fig3_AllTcCorr')

%% RNA plots
GP=load('group_imp.mat');
figure
geneGroups={GP.groups{6}{2}{1};
    GP.groups{7}{2}{66};
    find(expClu==1 & relChg>1);
    find(expClu==2 & relChg<-1 & geneLvl<10);
    find(expClu==2 & relChg<-1 & geneLvl>10)
    };
groupNames={'stress','ribi','ind.','rep (low exp)','rep (high exp)'}
for g=1:numel(geneGroups)
    plot(ExpData(1).time,median(dynExp{1}(geneGroups{g},:),'omitnan'),'-','LineWidth',2,'DisplayName',groupNames{g})
    hold on
end
xlabel('Time in H2o2')
ylabel('median expression change (log2)')
legend()
axis('tight')
save_gf(gcf,'Fig3_mRNA')
%% profiles individual
clearvars -except expclu relChg dynGenes expClu
[data,smeta]=getNaamaData('figure',3);
smeta.tp=str2double(regexp(smeta.con,'(?<=h202T)\d+','match','once'));
smeta.tp(isnan(smeta.tp),:)=-1;
GP=load('group_imp.mat');
[~,idx]=sortrows(smeta,{'tag2','ab','tp'});
smeta=smeta(idx,:);
data=data(:,idx);
metaProfile=meta_profile(data,'promoter',700,'useORF',false,'afterTSS',500,'scaled',false,'disCDS',false);
promLvl=squeeze(sum(metaProfile.cube(:,1:700,:),2));
orfLvl=squeeze(sum(metaProfile.cube(:,702:end,:),2));

geneIdx=find(relChg>1 & std(promLvl(:,strcmp(smeta.ab,'myc')&strcmp(smeta.tag2,'h3')),[],2,'omitnan')>2000)
close all
for i=1:numel(geneIdx)
    if mod(i-1,20)==0
        figure
    end
    subplot(5,4,mod(i-1,20)+1)
    imagesc(squeeze(metaProfile.cube(geneIdx(i),:,smeta.bad==0&strcmp(smeta.ab,'myc')&strcmp(smeta.tag2,'h3')))')
    title(GP.gene_infoR64.name(geneIdx(i)))
end
tssIdx=nan(6701,1);
hasTss=~isnan(GP.gene_infoR64.felixTss(:,1));
tssIdx(hasTss)=GP.chrIdx(GP.gene_infoR64.felixTss(hasTss,1))+GP.gene_infoR64.felixTss(hasTss,2);%
intGene=mat2cell(ans,ones(size(ans)),1)%{[GP.gene_table.HSP12]}
%intGene={GP.gene_table.PWP2,GP.gene_table.SPB1,GP.gene_table.URA7}
close all
intCmb={[16:21,27],[47,52,58,59,60,61,62]};
lineColor=lines(7);
c=0;
for g=1:numel(intGene)
    figure
    c=0;
    selRegion=quantile(acol(tssIdx(intGene{g})+[-1000,1000]),[0 1]);
    for ce=1:numel(intCmb)
        intExp=intCmb{ce}     
        c=c+1;
        subplot(2,2,c)
        hold off        
        for e=1:numel(intExp)
            plot(selRegion(1):selRegion(2),conv2(data(selRegion(1):selRegion(2),intExp(e)),[1:50 49:-1:1]'.^2./sum([1:50 49:-1:1].^2),'same'),'LineWidth',2,'DisplayName',strjoin(smeta{intExp(e),[1]}),'Color',lineColor(e,:))
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
%% summarize nucleosome dynamics in 
%% look at individual profiles
clear all
GP=load('group_imp.mat');
load('forFigure3add.mat', 'dynGenes','expClu','relChg')
[data,smeta]=getNaamaData('figure',3);
smeta.tp=str2double(regexp(smeta.con,'\d+$','match','once'));
smeta.tp(isnan(smeta.tp),:)=-1;
[~,idx]=sortrows(smeta,{'tag2','ab','tp'});
smeta=smeta(idx,:);
data=data(:,idx);
[tcs,~,smeta.tcId]=unique(smeta(:,[14 7]),'stable');
metaProfile=meta_profile(data,'promoter',600,'useORF',false,'afterTSS',600,'scaled',false,'disCDS',false);
geneGroups={GP.groups{7}{2}{66}(expClu(GP.groups{7}{2}{66})==2 & relChg(GP.groups{7}{2}{66})<-1),GP.groups{6}{2}{1}(expClu(GP.groups{6}{2}{1})==1 & relChg(GP.groups{6}{2}{1})>1)}
intTcs=[1,2,3,4]
groupNames={'RiBi','Stress'}
close all
figure
for g=1:numel(geneGroups)    
    for i=1:numel(intTcs)
        subplot(numel(geneGroups),numel(intTcs),(g-1)*numel(intTcs)+i)
        selSmp=smeta.tcId==i & smeta.bad==0;
        meanProfile=squeeze(mean(metaProfile.cube(geneGroups{g},:,selSmp),'omitnan'));
        imagesc(meanProfile','Xdata',[-600 600])
        colorbar()
        yticks(1:sum(selSmp))
        yticklabels(smeta.tp(selSmp))
        title(strjoin([tcs{i,:} groupNames{g}]))
    end
end
save_gf(gcf,'Fig3_GroupProfile')

    
end