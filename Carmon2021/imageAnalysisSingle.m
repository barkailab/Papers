function imageAnalysisSingle(cziFile,varargin)
ip=inputParser();
ip.addParameter('fixedNuclei',false);
ip.addParameter('mrna',false)
ip.parse(varargin{:});
probeTable=readtable('probes.xlsx');
fileIdx=dir(cziFile);
fileIdx=find(contains(probeTable.Var1,fileIdx.name));
reader=bfGetReader(cziFile);
meta=reader.getMetadataStore;
xml=meta.dumpXML;
xml=xml.toCharArray';
if contains(cziFile,'czi')
    channels=regexp(xml,'(?<=Fluor=").*?(?=")','match');
else
    channels=regexp(xml,'(?<=\<Channel).*?(?=/Channel)','match') ;
    channels=strcat('E',regexp(channels,'(?<=EmissionWavelength\=").*?(?=")','match','once'));
    if ~contains('DAPI',channels)
        channels=strrep(channels,'TMR','DAPI')
    end
end
[~,channelsIdx]=ismember(genvarname(channels),probeTable.Properties.VariableNames)
probes=cell(size(channels));
if numel(fileIdx)==1
    probes(channelsIdx>0)=probeTable{fileIdx,channelsIdx(channelsIdx>0)}
    probes(channelsIdx==0)={''}
else
    probes(channelsIdx>0)=channels(channelsIdx(channelsIdx>0));
    probes(channelsIdx==0)={''}
end
th=cell(size(channels));
ts=cell(size(channels));
%% define nuclei
tic
dapiCh=find(contains(channels,'DAPI')|contains(probes,'DAPI'));
dapi=double(getMaxProject(reader,dapiCh));
smoothDapi=imfilter(dapi,fspecial('disk',3),'circular');
bg=imfilter(dapi,fspecial('disk',100),'circular');
sD=smoothDapi-bg;
nObj=zeros(100,1);
nGood=zeros(100,1);
if ip.Results.fixedNuclei
    if ismember(reader.getSizeX,[3104,2084,2048])
        nucArea=25.^2*pi;
    elseif reader.getSizeX==2443
        nucArea=21.^2*pi;
    elseif reader.getSizeX==1024
        nucArea=13.^2.*pi;
    end
else
    nucAreaColl=zeros(100,1);
end
for j=1:100
    [lab,nObj(j)]=bwlabel(sD>(j/500.*quantile(sD(:),0.98)));
    labArea=accumarray(lab(lab>0),1);
    if ~ip.Results.fixedNuclei
        [~,idx]=max(movmean(histcounts(log(labArea),4:0.2:12,'Normalization','countdensity'),3));
        nucAreaColl(j)=exp(3.9+idx*0.2);
        nucArea=nucAreaColl(j);
    end    
    nGood(j)=sum(labArea>nucArea*0.66 & labArea<nucArea*1.5);      
end
[~,bestTh]=max(log(nGood)-log(nObj));
th{dapiCh}=bestTh;
if ~ip.Results.fixedNuclei
        nucArea=nucAreaColl(bestTh);
end
bestTh=bestTh/500.*quantile(sD(:),0.98);

nucleiMask=bwlabel(sD>bestTh);
nucleiMaskOld=nucleiMask;
clear bg sD
nuclei=regionprops('table',nucleiMask,'ConvexHull','ConvexArea','Centroid','Area');
% %% divide big nuclei
% for j=find(nuclei.Area>nucArea*1.5)'
%     temp=smoothDapi;
%     temp(nucleiMask~=j)=NaN;
%     clear nObj nGood
%     for k=1:50
%         lab=bwlabel(temp>(k/100*range(temp(temp>0))+min(temp(temp>0))));
%         objArea=accumarray(lab(lab>0),1);
%         nGood(k)=sum(objArea>nucArea/2 & objArea<nucArea*1.5);
%         if k>1
%             if nGood(k-1)>nGood(k)
%                 break
%             end
%         end
%     end
%     bestTh=find(nGood==max(nGood),1)/100*range(temp(temp>0))+min(temp(temp>0));
%     lab=bwlabel(temp>bestTh)+max(nucleiMask(:));
%     objArea=accumarray(lab(lab>0),1);
%     lab(~ismember(lab,find(objArea>nucArea/2 & objArea<nucArea*1.5)))=0;
%     nucleiMask(nucleiMask==j)=lab(nucleiMask==j);
% end
% nuclei=regionprops('table',nucleiMask,'ConvexHull','ConvexArea','Centroid','Area');
%badNuclei=nuclei.ConvexArea<nucArea/2|nuclei.ConvexArea>nucArea*1.5;
%nucleiMask(ismember(nucleiMask,find(badNuclei)))=0;
%nucleiMask=changem(nucleiMask,1:sum(~badNuclei),find(~badNuclei));
toc
%nuclei=regionprops('table',nucleiMask,'ConvexHull','ConvexArea','Centroid','Area');
clear badNuclei;
for j=setdiff(1:numel(channels),dapiCh)
    %% define spots
    if ip.Results.mrna
        upperTh=min([0.2,str2double(regexp(probes{j},'\d+.\d+','match','once'))]);
        manualTh=str2double(regexp(probes{j},'(?<=_mTh)0.\d+','match','once'));
        manualTh=manualTh(~isnan(manualTh));
        if numel(probeTable.Zlimit{fileIdx})>0
            [~,zLimit,~]=regexp(probeTable.Zlimit{fileIdx},'(\d+)','tokens','match');
            zLimit=str2double(zLimit);
            [locs,th{j}]=mrnaAnalyse(reader,j,upperTh,'zLimit',zLimit)
        else
            [locs,th{j}]=mrnaAnalyse(reader,j,upperTh,'manualTh',manualTh)
        end
        
        
    else
        [locs,th{j}]=tssAnalyse(reader,j)
    end
    %%assign spots to nuclei
    distMat=pdist2(locs.Centroid,nuclei.Centroid);
    for k=1:size(locs,1)
        [~,closeNuclei]=mink(distMat(k,:),5);
        minDis=nan(5,1);
        for kk=1:numel(closeNuclei)
            [y,x]=ind2sub(size(nucleiMask),find(nucleiMask==closeNuclei(kk)));
            minDis(kk)=min(pdist2(locs.Centroid(k,:),[x,y]));
        end
        bestNuc=closeNuclei(find(minDis==min(minDis),1));
        try
            locs.nuc(k)=bestNuc;
        catch
            k
            kk
        end
        if inpolygon(locs.Centroid(k,1),locs.Centroid(k,2),nuclei.ConvexHull{bestNuc}(:,1),nuclei.ConvexHull{bestNuc}(:,2))
            locs.nucDis(k)=0;
        else
            locs.nucDis(k)=min(minDis);
        end
    end
    ts{j}=locs;
end
out=[];
out.nuclei=nuclei;
out.nucleiMask=nucleiMask;
out.ts=ts;
out.th=th;
out.nucArea=nucArea;
try
    out.midLine=load(strrep(cziFile,'.czi','f4.mat'),'midLine');
end
if contains(cziFile,'nd2')
    saveStruct(strrep(cziFile,'.nd2','f6.mat'),out)
else
    saveStruct(strrep(cziFile,'.czi','f6.mat'),out)
end
end


function [locs,th]=tssAnalyse(reader,j)
pic=double(getMaxProject(reader,j)); %double(bfGetPlane(reader,j));
H = -fspecial('log',15,2.5);
ims=imfilter(pic,H,'circular');
ims=ims/max(ims(:));
clear nObj nGood
for t=1:50    
    [lab,nObj(t)]=bwlabel(ims>(t/100));
    objArea=accumarray(lab(lab>0),1);
    nGood(t)=sum(objArea>20&objArea<70);
    bgMedian(t)=median(pic(imdilate(lab>0,strel('disk',30)) & lab==0));
    bgStd(t)=std(pic(imdilate(lab>0,strel('disk',30)) & lab==0));
    sigMedian(t)=median(pic(lab>0));
    sigStd(t)=std(pic(lab>0));
end
nObj(1:find(nObj==max(nObj))-1)=NaN;
logChg=diff(movmean(log(nObj),3,'Endpoints','discard'),1);
logChg(1:find(logChg==min(logChg)))=NaN;

idx1=find((logChg(1:end-5)./logChg(6:end)<1.5)&logChg(1:end-5)>min(logChg/2)&sigMedian(2:end-7)>(bgMedian(2:end-7)+2*bgStd(2:end-7)),1)+1;
%[~,idx]=min(log(nObj)-log(nGood));
[~,idx2]=max(movmean(nObj,5,'omitnan','Endpoints','discard')./(movstd(nObj,5,'omitnan','Endpoints','discard')+10));
idx2=idx2+2;
if numel(idx1)>0
    idx=idx1;
else
    idx=idx2;
end
th=idx/100;
lab=bwlabel(ims>(idx/100));
locs=regionprops('table',lab);
locs.int=accumarray(lab(lab>0),pic(lab>0));
end

function [locs,th]=mrnaAnalyse(reader,channel,upperTh,varargin);
ip=inputParser;
ip.addParameter('zLimit',[1 35]);
ip.addParameter('manualTh',[]);
ip.addParameter('thFigure',false);

ip.parse(varargin{:});
zLimit=ip.Results.zLimit;
zLimit(2)=min(zLimit(2),reader.getSizeZ);
imageStack=zeros(reader.getSizeY,reader.getSizeX,range(zLimit)+1);
c=0;
for z=zLimit(1):zLimit(2)
    c=c+1;
    imageStack(:,:,c)=bfGetPlane(reader,reader.getIndex(z-1,channel-1,0)+1); 
end
bg=imboxfilt3(imageStack,[51,51,5]); %%background correction
imageStackBg=imageStack-bg;
ims=LOG_filter(double(imageStackBg),'fspecial',3);
if numel(ip.Results.manualTh)==0
    th=count_mrna(ims,upperTh,'connect',18,'showFigure',true);
else
    th=ip.Results.manualTh;
end
lab = threshstack(ims,th,'connect',18);
maxVecs=accumarray(lab(lab>0),imageStackBg(lab>0),[],@max);
maxImage=zeros(size(lab));
maxImage(lab>0)=maxVecs(lab(lab>0));
lab=lab.*(imageStackBg>(maxImage*0.1));
clear maxImage maxVecs
locs=regionprops3(lab,imageStackBg,'Volume','MaxIntensity','VoxelValues','MinIntensity','WeightedCentroid');
locs.Properties.VariableNames{contains(locs.Properties.VariableNames,'Centroid','IgnoreCase',true)}='Centroid';
locs.zPos=locs.Centroid(:,3)+zLimit(1)-1;
locs.Centroid=locs.Centroid(:,1:2);
if ip.Results.thFigure
    figure;
    subplot(1,3,1)
    imagesc(max(ims,[],3)>0.02*max(ims,[],'all'))
    subplot(1,3,2)
    imagesc(max(ims,[],3))
    caxis([0 5])
    subplot(1,3,3)
    imagesc(max(imageStackBg,[],3))
    linkaxes
    caxis([0 5])
    hold on
    scatter(locs.Centroid(:,1),locs.Centroid(:,2),'r.')
end
locs.int=accumarray(lab(lab>0),imageStackBg(lab>0));
locs.intRaw=accumarray(lab(lab>0),imageStack(lab>0));

locs=locs(locs.int>0&locs.Volume>1,:);
if iscell(locs.MaxIntensity)
    locs.MaxIntensity=cat(1,locs.MaxIntensity{:});
end

end

function checkMrna(pictype)
% figure;
% zSlide=uicontrol('Style','slider','Units','normalized','Position',[0.1 0.05 0.8 0.05],'Min',1,'Max',size(imageStackBg,3),'Value',1)
close all
c=0;
for i=1:size(imageStackBg,3)
        if c==9
        figure;
        c=0;
        end
        
        c=c+1;

    subplot(3,3,c)
    imagesc(imageStackBg(450:550,450:550,i),'XData',[450 550],'YData',[450 550])
    caxis([0 500])
    hold on
    lowDot=ismember(round(locs.zPos),i-3:i-1);
    highDot=ismember(round(locs.zPos),i+1:i+3);
    scatter(locs.Centroid(lowDot,1),locs.Centroid(lowDot,2),'k.')
    scatter(locs.Centroid(round(locs.zPos)==i,1),locs.Centroid(round(locs.zPos)==i,2),'ro')
    scatter(locs.Centroid(highDot,1),locs.Centroid(highDot,2),'r.')    
end
end