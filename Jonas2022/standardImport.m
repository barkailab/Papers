function [expTypes,medNuc,smeta,nucData,nucMeta]=standardImport(n,varargin)
ip=inputParser;
ip.addParameter('save',false);
ip.addParameter('factor',true);
ip.addParameter('Pugh',[]);
ip.addParameter('Chec',[]);
ip.parse(varargin{:});
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/nucMeta015AddFried.mat')
[data,smeta]=getNaamaData('figure',n);
h3B={'fast1b-wt-async-x18','fast2b-wt-async-x18','fast2a-wt-async-x22','fast2b-wt-async-x22'};

%oldH3=readtable('oldH3.txt','ReadVariableNames',false);
%smeta.bad(strcmp(smeta.tag2,'h3')&contains(smeta.gt,'wt')&contains(smeta.ab,{'myc','ha'})&contains(smeta.con,'async')&~contains(smeta.name,oldH3.Var1))=1;
smeta.gt(contains(smeta.gt,'rlf2'))={'rlf2'};
if numel(ip.Results.Pugh)>0
    allPugh=dir('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Rossi2021/*.mat');
    selProfiles=contains({allPugh.name},ip.Results.Pugh,'IgnoreCase',true);
    allPugh=allPugh(selProfiles,:)
    addPugh=zeros(size(data,1),numel(allPugh));
    for i=1:numel(allPugh)
        temp=load([allPugh(i).folder,'/' allPugh(i).name]);
        addPugh(:,i)=temp.normProfile(1:12157105);
        totPugh(i,1)=temp.totalReads;
    end
    allPugh=struct2table(allPugh);
    allPugh.totReads=totPugh;
    allPugh.tag2(:)={'none'};
    allPugh.tag(:)={'none'};
    allPugh.tagid(:)=0;
    allPugh.bad(:)=false;
    allPugh.gt(:)={'wt'};
    allPugh.con(:)={'async'};
    allPugh.ab=regexp(allPugh.name,'(?<=_).*(?=.mat)','match','once');
    allPugh.exp=regexp(allPugh.name,'^\d+(?=_)','match','once');
    smeta=[smeta;allPugh];
    data=[data,addPugh];
    clearvars -except data nucMeta smeta ip n
end
nBase=size(data,1);
nucExt=50;
nucPos=zeros(nBase,1);
for i=1:size(nucMeta,1)
    posi=nucMeta.pos(i)-nucExt:nucMeta.pos(i)+nucExt;
    if i>1
        posi=posi(posi>0 & posi>mean(nucMeta.pos([i-1,i])) & posi<=nBase);
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
clear data i posi h3B nucPos nucExt
smeta.tag2(cellfun('isempty',smeta.tag2))={'notag'}
highOcc=find(mean(movmean(nucData(:,strcmp(smeta.ab,'ha')& ~contains(smeta.tag2,{'gla','par'})& strcmp(smeta.gt,'wt') &smeta.bad==0),5)>30,2)>0.05);
highOcc=highOcc(highOcc>2);
nucData(unique(acol(highOcc+[-2:2])),:)=NaN;
lowOcc=find(mean(movmean(nucData(:,strcmp(smeta.ab,'ha') & ~contains(smeta.tag2,{'gla','par'})&strcmp(smeta.gt,'wt') &smeta.bad==0),3)<0.01,2)>0.05);
lowOcc=unique(acol(lowOcc+[-1:1]));
lowOcc=lowOcc(lowOcc>0);
nucData(lowOcc,:)=NaN;
clearvars  highOcc lowOcc
[smeta,idx]=sortrows(smeta,[7,14,15,9,10]);
nucData=nucData(:,idx);
clear idx
smeta.name(smeta.bad==1)=strcat(smeta.name(smeta.bad==1),'*');

figure;imagesc(corr(nucData,'rows','pairwise'))
labelSmp=1:size(nucData,2)%find(contains(smeta.tag2,'hht2'))%find(contains(smeta.con,{'async','alpha'}));%true(size(nucData,2),1))%find(contains(smeta.tag2,{'h3x2','h3y2'})&contains(smeta.con,{'async','alpha'}));
yticks(labelSmp)
yticklabels(strcat(smeta.ab(labelSmp),';',smeta.tag2(labelSmp),';',smeta.con(labelSmp),':',extractBefore(smeta.exp(labelSmp),'.'))) 


xticks(labelSmp)
xticklabels(strcat(smeta.ab(labelSmp),';',smeta.tag2(labelSmp),';',smeta.gt(labelSmp),';',smeta.con(labelSmp))) 
xtickangle(90)

smeta.bad(ismember(smeta.exp,{'x27.mat'})&contains(smeta.gt,'spt'))=1
[expTypes,~,smeta.expId]=unique(smeta(:,[7,9,10,14]),'stable')
[~,~,smeta.sid]=unique(smeta(:,[8:11]),'stable')
medNuc=nan(size(nucData,1),max(smeta.expId));
for i=unique(smeta.expId(smeta.bad==0 ))'   
    selRpt=find(smeta.expId==i & ~smeta.bad);
    selNucs=all(nucData(:,selRpt)>1,2);
    [~,bestRpt]=max(sum(corr(log(nucData(selNucs,selRpt))),2));
    bestRpt=selRpt(bestRpt);
    for j=selRpt'
        pj=median(diff(log(nucData(selNucs,[j bestRpt])),1,2));
        fac(j)=exp(pj);
    end
    if ip.Results.factor
        medNuc(:,i)=median(nucData(:,selRpt).*fac(selRpt),2);
    else
        medNuc(:,i)=median(nucData(:,selRpt),2);
    end
    %stdNuc(:,i)=std(log(nucData(:,selRpt).*fac(selRpt)+1),[],2);  
end
clear fac bestRpt i j nBase pj selNucs selRpt oldH3
if ip.Results.save
    if n<20
        save(sprintf('~/Documents/MATLAB/projects/Naama/forFigure%d.mat',n))
    else
        save(sprintf('~/Documents/MATLAB/projects/Gilad-CC/forFigure%d.mat',n))
    end
end