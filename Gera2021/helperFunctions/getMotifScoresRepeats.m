function out=getMotifScoresRepeats(forProm,TFNames)
    fileNames = dir('/home/labs/barkailab/tamarge/Master/TF ohnologs new Table*summary.xlsx');
    keepCol = 'goodSamples';
    for i = 1:numel(fileNames)
        temp = readtable([fileNames(i).folder,'/',fileNames(i).name]);
        goodRep{i} = temp.(keepCol);
    end
    goodRep = cat(1,goodRep{:});
    goodRep = goodRep(~cellfun('isempty', goodRep));
    goodRep = cellfun(@(x) strsplit(x,{',',' '}), goodRep,'UniformOutput',false);
    goodRep = cat(2,goodRep{:});
    clearvars -except forProm TFNames goodRep
    load('checAllver5.mat')
    FN=fieldnames(checAllver5.norm);
    for tf=1:numel(TFNames)
        TFRep{tf}=goodRep(contains(goodRep,[TFNames{tf},'_rep'])&ismember(goodRep,FN));
        if numel(TFRep{tf})==0
            TFRep{tf}=goodRep(cellfun('prodofsize',regexp(lower(goodRep),[lower(TFNames{tf}),'_\d'],'once'))>0&ismember(goodRep,FN));
        end
    end
    maxSize=lcm(numel(TFRep{1}),numel(TFRep{2}));
    if 1==1
        adjRep=maxSize./cellfun('prodofsize',TFRep)
    else
        adjRep=repmat(1,1,numel(TFNames));
    end
    for tf=1:numel(TFNames)    
        normVec=nan(12071326,numel(TFRep{tf}));
        for r=1:numel(TFRep{tf})
            normVec(:,r)=cat(1,checAllver5.norm.(TFRep{tf}{r}){:});
        end
        for g=1:numel(forProm)
            for p=1:numel(forProm(g).pos)
                if size(normVec,2)>0
                out{g,tf}{p}=reshape(repmat(sum(normVec(forProm(g).intPos{p},:))./sum(normVec(forProm(g).pos{p},:)),adjRep(tf),1),1,[]);     
                else
                    out{g,tf}{p}=[];
                end
            end
        end
    end
end