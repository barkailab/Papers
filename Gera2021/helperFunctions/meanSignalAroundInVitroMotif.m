function out = meanSignalAroundInVitroMotif(checStruct, TF, varargin)
ip = inputParser;
ip.addParameter('nmer',5);
ip.addParameter('windowSize',150);
ip.addParameter('motif',{});
ip.addParameter('fimo',false);
fimoTh=10^-4;
ip.parse(varargin{:});
TFtable = readtable('/home/labs/barkailab/tamarge/Master/TF ohnologs new table.xlsx');
TFtable = TFtable(1:42,:);
load('bestCISBPInVitroWM.mat');
 load('SC_genome.mat');
 SC_genome = upper(cat(2,SC_genome(1:16).Sequence));
 load('promoterIDXvec.mat');
 load('/home/labs/barkailab/tamarge/Master/mat Files/allOccFimo.mat')
 
 promoterIDXvec = logical(promoterIDXvec);
%%

nmer = ip.Results.nmer;
windowSize = ip.Results.windowSize;
nt = {'A','C','G','T'};
[rowIdx, paraIdx] = ind2sub([42,2],find(contains(TFtable{:,4:5}, upper(TF))));
if ~iscell(TF)
    TFpara = TFtable{rowIdx, 6-paraIdx}{1};
    TFs = {[TFtable{rowIdx, paraIdx+3}{1}(1), lower(TFtable{rowIdx, paraIdx+3}{1}(2:end))], [TFpara(1), lower(TFpara(2:end))]};
else
    TFs = TF;
end
for p = 1:2
    TF = TFs{p};
    TFpara = TFs{3-p};
    if isempty(ip.Results.motif)
        if isfield(bestCISBPInVitroWM, upper(TF))
            currPwm1 = bestCISBPInVitroWM.(upper(TF));
            flag = '';
        elseif isfield(bestCISBPInVitroWM, upper(TFpara))
            currPwm1 = bestCISBPInVitroWM.(upper(TFpara));
            flag = '*';
        else
            currPwm1 = [];
        end
    else
        currPwm1 = [];
    end
    
    if numel(currPwm1) > 0 | ~isempty(ip.Results.motif)
        if isempty(ip.Results.motif)
            [vBase,idxBase] = max(currPwm1);
            fullLogo = nt(idxBase);
            [vPos, idxPos] = maxk(vBase,nmer);
            fullLogo(setdiff(1:length(fullLogo), idxPos)) = {'N'};
            infLogo = fullLogo(min(idxPos):max(idxPos));
            infLogo = cat(2, infLogo{:});
            rcInfLogo = seqrcomplement(infLogo);
            infLogo = strrep(infLogo, 'N','.');
            rcInfLogo = strrep(rcInfLogo, 'N','.');
            if strcmp(flag, '')
                out.motif{p} = infLogo;
            else
                out.motif{p} = '';
            end
        else
            infLogo = ip.Results.motif{1};
            rcInfLogo = seqrcomplement(infLogo);
        end
        motifOcc = regexp(SC_genome,  sprintf('%s', infLogo))+floor(length(rcInfLogo)/2);
        motifOccRC = regexp(SC_genome,  sprintf('%s',rcInfLogo))+floor(length(rcInfLogo)/2);
        
        motifOcc = motifOcc(promoterIDXvec(motifOcc));        
        motifOccRC = motifOccRC(promoterIDXvec(motifOccRC));
        if ip.Results.fimo
            motifOcc=allOccFimo(contains(allOccFimo.motif_id,upper(TF))&allOccFimo.chrId<17&contains(allOccFimo.strand,'+')&allOccFimo.p_value<=fimoTh,:);
            motifOcc=motifOcc.center(promoterIDXvec(motifOcc.center))';
            motifOccRC=allOccFimo(contains(allOccFimo.motif_id,upper(TF))&allOccFimo.chrId<17&contains(allOccFimo.strand,'-')&allOccFimo.p_value<=fimoTh,:);
            motifOccRC=motifOccRC.center(promoterIDXvec(motifOccRC.center))';
        end
        for tf = 1:numel(TFs)
            fullProfile = chromosomes2fullProfile(checStruct, TFs(tf))';
            signalMat = fullProfile([motifOcc'+[-windowSize:windowSize]; motifOccRC'-[-windowSize:windowSize]]);
            out.aroundMotif{p}(tf,:) =  movmean(mean(signalMat),10);
        end
    else
        out = [];
    end
    
end
end


