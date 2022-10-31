function [seq, relBases] = findPromoterSeq(gene, GP, promoterLengthVec, varargin)

load('SC_genome.mat');

if ischar(gene)
    geneIdx = GP.gene_table.(upper(gene));
else
    geneIdx = gene;
end

ip = inputParser;
ip.addParameter('tss', {'Stein'});
ip.addParameter('intoGene', 0);

ip.parse(varargin{:});
tss = ip.Results.tss;
intoGene =  ip.Results.intoGene;

promL = promoterLengthVec(geneIdx);
L = cumsum([1,GP.chr_len]);
gdir = sign(diff(GP.gene_infoR64.position(geneIdx,[2,3]))) ;

if isnan(promL)
    seq = [];
else
    if strcmp(tss, 'Stein')
        TssRelPos = GP.gene_infoR64.stein_tss(geneIdx,:);
    elseif strcmpi(tss, 'tamar')
        TssRelPos = GP.gene_infoR64.tamarTss(geneIdx,:);
    elseif strcmpi(tss, 'ORF')
        TssRelPos = GP.gene_infoR64.position(geneIdx,:);
    end
    
    if gdir > 0
        relBases = TssRelPos(2)-promL: TssRelPos(2)-1+intoGene;
        seq = SC_genome(TssRelPos(1)).Sequence(TssRelPos(2)-promL: TssRelPos(2)-1+intoGene);
    else
        relBases = flipud(TssRelPos(2)+1-intoGene: TssRelPos(2)+promL);
        seqTemp = SC_genome(TssRelPos(1)).Sequence(TssRelPos(2)+1-intoGene: TssRelPos(2)+promL);
        seq = seqrcomplement(seqTemp);
    end
end
end
