function out=mer_occupancy(profile,nmer,varargin)
ip=inputParser;
ip.addParameter('region','promoter');
ip.addParameter('window',21);
ip.addParameter('N', 0);
ip.addParameter('intBases', []);
ip.addParameter('method', 'normal');

ip.parse(varargin{:})
load('SC_genome.mat');
if strcmp(class(profile), 'cell')
    profile = cellfun(@(x)x(:), profile, 'UniformOutput', false);
    profile = [cat(1, profile{:}); nan(85779,1)];
end

% nmersFiles = struct2table(dir('/home/labs/barkailab/tamarge/Master/mat Files/genomeMers/*.mat'));
% nmersFiles.N = str2double(regexp(nmersFiles.name, '(?<=N)\d+', 'match', 'once'));
% nmersFiles.nmer = str2double(regexp(nmersFiles.name, '(?<=nmer)\d+', 'match', 'once'));

if 1==1% any(ismember(nmersFiles.N, ip.Results.N) & ismember(nmersFiles.nmer, nmer))
    %fileIdx = find(ismember(nmersFiles.N, ip.Results.N) & ismember(nmersFiles.nmer, nmer));
    load('./nmer7N0.mat','nmerRed','sc4red')
    %load([nmersFiles.folder{fileIdx}, '/', nmersFiles.name{fileIdx}]);
else
    clear sc4base
    ntBases={'A','C','G','T'};
    for i=1:17
        [~,sc4base{i}]=ismember(upper(SC_genome(i).Sequence'),ntBases);
        sc4base{i}=sc4base{i}-1;
    end
    clear SC_genome
    sc4base=cat(1,sc4base{:});
    kernel=4.^[0:nmer-1];
    %% create mers table
    motifVal=cellfun(@str2num,mat2cell(dec2base([0:4^nmer-1]',4,nmer),ones(4^nmer,1),ones(nmer,1)));
    nmerTable=table(mat2cell(cell2mat(ntBases(motifVal+1)),ones(4^nmer,1),nmer),[0:4^nmer-1]'+1,'VariableNames',{'seq','value'});
    nmerTable.rcValue=sum((3-motifVal).*kernel,2)+1
    nmerTable.rcSeq=nmerTable.seq(nmerTable.rcValue)
    N = ip.Results.N;
    if ip.Results.N>0 & round(nmer/2) == nmer/2
        kernel = [kernel(1:nmer/2), zeros(1,N), kernel(nmer/2+1:end)];
    end
    sc4nmer=conv(sc4base,kernel','same')+1;
    
    nmerRed=nmerTable(nmerTable.value<=nmerTable.rcValue,:);
    sc4red=changem(sc4nmer,nmerRed.value,nmerRed.rcValue);
    sc4red=changem(sc4red,1:size(nmerRed,1),nmerRed.value);
    clear sc4nmer nmerTable sc4base
    % load('sc4red7mer.mat');
    save(sprintf('./nmer%dN%d.mat',nmer,N),  'sc4red', 'nmerRed');
end
GP=load('/home/labs/barkailab/felixj/Documents/MATLAB/scripts/gene/group_imp.mat');
chr_idx=[0 cumsum(GP.chr_len)];
promLen=700;
GP.gene_infoR64.status(cellfun('isempty',GP.gene_infoR64.status))={'Dubious'};
tssData='felixTss';
if strcmpi(ip.Results.region,'promoter')
    gb=false(sum(GP.chr_len),1);
    for i=1:6701
        if all(~isnan(GP.gene_infoR64.(tssData)(i,:)))
            if diff(GP.gene_infoR64.(tssData)(i,2:3))>0
                prom=chr_idx(GP.gene_infoR64.(tssData)(i,1))+[GP.gene_infoR64.(tssData)(i,2)-promLen:GP.gene_infoR64.(tssData)(i,2)];
            else
                prom=chr_idx(GP.gene_infoR64.(tssData)(i,1))+[GP.gene_infoR64.(tssData)(i,2):GP.gene_infoR64.(tssData)(i,2)+promLen];
            end
            gb(prom)=true;
        end
    end
    clear prom
    notGb=false(size(gb));
    for i=find(~contains(GP.gene_infoR64.status,'Dubious')| all(~isnan(GP.gene_infoR64.(tssData)),2))'
        if all(~isnan(GP.gene_infoR64.(tssData)(i,:)))
            cds=chr_idx(GP.gene_infoR64.(tssData)(i,1))+[min(GP.gene_infoR64.(tssData)(i,2:3)):max(GP.gene_infoR64.(tssData)(i,2:3))];
        else
            cds=chr_idx(GP.gene_infoR64.position(i,1))+[min(GP.gene_infoR64.position(i,2:3)):max(GP.gene_infoR64.position(i,2:3))];
        end
        notGb(cds)=true;
    end
else
    gb=true(sum(GP.chr_len),1);
end
gb(notGb)=false;
if ~isempty(ip.Results.intBases)
    gb = logical(ip.Results.intBases);
    if length(gb) < size(sc4red,1)
        gb(size(sc4red,1)) = false;
    end
end
%occ=movmean(profile,2*ip.Results.window+1,1);
score=zeros(size(nmerRed,1),size(profile,2));

for i=1:size(profile,2)
    if strcmp(ip.Results.method, 'normal')
        occ_i=movmean(profile(:,i),ip.Results.window,1);
    else
        %occ_i= max(movmean(profile(:,i),ip.Results.window,1) - movmean(profile(:,i),nmer,1),0);
        occ_i= max(movmean(profile(:,i),ip.Results.window,1) - movmean(profile(:,i),5,1).*5/ip.Results.window ,0);
    end
    score(:,i)=accumarray(sc4red(gb),occ_i(gb),[size(nmerRed,1) 1],@(x)mean(x,'omitnan'));
end
out.mers=nmerRed;
out.mers.occ = accumarray(sc4red(gb),1,[size(nmerRed,1) 1],@(x)sum(x,'omitnan'));
out.score=score;

end


