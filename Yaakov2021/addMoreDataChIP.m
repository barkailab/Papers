%%% add additional profile to the data matrix %%%%
clear all
Exp_folder{14} =  '/home/labs/barkailab/LAB/data/RAWSEQ/20190502_gilad_michal_offir/190502_NB551168_0317_AH5NNLBGX9/out' ; 
Exp_folder{15} =  '/home/labs/barkailab/LAB/data/RAWSEQ/20190502_gilad_michal_offir/190502_NB551168_0317_AH5NNLBGX9/out' ; 
Exp_folder{16} =  '/home/labs/barkailab/LAB/data/RAWSEQ/20190519_Gilad/190519_NB501465_0502_AH5NM2BGX9/out' ; 
Exp_folder{17} =  '/home/labs/barkailab/LAB/data/RAWSEQ/20190529_gilad/190529_NB551168_0324_AHGLV7BGXB/out' ; 
Exp_folder{18} =  '/home/labs/barkailab/LAB/data/RAWSEQ/20190616_Gilad/190616_NB551168_0331_AHGKM7BGXB/out' ; 
Exp_folder{19} =  '/home/labs/barkailab/LAB/data/RAWSEQ/20190711_Gilad_Michal_Nelly/out' ; 
Exp_folder{20} =  '/home/labs/barkailab/LAB/data/RAWSEQ/20190829_GiladMiri/out' ; 
Exp_folder{21} =  '/home/labs/barkailab/LAB/data/RAWSEQ//20190904_gilad/190904_A00929_0006_AHL3VHDMXX/out' ; 
Exp_folder{22} =  '/home/labs/barkailab/LAB/data/RAWSEQ/20191113_NovaSeq/out'
Exp_folder{23} =  '/home/labs/barkailab/LAB/data/RAWSEQ/20200101_Gilad_Sagie_Offir/out'
Exp_folder{24} =  '/home/labs/barkailab/LAB/data/RAWSEQ/20200109_gilad_felix/out'
Exp_folder{25} =  '/home/labs/barkailab/LAB/data/RAWSEQ/20200216_GiladGatDivyaTamarG/out_bcs11/'; 
Exp_folder{27} =  '/home/labs/barkailab/LAB/data/RAWSEQ/20200907_GiladTamarFelix/outm0/chip/std/'; 
Exp_folder{28} =  '/home/labs/barkailab/LAB/data/SEQ/Gilad_X28_200923/std'; 
Exp_folder{30} =  '/home/labs/barkailab/LAB/data/SEQ/Gilad_X30_20201108/naama'; 
Exp_folder{31} =  '/home/labs/barkailab/LAB/data/SEQ/Gilad_X31_20201127/naama'; 
Exp_folder{32} =  '/home/labs/barkailab/LAB/data/SEQ/Gilad_X32_20210131/naama'; 
Exp_folder{34} =  '/home/labs/barkailab/LAB/data/SEQ/Gilad_X33X34_210218/naama';
Exp_folder{35} =  '/home/labs/barkailab/LAB/data/SEQ/Gilad_X35_20210322/naama';
Exp_folder{36} = '/home/labs/barkailab/LAB/data/SEQ/Gilad_X36_20210426/naama';
Exp_folder{36} = '/home/labs/barkailab/LAB/data/SEQ/Gilad_X36Spike_20210426/spikeInPom';
RmInd = [7695000:7750000] ; 
% nNuc{1} = '/*_*_*_*220.out' ; nNuc{2} = '/*_*_*_*124.out' ;nNuc{3} = '/*_*_*_*89.out' ; nNuc{4} = '/*_*_*_*500.out'
% SDir{1} = 'data_Nuc220' ;  SDir{2} = 'data_Nuc124' ;SDir{3} = 'data_Nuc89'; SDir{4} ='data_NucAll' ; 

nNuc{1} = '/*.out'
SDir{1} ='data' ;

load '/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/GeneralPara/GPNaama.mat' ;
global GPNaama
innerBCs=readtable('/home/labs/barkailab/LAB/scripts/felix_bcl2fastq/dna_barcodes.txt','ReadVariableNames',false);
innerBCs=innerBCs(:,[5 6]);
outerBCs=table({'Enr1','Enr2','Enr3','Enr4'}',{'TCTACTCT','CTCCTTAC','TATGCAGT','TACTCCTT'}')
SaveDir = '~/gilad/1pike/';
sd=1


for num_exp =  36
    XLSName = [Exp_folder{num_exp},'/X',int2str(num_exp),'_DNA_index.xlsx']
    
    [num,NameIndex_list] = xlsread(XLSName);
    [~,bcOutIdx]=ismember(regexp(NameIndex_list(:,1),'^[ACGT]{8}','match','once'),outerBCs.Var2)
    [~,bcInIdx]=ismember(regexp(NameIndex_list(:,1),'[ACGT]{8}$','match','once'),innerBCs.Var6)
    NameIndex_list(bcOutIdx>0 & bcInIdx>0,1)=strcat(outerBCs.Var1(bcOutIdx(bcOutIdx>0 & bcInIdx>0)),'_',innerBCs.Var5(bcInIdx(bcOutIdx>0 & bcInIdx>0)))
    
    %PlateInd = cellfun(@(x) char(x) ,NameIndex_list(:,1),'UniformOutput',false)' ;
    
    DirName = Exp_folder{num_exp} ;
    ProfileToAdd = dir([DirName,nNuc{sd}]) ;
    if size(NameIndex_list)>2
        NameIndex_list=NameIndex_list(:,[1,3]);
    end
    parfor j=1:size(ProfileToAdd,1)
        try
            file_name = [DirName,'/',ProfileToAdd(j).name]
            %ProfID = regexp(ProfileToAdd(j).name,'(\w*_\w*)_\w*_\w*','tokens') ; ProfID = char(ProfID{1});
            ProfID = regexp(ProfileToAdd(j).name,'Enr\d_\d{1,2}[A-H]','match') ;
            l = find(contains(NameIndex_list(:,1), ProfID))            
            if (l) 
                saveFileName= [SaveDir, '/' NameIndex_list{l,2} '.mat']
                if ~exist(saveFileName,'file')
                    profile=[];
                    fid = fopen(file_name) ;
                    profile.data = fscanf(fid, '%g\n');  profile.data(RmInd) = 0 ;
                    profile.SumReads = sum(profile.data)/100;
                    profile.data  =  profile.data/sum(profile.data)* 10^6 *100  ;
                    profile.name = NameIndex_list{l,2}
                    profile.full_name = file_name ;
                    fclose(fid) ;
                    saveProfile(profile, SaveDir) ;
                end
            end
            %system(['rm ', file_name]) ;
            bad(j)=false;
        catch
            bad(j)=true;
        end
        
    end
    find(bad)
end


