function [data,smeta]=getNaamaData(varargin)
ip=inputParser;
ip.addParameter('figure',1)
ip.addParameter('nucPos',[])
ip.parse(varargin{:});
if ismember(ip.Results.figure,1) 
        allProfiles=struct2table([dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*wt*async*.mat');
            dir('/home/labs/barkailab/felixj/gilad/*wt*async*.mat')
            dir('/home/labs/barkailab/felixj/gilad/glaSpike200/*wt*async*.mat')])
elseif ismember(ip.Results.figure,[1,2]) 
    allProfiles=struct2table([dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*wt*async*.mat');dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*spt6oe*async*.mat')]);
elseif ip.Results.figure==3
    allProfiles=struct2table([
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*fast1a*wt*async*x18.mat');
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*fast1a*wt*h202*x18.mat');
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*h2a*wt*async*x19.mat');
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*h2a*wt*h202*x19.mat')])
elseif ip.Results.figure==4
    allProfiles=struct2table([dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*.mat');dir('/home/labs/barkailab/felixj/gilad/*.mat')]);
    keep=(startsWith(allProfiles.name,{'myc-','ha-','k'})& contains(allProfiles.name,{'-fast','-h2a','-z'}) &contains(allProfiles.name,{'wt','hir1','asf','hir2','rtt','mec','K56','htz1'})&contains(allProfiles.name,{'async','alpha'}));
    keep=keep| (contains(allProfiles.name,{'-x20','-x21'}) & startsWith(allProfiles.name,{'myc-','ha-'}) & contains(allProfiles.name,{'-fast','-h2a'}) )
    allProfiles=allProfiles(keep,:);
elseif ip.Results.figure==5
    allProfiles=struct2table(dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*.mat'));
    keep=contains(allProfiles.name,{'myc-','ha-'})& contains(allProfiles.name,{'fast','h2a'}) &contains(allProfiles.name,{'wt','hir1','asf','hir2','-K56','rtt','mec'})&contains(allProfiles.name,{'async','alpha'});
    allProfiles=allProfiles(keep,:); 
elseif ip.Results.figure==6
    allProfiles=struct2table([dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*fast*wt*x17*.mat');dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*fast*wt*x18*.mat');dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*fast*wt*x23*.mat')])  
elseif ip.Results.figure==7
    allProfiles=struct2table([dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*fast*wt*alpha*x23*.mat')]);
elseif ip.Results.figure==8
    allProfiles=struct2table([dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*wt*async*.mat');dir('/home/labs/barkailab/felixj/gilad/*wt*async*.mat');dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*htz1*async*.mat')]);
elseif ip.Results.figure==9
    allProfiles=struct2table([dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*async*.mat');dir('/home/labs/barkailab/felixj/gilad/*async*.mat')]);
elseif ip.Results.figure==10
    allProfiles=struct2table([dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*-fast*wt*x16*.mat');dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*-fast*wt*alpha*.mat');dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*-fast*wt*async*.mat')])
elseif ip.Results.figure==11
    allProfiles=struct2table([
        %dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*-fast*cc*.mat');
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*-fast*alpha*.mat');
        dir('/hoSEQme/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*-fast*async*.mat')
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*-fast*HU*.mat')])
elseif ip.Results.figure==12
        allProfiles=struct2table([
        dir('/home/labs/barkailab/felixj/gilad/*-fast*async*.mat');
        dir('/home/labs/barkailab/felixj/gilad/*-fast*alpha*.mat');
        dir('/home/labs/barkailab/felixj/gilad/*-h3*async*.mat');
        dir('/home/labs/barkailab/felixj/gilad/*-h3*alpha*.mat');
        dir('/home/labs/barkailab/felixj/gilad/*-h4*async*.mat');
        dir('/home/labs/barkailab/felixj/gilad/*-h4*alpha*.mat');
        dir('/home/labs/barkailab/felixj/gilad/*-hh*alpha*.mat');
        dir('/home/labs/barkailab/felixj/gilad/*-hh*async*.mat');
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*-fast*async*.mat');
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*-fast*alpha*.mat');
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*h4*async*.mat')]);
elseif ip.Results.figure==13
        allProfiles=struct2table([
        dir('/home/labs/barkailab/felixj/gilad/*-fast*async*.mat');
        dir('/home/labs/barkailab/felixj/gilad/*-fast*alpha*.mat');
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*-fast*async*.mat');
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*-fast*alpha*.mat');
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*-fast*HU*.mat');
        dir('/home/labs/barkailab/felixj/gilad/*-h32x*async*.mat');
        dir('/home/labs/barkailab/felixj/gilad/*-h32x*alpha*.mat');
        dir('/home/labs/barkailab/felixj/gilad/*-h32x*HU*.mat')]);
elseif ip.Results.figure==14 %% for analyusis of NaCl timecourse
        allProfiles=struct2table([
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*-fast*wt*async*.mat');
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*-fast*wt*nacl*.mat');
        ]);
elseif ip.Results.figure==15 %%
        allProfiles=struct2table([
        dir('/home/labs/barkailab/felixj/gilad/*async*.mat');
        dir('/home/labs/barkailab/felixj/gilad/*alpha*.mat');
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*async*.mat');
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*alpha*.mat');
        ])
elseif ip.Results.figure==16 %% only async strains
        allProfiles=struct2table([
        dir('/home/labs/barkailab/felixj/gilad/*wt*async*.mat');
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*wt*async*.mat');
        ])
elseif ip.Results.figure==17 %% heatshock timecourse
        allProfiles=struct2table([
        dir('/home/labs/barkailab/felixj/gilad/*t*x33.mat');
        dir('/home/labs/barkailab/felixj/gilad/*t*x34.mat');
        ])
elseif ip.Results.figure==18 %% only async strains
        allProfiles=struct2table([
        dir('/home/labs/barkailab/felixj/gilad/*x33.mat');
        dir('/home/labs/barkailab/felixj/gilad/*x34.mat');        
        dir('/home/labs/barkailab/felixj/gilad/*wt*async*.mat');
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*wt*async*.mat');
        ])
elseif ip.Results.figure==19 %% only async strains
    allProfiles=struct2table([
        dir('/home/labs/barkailab/felixj/gilad/*-fast*wt*.mat');
        dir('/home/labs/barkailab/felixj/gilad/*-h3*wt*.mat');
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*-fast*wt*.mat');
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*-h3*wt*.mat');

        ])
elseif ip.Results.figure==20 %% only async strains
    allProfiles=struct2table([
        dir('/home/labs/barkailab/felixj/gilad/*fast*.mat');
        dir('/home/labs/barkailab/felixj/gilad/*h3*.mat');
        dir('/home/labs/barkailab/felixj/gilad/*hht2*.mat');
        dir('/home/labs/barkailab/felixj/gilad/*h3x2*.mat');
        dir('/home/labs/barkailab/felixj/gilad/*h32x*.mat');
        
        dir('/home/labs/barkailab/felixj/gilad/*200/*fast*.mat');
        dir('/home/labs/barkailab/felixj/gilad/*200/*h3*.mat');
        dir('/home/labs/barkailab/felixj/gilad/*200/*hht2*.mat');
        dir('/home/labs/barkailab/felixj/gilad/*200/*h3x2*.mat');
        dir('/home/labs/barkailab/felixj/gilad/*200/*h32x*.mat');
        
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*h3*.mat');
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*fast*.mat');
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*h3x2*.mat');
        ])
elseif ip.Results.figure==21 %% only async strains
    allProfiles=struct2table([
        dir('/home/labs/barkailab/felixj/gilad/*-fast*async*.mat');
        dir('/home/labs/barkailab/felixj/gilad/*-fast*alpha*.mat');
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*-fast*alpha*.mat');
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*-fast*async*.mat');
        ])
elseif ip.Results.figure==22 %% look at spike in experiments
    allProfiles=struct2table([
        dir('/home/labs/barkailab/felixj/gilad/*rlf*-x36*.mat');
        dir('/home/labs/barkailab/felixj/gilad/*/*rlf*-x36*.mat');
        dir('/home/labs/barkailab/felixj/gilad/*para*-x36*.mat');
        dir('/home/labs/barkailab/felixj/gilad/*/*para*-x36*.mat');
        ])
elseif ip.Results.figure==23 %% look cell cycke h3x2
    allProfiles=struct2table([
        dir('/home/labs/barkailab/felixj/gilad/*h3x2*.mat');
        dir('/home/labs/barkailab/felixj/gilad/*h32x*.mat');
        ])
elseif  ip.Results.figure==24
    allProfiles=struct2table([
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*-fast*.mat');
        dir('/home/labs/barkailab/felixj/gilad/*-fast*.mat');
        dir('/home/labs/barkailab/felixj/gilad/*h32x*.mat');
        dir('/home/labs/barkailab/felixj/gilad/*h3x2*.mat');
        ])    
elseif  ip.Results.figure==25
    allProfiles=struct2table([
        dir('/home/labs/barkailab/felixj/gilad/glaSpike200/*.mat')]);
elseif ip.Results.figure==26
    allProfiles=struct2table([
        dir('/home/labs/barkailab/felixj/gilad/glaSpike200/*h3x2*.mat');    
        %dir('/home/labs/barkailab/felixj/gilad/glaSpike200/*22*.mat');    
        %dir('/home/labs/barkailab/felixj/gilad/*h32x*.mat');    
        dir('/home/labs/barkailab/felixj/gilad/*h3x2*.mat')]);    
        %dir('/home/labs/barkailab/felixj/gilad/yoav/*.mat')]);    
elseif ip.Results.figure==27 % anchor away strains
    allProfiles=struct2table([
        dir('/home/labs/barkailab/felixj/gilad/glaSpike200/*hht2*.mat');    
        dir('/home/labs/barkailab/felixj/gilad/*hht2*.mat');
        dir('/home/labs/barkailab/felixj/gilad/*h2b*.mat');
        dir('/home/labs/barkailab/felixj/gilad/glaSpike200/*h2b*.mat');    
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*h2b*.mat');
        ]);
elseif ip.Results.figure==28
    allProfiles=struct2table([   
        dir('/home/labs/barkailab/felixj/gilad/yoav/*.mat')]);    
elseif ip.Results.figure==29
    allProfiles=struct2table([   
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*h1*.mat')
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*hho1*.mat')
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*h2*wt-async*.mat')
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*fast*wt-async*.mat')
        dir('/home/labs/barkailab/felixj/gilad/*h2*wt-async*.mat')
        dir('/home/labs/barkailab/felixj/gilad/*fast*wt-async*.mat')
        dir('/home/labs/barkailab/felixj/gilad/glaSpike200/*d1.mat')
                ]);
elseif ip.Results.figure==30
    allProfiles=struct2table([   
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*wt*async*.mat')
        dir('/home/labs/barkailab/felixj/gilad/*wt*async*.mat')
        dir('/home/labs/barkailab/felixj/gilad/glaSpike200/*wt*async*.mat')
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*wt*alpha*.mat')
        dir('/home/labs/barkailab/felixj/gilad/*wt*alpha*.mat')
        dir('/home/labs/barkailab/felixj/gilad/glaSpike200/*wt*alpha*.mat')
        ]);
elseif ip.Results.figure==31
    allProfiles=struct2table([
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*hht2*.mat')
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*wt*x18*.mat')
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/k79*.mat')
        dir('/home/labs/barkailab/LAB/data/SEQ/gilad_20180313_histones/MatLB/data/*hht2*.mat')
        dir('/home/labs/barkailab/felixj/gilad/*hht2*.mat')
        dir('/home/labs/barkailab/felixj/gilad/glaSpike200/*hht2*.mat')
        ]);
elseif ip.Results.figure==32
    allProfiles=struct2table([
        dir('/home/labs/barkailab/felixj/gilad/glaSpikeFS/*.mat')
        ]);
end
allProfiles=unique(allProfiles)
expDesc=regexp(allProfiles.name,'-','split');
expDesc=cat(1,expDesc{:});
allProfiles.ab=expDesc(:,1);
allProfiles.tag=expDesc(:,2);
allProfiles.gt=expDesc(:,3);
allProfiles.con=expDesc(:,4);
allProfiles.exp=expDesc(:,5);
[allTags,~,allProfiles.tagid]=unique(allProfiles.tag);
if ismember(ip.Results.figure,[1,2]) 
    intTags=find(true(size(allTags)));%find(~contains(allTags,{'h1','z','sup'},'IgnoreCase',true))
    intProfiles=ismember(allProfiles.ab,{'myc','ha','myc05','myc25'})  & ismember(allProfiles.tagid,intTags);
    %intProfiles=ismember(allProfiles.ab,{'myc','ha'}) & ismember(allProfiles.tagid,[1:4,15:35,39:40]);
elseif ip.Results.figure==2
    intProfiles=ismember(allProfiles.ab,{'myc','ha'}) & ismember(allProfiles.tagid,[1:4,15:35,39:40]);    
elseif ip.Results.figure==3
    intProfiles=true(size(allProfiles,1),1);
elseif ip.Results.figure==4
    intProfiles=ismember(allProfiles.tagid,[1:9,15:18]) & ~contains(allProfiles.con,'h2o2') &contains(allProfiles.ab,{'myc','ha','k'})
elseif ip.Results.figure==5
    intProfiles=allProfiles.tagid<10 & ~contains(allProfiles.con,'h2o2') &contains(allProfiles.ab,{'myc','ha'})
elseif ip.Results.figure==6      
    intProfiles=allProfiles.tagid<5 & contains(allProfiles.name,'wt') & ~contains(allProfiles.name,'HU') &contains(allProfiles.ab,{'myc','ha'})
elseif ip.Results.figure==7
    intProfiles=true(size(allProfiles,1),1);
elseif ip.Results.figure==8
    intProfiles=contains(allProfiles.ab,{'myc','ha'}) & ~contains(allProfiles.tag,'sup');    
elseif ip.Results.figure==9
    intProfiles=contains(allProfiles.ab,{'myc','ha'}) & ~contains(allProfiles.tag,'sup')& contains(allProfiles.tag,{'fast','h2a'},'IgnoreCase',true);    
elseif ip.Results.figure==10
    intProfiles=contains(allProfiles.ab,{'myc','ha'}) & ~contains(allProfiles.con,{'h2o2'});
elseif ip.Results.figure==11
    intProfiles=contains(allProfiles.ab,{'myc','ha'}) & ~contains(allProfiles.con,{'h2o2'}) & contains(allProfiles.gt,{'wt','rtt','mec1','asf'});
elseif ip.Results.figure==12
    intProfiles=contains(allProfiles.gt,{'wt','asf','hir','spt21','k27'}) &~contains(allProfiles.con,'h2o2') & contains(allProfiles.ab,{'myc','h3','ha'});
elseif ip.Results.figure==13
    intProfiles=contains(allProfiles.gt,{'wt','asf','rtt','mec'}) &~contains(allProfiles.con,'h2o2') & contains(allProfiles.ab,{'myc','ha'});
elseif ip.Results.figure==14
    intProfiles=contains(allProfiles.ab,{'HA','myc'},'IgnoreCase',true);
elseif ip.Results.figure==15
    intProfiles=contains(allProfiles.gt,{'wt','asf','hir','rlf'},'IgnoreCase',true)&contains(allProfiles.ab,{'HA','myc','h3'},'IgnoreCase',true)&~contains(allProfiles.tag,{'h1','sup','z'},'IgnoreCase',true)&~contains(allProfiles.con,{'h2o2'},'IgnoreCase',true);
elseif ip.Results.figure==16
    intProfiles=contains(allProfiles.gt,{'wt','asf','hir','rlf'},'IgnoreCase',true)&contains(allProfiles.ab,{'HA','myc','h3'},'IgnoreCase',true)&~contains(allProfiles.tag,{'h1','sup','z'},'IgnoreCase',true)&~contains(allProfiles.con,{'h2o2'},'IgnoreCase',true);
elseif ip.Results.figure==17
    intProfiles=contains(allProfiles.gt,{'wt','rpb1ts'},'IgnoreCase',true)&contains(allProfiles.tag,{'h3'},'IgnoreCase',true);
elseif ip.Results.figure==18
    intProfiles=contains(allProfiles.tag,{'fast','h3','notag','nc1'},'IgnoreCase',true)&contains(allProfiles.ab,{'myc','ha'},'IgnoreCase',true);    
elseif ip.Results.figure==19
    intProfiles=~contains(allProfiles.con,{'cc','HU'},'IgnoreCase',true)&contains(allProfiles.ab,{'myc','ha'},'IgnoreCase',true);
elseif ip.Results.figure==20
    intProfiles=~contains(allProfiles.tag,{'tev'},'IgnoreCase',true)&contains(allProfiles.ab,{'myc','ha','h3'})&(cellfun('prodofsize',regexp(allProfiles.tag,'^hht2|^h3x2|^h32x','once')) | ...
        (contains(allProfiles.con,{'async','alpha','cc'},'IgnoreCase',true)&~contains(allProfiles.con,{'h2o2'},'IgnoreCase',true) &cellfun('prodofsize',regexp(allProfiles.tag,'^fast|^h3','once'))&contains(allProfiles.gt,{'wt','rtt'})));
elseif ip.Results.figure==21
    intProfiles=contains(allProfiles.ab,{'myc','ha'}) &~contains(allProfiles.con,{'h2o2'})
elseif ip.Results.figure==22
    intProfiles=true(size(allProfiles.ab));
elseif ip.Results.figure==23
   stand intProfiles=contains(allProfiles.gt,{'wt','rtt','sml'})&contains(allProfiles.ab,{'ha','myc'})&contains(allProfiles.exp,{'31','35','36'})&~contains(allProfiles.ab,{'f'});%true(size(allProfiles.ab));
elseif ip.Results.figure==24
   intProfiles= contains(allProfiles.con,{'async','alpha','CC','HU'},'IgnoreCase',true)&~contains(allProfiles.con,{'h2o2'},'IgnoreCase',true)&contains(allProfiles.ab,{'myc','ha'},'IgnoreCase',true)&contains(allProfiles.gt,{'asf','mec','mec1','mecsml','rtt','wt'});
elseif ip.Results.figure==25
   %intProfiles=contains(allProfiles.con,'spike','IgnoreCase',true)&contains(allProfiles.tag,'h3x2','IgnoreCase',true);%contains(all.Profiles.tag)
   intProfiles=contains(allProfiles.con,'spike','IgnoreCase',true)%&contains(allProfiles.tag,'h3x2','IgnoreCase',true);%contains(all.Profiles.tag)
   %intProfiles= contains(allProfiles.con,{'alpha','async'},'IgnoreCase',true);%~cellfun('isempty',regexp(allProfiles.ab,'^(myc|ha)$'))&
elseif ip.Results.figure==26
    intProfiles=(startsWith(allProfiles.name,{'myc-','ha-','k56ac','k9ac','DNA'})&contains(allProfiles.gt,{'wt','rtt','mec','hst4','vps'})) | contains(allProfiles.exp,'k1')&contains(allProfiles.gt,'WT','IgnoreCase',true);
elseif ip.Results.figure==27
    intProfiles=startsWith(allProfiles.name,{'ha-','myc-','h3-','h2b-'})&~contains(allProfiles.tag,{'nc','int'},'IgnoreCase',true)&~contains(allProfiles.con,'alpha')&(contains(allProfiles.exp,'s','ignoreCase',true)|contains(allProfiles.con,'spike','ignoreCase',true));
elseif ip.Results.figure==28
    intProfiles=(startsWith(allProfiles.name,{'myc-','ha-'})&contains(allProfiles.gt,{'wt','rtt','mec'}))|contains(allProfiles.name,'-k2');
elseif ip.Results.figure==29
    intProfiles=(startsWith(allProfiles.name,{'myc-','ha-'})&~contains(allProfiles.tag,{'inth2','sup'}))|contains(allProfiles.name,'-d1');    
    
elseif ip.Results.figure==30
    intProfiles=contains(allProfiles.tag,{'fast','h3x2'})&~contains(allProfiles.tag,{'sup'})&contains(allProfiles.exp,{'x','s','X','S'})&contains(allProfiles.ab,{'myc','ha','h3'});
elseif ip.Results.figure==31
    intProfiles=~contains(allProfiles.exp,'f')&~contains(allProfiles.con,'HU')&~contains(allProfiles.name,'mTMPyc')   ;
elseif ip.Results.figure==32
    intProfiles=contains(allProfiles.name,{'h3x2','edu'})%~contains(allProfiles.exp,'f')&~contains(allProfiles.con,'HU')&~contains(allProfiles.name,'mTMPyc')   ;
end


allProfiles=allProfiles(intProfiles,:);
nProfiles=sum(intProfiles);
GP=load('group_imp.mat');
if numel(ip.Results.nucPos)>0
    data=zeros(max(nucPos),nProfiles);
else    
    data=zeros(sum(GP.chr_len),nProfiles);
end
totReads=zeros(nProfiles,1);
for i=1:nProfiles
    temp=load([allProfiles.folder{i} '/' allProfiles.name{i}]);
    if numel(ip.Results.nucPos)>0
        data(:,i)=accumarray(ip.Results.nucPos(ip.Results.nucPos>0),data(ip.Results.nucPos>0,i),[],@mean);
    else
        data(:,i)=temp.profile.data;
    end   
    if isfield(temp.profile,'readCount')
        allProfiles.folder{i}=temp.profile.readCount;
    end
    totReads(i)=temp.profile.SumReads;
end
clear temp
allProfiles.totReads=totReads;
smeta=allProfiles;
smeta.tag2(contains(smeta.tag,'fast'))={'h3'};
smeta.tag2(contains(smeta.tag,'h3')&~contains(smeta.tag,'h32x')&~contains(smeta.tag,'h3tev'))={'h3clean'};
smeta.tag2(contains(smeta.tag,'nc'))={'nc_h3'};
smeta.tag2(contains(smeta.tag,'h2a'))={'h2a'};
smeta.tag2(contains(smeta.tag,'h2b')&~contains(smeta.tag,'int'))={'h2b'};
smeta.tag2(contains(smeta.tag,'h4'))={'h4'};
smeta.tag2(contains(smeta.tag,'z') & ~contains(smeta.tag,'znc'))={'h2z'};
smeta.tag2(contains(smeta.tag,'znc'))={'nc_h2z'};
smeta.tag2(contains(smeta.tag,'notag'))={'none'};
smeta.tag2(contains(smeta.tag,'h1NCon'))={'nc_h1'};
smeta.tag2(contains(smeta.tag,'h1on'))={'h1'};
smeta.tag2(contains(smeta.tag,'h32x'))={'h3y2'};
smeta.tag2(contains(smeta.tag,'h3x2'))={'h3x2'};

smeta.tag2(contains(smeta.tag,'hht2'))={'hht2'};
smeta.tag2(contains(smeta.tag,'hhf1'))={'hhf1'};

smeta.tag2(contains(smeta.tag,'h2bNC'))={'nc_h2b'};
smeta.tag2(contains(smeta.tag,'h2anc'))={'nc_h2a'};
smeta.tag2(contains(smeta.tag,'slow'))={'slow_h3'};
smeta.tag2(contains(smeta.tag,'h3tev'))={'h3Tev2'};
smeta.tag2(contains(smeta.tag,'h2b')&contains(smeta.tag,'int'))={'h2Int'};
smeta.tag2(contains(smeta.tag,'paradoxus'))={'parNC'};
smeta.tag2(contains(smeta.tag,'glabrata'))={'glaNC'};
smeta.tag2(contains(smeta.tag,'DNA'))={'DNA'};

%% highlight bad samples
badList=readtable('badSamples.txt','ReadVariableNames',false);


smeta.bad=contains(smeta.name,badList.Var1);

end