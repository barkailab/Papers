descTF = readtable('/home/labs/barkailab/divyakr/matlab/FinalForThesis/descTF.xlsx');
myCurrTable = descTF;
myCurrTable(cellfun(@isempty,myCurrTable.DBDname),:)=[];

allfamilies = readtable('/home/labs/barkailab/divyakr/matlab/FinalForThesis/YETFASCO_allProteins.xlsx', 'Sheet',  'CISBP + added')       ;
familyList(:,1) = unique(myCurrTable.FAMILY);
k=1;n1=0;n2=0;
for famNo = 1:length(familyList)
    
    if sum( ismember(myCurrTable.FAMILY, familyList(famNo,1)))<=2
       
        n1=n1+sum( ismember(myCurrTable.FAMILY, familyList(famNo,1)));
        n2=n2+sum( ismember(allfamilies.Family_Name, familyList(famNo,1)));
        familyPie(1,2)= {n2}; familyPie(1,3)={n1}; familyPie(1,1)={'Others'};
        continue
    end
    k=k+1;
    familyPie(k,2)={sum( ismember(allfamilies.Family_Name , familyList(famNo,1)))};
    familyPie(k,3)= {sum( ismember(myCurrTable.FAMILY, familyList(famNo,1)))};
    familyPie(k,1)= familyList(famNo,1);
end

familyPie(contains(familyPie(:,1), 'Others'), 2)= {familyPie{contains(familyPie(:,1), 'APSES'), 2}+familyPie{contains(familyPie(:,1), 'Others'), 2}};
familyPie(contains(familyPie(:,1), 'Others'), 3)= {familyPie{contains(familyPie(:,1), 'APSES'), 3}+familyPie{contains(familyPie(:,1), 'Others'), 3}};
familyPie(contains(familyPie(:,1), 'APSES'),:) = [];

familyPieFrac.famName = familyPie(:,1);
familyPieFrac.allTFs = cell2mat(familyPie(:,2))/sum(cell2mat(familyPie(:,2)));
familyPieFrac.myTFs= cell2mat(familyPie(:,3))/sum(cell2mat(familyPie(:,2)));
% familyPieFrac=struct2table(familyPieFrac);
familyPieFrac.allTFsubt = familyPieFrac.allTFs-familyPieFrac.myTFs;
a=[familyPieFrac.allTFsubt' ; familyPieFrac.myTFs'];
figure('color' , [1 1 1], 'Renderer', 'painters' )
labels = [familyPie(:,1) {''; ''; ''; ''; '';''}  ]';
pie(a, labels)
colormap(brewermap(12, 'paired'));
 clearvars familyPieFrac a familyPie n1 n2 k allfamilies familyList  famNo labels