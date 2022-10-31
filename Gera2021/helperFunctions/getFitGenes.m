function fitGenes= getFitGenes(sumProm,zth,qth)
    zSc=nanZscore(sumProm);
    fitGenes=find(zSc>min(zth,quantile(zSc,qth)));
end
