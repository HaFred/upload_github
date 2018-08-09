function [k1,k2]=kSel_subsample(projected)
krange=[2:7];
prop=0.8;
reps=200;
FreqMatrixX=projected;
% clear S mS coeff COEFF DataForPCA vectorMat

% FreqMatrix is MxN matrix where M is number of observations, N is number of
% features
num_samps=round(prop*size(FreqMatrixX,1));
S=zeros(length(krange),reps);
R=zeros(length(krange),reps);

% generate 2 subsamples
for k=krange
    for m=1:reps
    % generating sub-sample indices    
    r1=randperm(size(FreqMatrixX,1))';
    ss1=r1(1:num_samps,1);
    r2=randperm(size(FreqMatrixX,1))';
    ss2=r2(1:num_samps,1);
    % performing Kmeans on subsamples
    [idx1]=kmeans(FreqMatrixX(ss1,:),k,'Distance','correlation','Replicates',10);
    [idx2]=kmeans(FreqMatrixX(ss2,:),k,'Distance','correlation','Replicates',10);
    % find set intersection, cluster distance using Rand Index or Jaccard
    % Jaccard index
    [C,iss1,iss2]=intersect(ss1,ss2);
    D1=idx1(iss1,1);
    D2=idx2(iss2,1);    
    % partition similarity
    S(k-1,m)=partsim(D1,D2);
    
    % rand Index
    [R(k-1,m),output]=randIndex(D1,D2);
%     if output==1
%         disp('function randIndex works...');  
%     else
%         disp('A error in randIndex');
%     end
    end 
end
mS=mean(S,2)';
mR=mean(R,2)';
k1=krange(find(mS==max(mS),1,'first'));
k2=krange(find(mR==max(mR),1,'first'));
display(k1)
display(k2)
end
