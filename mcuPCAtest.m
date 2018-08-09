% ZStrucTMP=load('E:\DATA\Me\Project\20110512 DevelopFoundOf5Th EMC\Near Scan\Data\Data For FPGA\Result_FPGA_Chipsecurity_0824\9030A_Ez.txt');
clc
close all
tic
% clear;
% global constant For mcu
Ystep=150;Xstep=150;
Ylength=110;Xlength=110;
[XSerialN,YSerialN]=meshgrid(1:1:Ylength,1:1:Xlength); % set XN and YN to 75*75 square
XPosition=XSerialN*Ystep;
YPosition=YSerialN*Xstep;
FreqNum=4001;
FreqCent=500;
FreqSpan=1000;
CountSeri=1:4001;
FreqSeriOrigin=(FreqCent-FreqSpan/2)+CountSeri*(FreqSpan/(FreqNum-1));

% XStrucTMP=load('D:\Study\Me\HongZiyang\Data\NF_FPGA\N9030_1Hx.mat');
% YStrucTMP=load('D:\Study\Me\HongZiyang\Data\NF_FPGA\N9030_1Hy.mat');

% load('F:\Individual Information_Assignment\My Assignments\EN5\Hong\DATA-SCAN\X_98');
% DataSETX=X_98;
% load('F:\Individual Information_Assignment\My Assignments\EN5\Hong\DATA-SCAN\Y_98');
% DataSETY=Y_98;
% clear X_98 Y_98;
% DataSETX=Y98_E; % Ez data
% clear Y98_E X98_E

format longE;
toc
% clear XStrucTMP; clear YStrucTMP;
% ZeroMatrix=zeros(FreqNum-1,Ylength*Xlength);
% FreqMatrixExtract=zeros(FreqNum-1,Xlength,Ylength);

%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%   To extract the image
%        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for FrequencyNum=1:FreqNum
%     %FrequencyNum=915;    %0+80+1-10+100+100;  FrequencyNum(1:1000)
%     %Frequency=FrequencyCent-FrequencySpan/2+FrequencyNum*FrequencySpan/(FrequencyPoin-1);
%     %%show current frequency    
%     for j=0:(Ylength-1)
%     for i=1:Xlength
%         if rem(j,2)==0   %divided with no remainder
% %        FreqMatrixExtractX(FrequencyNum,i,j+1)=DataSETX(Xlength+1-i+j*Xlength,FrequencyNum);
%        FreqMatrixExtractX(FrequencyNum,i,j+1)=DataSETX(Xlength+1-i+j*Xlength,FrequencyNum);
%         else	
% %        FreqMatrixExtractX(FrequencyNum,i,j+1)=DataSETX(i+j*Xlength,FrequencyNum);
%        FreqMatrixExtractX(FrequencyNum,i,j+1)=DataSETX(i+j*Xlength,FrequencyNum);       
%         end
%     end
%     end
%     disp(FrequencyNum);
%     end
%     disp('successfully extracted')
%     toc
%        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% % FreqMatrixExtractHXHY = 10.*log10(10.^(FreqMatrixExtractX./10)+10.^(FreqMatrixExtractY./10));
% 
%        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           %%%%%%%%%%%%%%%%%%%%%%%%%%%%   To remove the ripple in the image
%    for FrequencyNum=1:FreqNum
%     %FrequencyNum=915;    %0+80+1-10+100+100;  FrequencyNum(1:1000)
%     %Frequency=FrequencyCent-FrequencySpan/2+FrequencyNum*FrequencySpan/(FrequencyPoin-1);
%     %%show current frequency    
%     for j=1:(Ylength-1)
%         if rem(j,2)==1   %divided with no remainder
%        FreqMatrixNoRip(FrequencyNum,:,j)=FreqMatrixExtractX(FrequencyNum,:,j); 
%         else
%        FreqMatrixNoRip(FrequencyNum,:,j)=FreqMatrixExtractX(FrequencyNum,:,j-1);
%         end
%     end
%     FreqMatrixNoRip(FrequencyNum,:,Ylength)=FreqMatrixExtractX(FrequencyNum,:,Ylength-1); 
%    end
%     FreqMatrixExtractX=FreqMatrixNoRip;
% 
    FreqMatrixExtractX=FreqMatrixExtractY;
    clear FreqMatrixExtractY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% the following can be deleted, but useful
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Use std to find the pattern
    for FrequencyNum=1:FreqNum
       STDImage(FrequencyNum)=std(std((FreqMatrixExtractX(FrequencyNum,:,:))));
%        SUMImage(FrequencyNum)=sum(sum(FreqMatrixExtractY(FrequencyNum,:,:)));
    end
  
%     subplot(2,1,1);
    STDImage(1)=[];
    FreqSeriOrigin(1)=[];
     figure(10); plot(FreqSeriOrigin,STDImage);
     xlabel('Frequency(MHz)');
    ylabel('Standard Deviation');
%     subplot(2,1,2);loglog(SUMImage(2:(FreqNum-1)));

    Counti=1;
    clear ImagePattern FreqNumPat;
    for FrequencyNum=1:FreqNum-1
        if  STDImage(FrequencyNum)>0.5 % 0.5 can be used for all in our paper | ori:0.74 for Hx, 0.74 for Hy,0.41 for Ez
            ImagePattern(Counti)=FrequencyNum;
            Counti=Counti+1;
        end
    end    

                    
   FreqMatrix=FreqMatrixExtractX;
%    FreqMatrix(1,:,:)=[];
                                          
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%% PCA+KMeans
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FreqNumPat=length(ImagePattern);

% implement the PCA with mean-operation for higher operating speed
DataForShow=zeros(FreqNumPat,Ylength*Xlength);
% DataForPCA=zeros(FreqNumPat,Ylength*Xlength);
 for i=1:FreqNumPat %loading one by one       
      TMPMatrix=reshape(FreqMatrix(ImagePattern(i)+1,:,:),[Ylength Xlength]); %t is trainingset
      DataForShow(i,:)=reshape(TMPMatrix,[1 Ylength*Xlength]);
      
%       % mean-operation
%       for a=1:(Ylength/2)
%           for b=1:(Xlength/2)
%              Data3D(i,a,b)=sum(sum(FreqMatrixNoRip(i,a*2-1:a*2,b*2-1:b*2)))/4;
%           end
%       end
%       TMPMatrix=reshape(Data3D(i,:,:),[Ylength/2 Xlength/2]);
%       DataForPCA(i,:)=reshape(TMPMatrix,[1 Ylength*Xlength/4]);
 end
DataForPCA=DataForShow;
disp('PCA starts')

[COEFF,SCORE,latent]=princomp(DataForPCA); 
%Coeff is principle component,aka basis vector
%score is input's interpretation in coeff basis
%latent is eigenvalues which have been ranked
%vectorMat represents PCA dimensional reduction matrix
coeff=COEFF';
vectorMat=zeros(FreqNumPat-1,Ylength*Xlength);
vectorMat(1,:)=coeff(1,:);
vectorMat(2:(FreqNumPat-1),:)=-1.*coeff(2:(FreqNumPat-1),:);

disp('we got vec...')

FeatureTotal=sum(latent)
PCAFeature=sum(latent(1:5))
PCAFeature/FeatureTotal

immean=(mean(DataForPCA))';

%kmeans preparation
DimentionNum=5; % decomposition pc's number
projected=zeros(FreqNumPat,DimentionNum);%mapped and decrease demension
for i=1:FreqNumPat %test number1001
    DataCentrized=DataForPCA(i,:)'-immean; %data centerize
    projected(i,:)=vectorMat(1:DimentionNum,:)*DataCentrized;%decomposition
end

disp('KMeans...')
%kmeans_solo
[k1,k2]=kSel_subsample(projected);
KmeansGroup=k1;
%projected=whiten(projected);
[idx,centroids]=kmeans(projected,KmeansGroup,'Distance','correlation','Replicates',10); %index,centroids
close all
% figure(10);plot(STDImage(2:(FreqNum)));
for k1=1:KmeansGroup
    ind=find(idx==k1);    
   	figure('NumberTitle','off','Name',strcat('Cluster',num2str(k1)));
%     gray();
    title(['cluster',num2str(k1)]);
        for i=1:min(length(ind),12)
            subplot(3,4,i);
            h=pcolor(YPosition,XPosition,rot90(rot90(reshape(DataForShow(ind(i),:,:),[Ylength Xlength])')));
            set(h,'edgecolor','none');   
            title([num2str((ImagePattern(ind(i)))/4),'MHz']);
            axis('off');
            colormap(jet);       
        end
end