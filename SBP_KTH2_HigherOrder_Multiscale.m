% SBP and CSBP model, evaluation on KTH-TIPS 2b dataset
% Author: Thanh Phuong NGUYEN, U2IS-ENSTA Paristech



if matlabpool('size')==0
	matlabpool open 12;
end

%Configuration
rootpic='Datasets/KTH-TIPS2-all/';
classNum=11;
list=dir([rootpic,'*.png']);
picNum=classNum*108*4;
filenames={};
for i=1:picNum
	filenames{i}=sprintf('%s%s', rootpic, list(i).name);
end
R1P1Set={[1 5; 2 8]};
%R1P1Set={[1 4; 2 8], [1 5; 2 6], [1 5; 2 10]};

fid=fopen('Result_KTH_TIP2b_detail_HigherOrder_Multiscale_new.txt','a');
fid1=fopen('Result_KTH_TIP2b_mean_HigherOrder_Multiscale_new.txt','a');

resultsTotal=cell(3,length(R1P1Set),3);
resultsMean=zeros(3,length(R1P1Set),3);
%Construct descriptor
	for rp=1:length(R1P1Set)
	 	
                label=['KTH_TIP2b_'];
		R1=R1P1Set{rp}(:,1); P1=R1P1Set{rp}(:,2);
		labelR1=['_R1']; 
		for j=1:length(R1) 
			labelR1=[labelR1 '_' num2str(R1(j)) ];
		end
		labelP1=['_P1']; 
		for j=1:length(P1) 
			labelP1=[labelP1 '_' num2str(P1(j)) ];
		end
		label=[label labelR1 labelP1];
		%Descriptor

nS=5;
for R=2:nS
	R
	P=16;
	if ~exist(['mapping_riu2_' num2str(P) '.mat'])
		patternMappingriu2 = getmapping(P,'riu2');	
		save(['mapping_riu2_' num2str(P) '.mat'],'patternMappingriu2');
        else
		load(['mapping_riu2_' num2str(P) '.mat'],'patternMappingriu2');	
	end

	
         if ~exist([label num2str(R) num2str(P) '_descr.mat'])
	parfor ii=1:length(filenames)
	  	    filename=filenames{ii};
		    Gray = imread(filename);    
                    Gray=im2double(Gray);
		    Gray = (Gray-mean(Gray(:)))/std(Gray(:))*20+128; % remove global illumination influence


		     [M1, M2, M3, M4]=centralmoment(Gray,R1,P1,1);

	    %CLBP
	    %[M1_S  M1_C] =sbp(M1,R,P,patternMappingriu2,'x');
	    %[M2_S  M2_C]=sbp(M2,R,P,patternMappingriu2,'x');
	    [M1_S M_M1  M1_C] =clbp(M1,R,P,patternMappingriu2,'x');
	    [M2_S M_M2 M2_C]=clbp(M2,R,P,patternMappingriu2,'x');

            %Normalized moment
            B3=M3./(M2.^(3/2));
            B4=M4./(M2.^2);
	    %[B3_S B3_C]=sbp(B3,R,P,patternMappingriu2,'x');
	    %[B4_S B4_C]=sbp(B4,R,P,patternMappingriu2,'x');

	    [B3_S M_B3 B3_C]=clbp(B3,R,P,patternMappingriu2,'x');
	    [B4_S M_B4 B4_C]=clbp(B4,R,P,patternMappingriu2,'x');

 	   

	    SBP_M1 = M1_S;
	    idx = find(M1_C);
	    SBP_M1(idx) = SBP_M1(idx)+patternMappingriu2.num;

            SBP_M2 = M2_S;
	    idx = find(M2_C);
	    SBP_M2(idx) = SBP_M2(idx)+patternMappingriu2.num;


            SBP_B3 = B3_S;
	    idx = find(B3_C);
	    SBP_B3(idx) = SBP_B3(idx)+patternMappingriu2.num;

            SBP_B4 = B4_S;
	    idx = find(B4_C);
	    SBP_B4(idx) = SBP_B4(idx)+patternMappingriu2.num;


	    SBP_M1M2=[SBP_M1(:) SBP_M2(:)];
	    SBP_M1M2=double(SBP_M1M2);
	    Hist3D=hist3(SBP_M1M2,[patternMappingriu2.num*2,patternMappingriu2.num*2]);
	    M1M2H(ii,:)=reshape(Hist3D(:),1,numel(Hist3D));

            SBP_M1B3=[SBP_M1(:) SBP_B3(:)];
	    SBP_M1B3=double(SBP_M1B3);
	    Hist3D=hist3(SBP_M1B3,[patternMappingriu2.num*2,patternMappingriu2.num*2]);
	    M1B3H(ii,:)=reshape(Hist3D(:),1,numel(Hist3D));

	    SBP_M2B3=[SBP_M2(:) SBP_B3(:)];
	    SBP_M2B3=double(SBP_M2B3);
	    Hist3D=hist3(SBP_M2B3,[patternMappingriu2.num*2,patternMappingriu2.num*2]);
	    M2B3H(ii,:)=reshape(Hist3D(:),1,numel(Hist3D));

	    SBP_M1B4=[SBP_M1(:) SBP_B4(:)];
	    SBP_M1B4=double(SBP_M1B4);
	    Hist3D=hist3(SBP_M1B4,[patternMappingriu2.num*2,patternMappingriu2.num*2]);
	    M1B4H(ii,:)=reshape(Hist3D(:),1,numel(Hist3D));

	    SBP_M2B4=[SBP_M2(:) SBP_B4(:)];
	    SBP_M2B4=double(SBP_M2B4);
	    Hist3D=hist3(SBP_M2B4,[patternMappingriu2.num*2,patternMappingriu2.num*2]);
	    M2B4H(ii,:)=reshape(Hist3D(:),1,numel(Hist3D));

	    SBP_B3B4=[SBP_B3(:) SBP_B4(:)];
	    SBP_B3B4=double(SBP_B3B4);
	    Hist3D=hist3(SBP_B3B4,[patternMappingriu2.num*2,patternMappingriu2.num*2]);
	    B3B4H(ii,:)=reshape(Hist3D(:),1,numel(Hist3D));
	
	    M_M1M2=[M_M1(:) M_M2(:)];
	    M_M1M2=double(M_M1M2);
	    Hist3D=hist3(M_M1M2,[patternMappingriu2.num*2,patternMappingriu2.num*2]);
	    M_M1M2H(ii,:)=reshape(Hist3D(:),1,numel(Hist3D));

	    M_M1B3=[M_M1(:) M_B3(:)];
	    M_M1B3=double(M_M1B3);
	    Hist3D=hist3(M_M1B3,[patternMappingriu2.num*2,patternMappingriu2.num*2]);
	    M_M1B3H(ii,:)=reshape(Hist3D(:),1,numel(Hist3D));

	    M_M2B3=[M_M2(:) M_B3(:)];
	    M_M2B3=double(M_M2B3);
	    Hist3D=hist3(M_M2B3,[patternMappingriu2.num*2,patternMappingriu2.num*2]);
	    M_M2B3H(ii,:)=reshape(Hist3D(:),1,numel(Hist3D));

	    M_M1B4=[M_M1(:) M_B4(:)];
	    M_M1B4=double(M_M1B4);
	    Hist3D=hist3(M_M1B4,[patternMappingriu2.num*2,patternMappingriu2.num*2]);
	    M_M1B4H(ii,:)=reshape(Hist3D(:),1,numel(Hist3D));

	    M_M2B4=[M_M2(:) M_B4(:)];
	    M_M2B4=double(M_M2B4);
	    Hist3D=hist3(M_M2B4,[patternMappingriu2.num*2,patternMappingriu2.num*2]);
	    M_M2B4H(ii,:)=reshape(Hist3D(:),1,numel(Hist3D));

	    M_B3B4=[M_B3(:) M_B4(:)];
	    M_B3B4=double(M_B3B4);
	    Hist3D=hist3(M_B3B4,[patternMappingriu2.num*2,patternMappingriu2.num*2]);
	    M_B3B4H(ii,:)=reshape(Hist3D(:),1,numel(Hist3D));

	    M1M2_M1B3_M2B3H(ii,:)=[M1M2H(ii,:) M1B3H(ii,:) M2B3H(ii,:)];
            M1M2_M1B3_M2B3_B3B4_M1B4_M2B4H(ii,:)=[M1M2_M1B3_M2B3H(ii,:) M1B4H(ii,:) M2B4H(ii,:) B3B4H(ii,:)];

	    M_M1M2_M1B3_M2B3H(ii,:)=[M_M1M2H(ii,:) M_M1B3H(ii,:) M_M2B3H(ii,:)];
            M_M1M2_M1B3_M2B3_B3B4_M1B4_M2B4H(ii,:)=[M_M1M2_M1B3_M2B3H(ii,:) M_M1B4H(ii,:) M_M2B4H(ii,:) M_B3B4H(ii,:)];
            		



		end %end parfor

	    save([label num2str(R) num2str(P) '_descr.mat'],'M1M2H','M1M2_M1B3_M2B3H','M1M2_M1B3_M2B3_B3B4_M1B4_M2B4H','M_M1M2H','M_M1M2_M1B3_M2B3H','M_M1M2_M1B3_M2B3_B3B4_M1B4_M2B4H');
     
	clear M1M2H;
        clear M1B3H;
        clear M2B3H;
        clear M1B4H;
        clear M2B4H;
	clear B3B4H;
        clear M1M2_M1B3_M2B3H;
	clear M1M2_M1B3_M2B3_B3B4_M1B4_M2B4H;

	clear M_M1M2H;
        clear M_M1B3H;
        clear M_M2B3H;
        clear M_M1B4H;
        clear M_M2B4H;
        clear M_B3B4H;


        clear M_M1M2_M1B3_M2B3H;
	clear M_M1M2_M1B3_M2B3_B3B4_M1B4_M2B4H;
       end %end if
end %end R

%Classification
SBP2=[];
SBP3=[];
SBP4=[];
CSBP2=[];
CSBP3=[];
CSBP4=[];

labelRP='_RP_';
for R=2:nS
	 P=16;
	 labelRP=[labelRP num2str(R) num2str(P) '-']
         load([label num2str(R) num2str(P) '_descr.mat'],'M1M2H','M1M2_M1B3_M2B3H','M1M2_M1B3_M2B3_B3B4_M1B4_M2B4H','M_M1M2H','M_M1M2_M1B3_M2B3H','M_M1M2_M1B3_M2B3_B3B4_M1B4_M2B4H');
         SBP2=[SBP2 M1M2H];
         SBP3=[SBP3 M1M2_M1B3_M2B3H];
         SBP4=[SBP4 M1M2_M1B3_M2B3_B3B4_M1B4_M2B4H];


         CSBP2=[SBP2 M_M1M2H];
         CSBP3=[SBP3 M_M1M2_M1B3_M2B3H];
         CSBP4=[SBP4 M_M1M2_M1B3_M2B3_B3B4_M1B4_M2B4H];


         clear M1M2H;
         clear M1M2_M1B3_M2B3H;
         clear M1M2_M1B3_M2B3_B3B4_M1B4_M2B4H;

         clear M_M1M2H;
         clear M_M1M2_M1B3_M2B3H;
         clear M_M1M2_M1B3_M2B3_B3B4_M1B4_M2B4H;

	 if(R>=3)
		for i=1:3
			[trainIDSet testIDSet trainClassIDSet testClassIDSet] = constructKTHTIPS2LearningSet(i);
			rsM1M2={};			
			rsM1M2_M1B3_M2B3={};
			rsM1M2_M1B3_M2B3_B3B4_M1B4_M2B4={};

			CrsM1M2={};			
			CrsM1M2_M1B3_M2B3={};
			CrsM1M2_M1B3_M2B3_B3B4_M1B4_M2B4={};
			
		
			for j=1:length(trainIDSet)
				[CP, fd]=ClassificationkNN(SBP2,trainIDSet{j},trainClassIDSet{j},testIDSet{j},testClassIDSet{j},fid,[label labelRP '_NumberTrainGroup_' num2str(i) '_iter_' num2str(j) '_M1M2']);
				rsM1M2{j}=CP;

				[CP, fd]=ClassificationkNN(SBP3,trainIDSet{j},trainClassIDSet{j},testIDSet{j},testClassIDSet{j},fid,[label labelRP '_NumberTrainGroup_' num2str(i) '_iter_' num2str(j) '_rsM1M2_M1B3_M2B3']);
				rsM1M2_M1B3_M2B3{j}=CP;


				[CP, fd]=ClassificationkNN(SBP4,trainIDSet{j},trainClassIDSet{j},testIDSet{j},testClassIDSet{j},fid,[label labelRP '_NumberTrainGroup_' num2str(i) '_iter_' num2str(j) '_M1M2_M1B3_M2B3_B3B4_M1B4_M2B4H']);
				rsM1M2_M1B3_M2B3_B3B4_M1B4_M2B4{j}=CP;


				[CP, fd]=ClassificationkNN(CSBP2,trainIDSet{j},trainClassIDSet{j},testIDSet{j},testClassIDSet{j},fid,[label '_Complement_' labelRP '_NumberTrainGroup_' num2str(i) '_iter_' num2str(j) '_M1M2']);
				CrsM1M2{j}=CP;

				[CP, fd]=ClassificationkNN(CSBP3,trainIDSet{j},trainClassIDSet{j},testIDSet{j},testClassIDSet{j},fid,[label '_Complement_' labelRP '_NumberTrainGroup_' num2str(i) '_iter_' num2str(j) '_rsM1M2_M1B3_M2B3']);
				CrsM1M2_M1B3_M2B3{j}=CP;


				[CP, fd]=ClassificationkNN(CSBP4,trainIDSet{j},trainClassIDSet{j},testIDSet{j},testClassIDSet{j},fid,[label '_Complement_' labelRP '_NumberTrainGroup_' num2str(i) '_iter_' num2str(j) '_M1M2_M1B3_M2B3_B3B4_M1B4_M2B4H']);
				CrsM1M2_M1B3_M2B3_B3B4_M1B4_M2B4{j}=CP;
			end %endj

			mM1M2=mean(cat(1,rsM1M2{:}));
			mM1M2_M1B3_M2B3=mean(cat(1,rsM1M2_M1B3_M2B3{:}));
			mM1M2_M1B3_M2B3_B3B4_M1B4_M2B4=mean(cat(1,rsM1M2_M1B3_M2B3_B3B4_M1B4_M2B4{:}));

			CmM1M2=mean(cat(1,CrsM1M2{:}));
			CmM1M2_M1B3_M2B3=mean(cat(1,CrsM1M2_M1B3_M2B3{:}));
			CmM1M2_M1B3_M2B3_B3B4_M1B4_M2B4=mean(cat(1,CrsM1M2_M1B3_M2B3_B3B4_M1B4_M2B4{:}));

			fprintf(fid,'%s_Multiscale_%s_NumberTrainGroup_%d_M1M2:   %f\n',label,labelRP,i,mM1M2);
			fprintf(fid,'%s_Multiscale_%s_NumberTrainGroup_%d_M1M2_M1B3_M2B3:   %f\n',label,labelRP,i,mM1M2_M1B3_M2B3);
			fprintf(fid,'%s_Multiscale_%s_NumberTrainGroup_%d_M1M2_M1B3_M2B3_B3B4_M1B4_M2B4:   %f\n',label,labelRP,i,mM1M2_M1B3_M2B3_B3B4_M1B4_M2B4);


			fprintf(fid1,'%s_Multiscale_%s_NumberTrainGroup_%d_M1M2:   %f\n',label,labelRP,i,mM1M2);
			fprintf(fid1,'%s_Multiscale_%s_NumberTrainGroup_%d_M1M2_M1B3_M2B3:   %f\n',label,labelRP,i,mM1M2_M1B3_M2B3);
			fprintf(fid1,'%s_Multiscale_%s_NumberTrainGroup_%d_M1M2_M1B3_M2B3_B3B4_M1B4_M2B4:   %f\n',label,labelRP,i,mM1M2_M1B3_M2B3_B3B4_M1B4_M2B4);

			fprintf(fid,'%s_Complement_Multiscale_%s_NumberTrainGroup_%d_M1M2:   %f\n',label,labelRP,i,CmM1M2);
			fprintf(fid,'%s_Complement_Multiscale_%s_NumberTrainGroup_%d_M1M2_M1B3_M2B3:   %f\n',label,labelRP,i,CmM1M2_M1B3_M2B3);
			fprintf(fid,'%s_Complement_Multiscale_%s_NumberTrainGroup_%d_M1M2_M1B3_M2B3_B3B4_M1B4_M2B4:   %f\n',label,labelRP,i,CmM1M2_M1B3_M2B3_B3B4_M1B4_M2B4);


			fprintf(fid1,'%s_Complement_Multiscale_%s_NumberTrainGroup_%d_M1M2:   %f\n',label,labelRP,i,CmM1M2);
			fprintf(fid1,'%s_Complement_Multiscale_%s_NumberTrainGroup_%d_M1M2_M1B3_M2B3:   %f\n',label,labelRP,i,CmM1M2_M1B3_M2B3);
			fprintf(fid1,'%s_Complement_Multiscale_%s_NumberTrainGroup_%d_M1M2_M1B3_M2B3_B3B4_M1B4_M2B4:   %f\n',label,labelRP,i,CmM1M2_M1B3_M2B3_B3B4_M1B4_M2B4);



		end %endi
	end %endif
end %end R

end %end rp




fclose(fid);	
fclose(fid1);

