% SBP and CSBP model, evaluation on KTH-TIPS 2b dataset
% Author: Thanh Phuong NGUYEN, U2IS-ENSTA Paristech

clear all

% read picture ID of training and test samples, and read class ID of
% training and test samples

if matlabpool('size')==0
	matlabpool open 12;
end
%Preparation
    R1P1Set={[1 5; 2 8], [1 5; 2 6], [1 6; 2 10], [1 6; 2 12] };
    %R1P1Set={[1 4], [1 6] , [1 8], [1 4; 2 4], [1 4; 2 6],[1 4; 2 8],[1 5; 2 4] [1 5; 2 6],[1 5; 2 8],[1 5; 2 10], [1 6; 2 10],[1 6; 2 12]};
    rootpic10 = 'Datasets/Outex/Outex_TC_00010/';
    picNum10 =4320;
    outfilename10='Result_Outex_TC10_SBP_HighOrder.txt'
    rootpic12 = 'Datasets/Outex/Outex_TC_00012/';
    % picture number of the database
    picNum12 =9120;
    outfilename12='Result_Outex_TC12_SBP_HighOrder.txt'

	trainTxt10 = sprintf('%s/000/train.txt', rootpic10);
	testTxt10 = sprintf('%s/000/test.txt', rootpic10);
	[trainIDs10, trainClassIDs10] = ReadOutexTxt(trainTxt10);
	[testIDs10, testClassIDs10] = ReadOutexTxt(testTxt10);
        for i=1:picNum10
		filenames10{i} = sprintf('%s//images//%06d.ras', rootpic10, i-1);
        end
        for i=1:picNum12
		filenames12{i} = sprintf('%s//images//%06d.ras', rootpic12, i-1);
        end

	
	trainTxt12_1 = sprintf('%s/000/train.txt', rootpic12);
	testTxt12_1 = sprintf('%s/000/test.txt', rootpic12);
	[trainIDs12_1, trainClassIDs12_1] = ReadOutexTxt(trainTxt12_1);
	[testIDs12_1, testClassIDs12_1] = ReadOutexTxt(testTxt12_1);


	trainTxt12_2 = sprintf('%s/001/train.txt', rootpic12);
	testTxt12_2 = sprintf('%s/001/test.txt', rootpic12);
	[trainIDs12_2, trainClassIDs12_2] = ReadOutexTxt(trainTxt12_2);
	[testIDs12_2, testClassIDs12_2] = ReadOutexTxt(testTxt12_2);


        save('param_Outex_SBP.mat','trainIDs10', 'testIDs10', 'trainClassIDs10', 'testClassIDs10', 'trainIDs12_1', 'testIDs12_1', 'trainClassIDs12_1', 'testClassIDs12_1','trainIDs12_2', 'testIDs12_2', 'trainClassIDs12_2', 'testClassIDs12_2','R1P1Set');

for tc=10:2:12
if tc==10
	filenames=filenames10;
	trainIDs={trainIDs10};
	testIDs={testIDs10};
	trainClassIDs={trainClassIDs10};
	testClassIDs={testClassIDs10};
	labelI=['TC10_'];
	fid=fopen(outfilename10,'a');
else
	filenames=filenames12;
	trainIDs={trainIDs12_1,trainIDs12_2};
	testIDs={testIDs12_1,testIDs12_2};
	trainClassIDs={trainClassIDs12_1,trainClassIDs12_2};
	testClassIDs={testClassIDs12_1,testClassIDs12_2};
	labelI=['TC12_'];
	fid=fopen(outfilename12,'a');
end
for R=1:3 
    P=8*R;
    if ~exist(['mapping_riu2_' num2str(P) '.mat'])
        patternMappingriu2 = getmapping(P,'riu2');	
	save(['mapping_riu2_' num2str(P) '.mat'],'patternMappingriu2');
    else
	load(['mapping_riu2_' num2str(P) '.mat'],'patternMappingriu2');	
    end
    %label=[labelI '_R_' num2str(R) '_P_' num2str(P)];
    for rs=1:length(R1P1Set)
	R1=R1P1Set{rs}(:,1); P1=R1P1Set{rs}(:,2);
	labelR1=['_R1']; 
	for j=1:length(R1) 
		labelR1=[labelR1 '_' num2str(R1(j)) ];
        end
	labelP1=['_P1']; 
	for j=1:length(P1) 
		labelP1=[labelP1 '_' num2str(P1(j)) ];
    end
        label=[labelI '_R_' num2str(R) '_P_' num2str(P) labelR1 labelP1];
        %Descriptor
         parfor ii=1:length(filenames)
            filename=filenames{ii};
            Gray = imread(filename);
	    Gray = im2double(Gray);
	    Gray = (Gray-mean(Gray(:)))/std(Gray(:))*20+128; %label image normalization, to remove global intensity

	    [M1, M2, M3, M4]=centralmoment(Gray,R1,P1,1);

	    %CLBP
	    [M1_S  M1_C] =sbp(M1,R,P,patternMappingriu2,'x');
	    [M2_S  M2_C]=sbp(M2,R,P,patternMappingriu2,'x');
	    [M3_S M3_C]=sbp(M3,R,P,patternMappingriu2,'x');
	    [M4_S M4_C]=sbp(M4,R,P,patternMappingriu2,'x');

            %Normalized moment
            B3=M3./(M2.^(3/2));
            B4=M4./(M2.^2);
	    [B3_S B3_C]=sbp(B3,R,P,patternMappingriu2,'x');
	    [B4_S B4_C]=sbp(B4,R,P,patternMappingriu2,'x');

 	   

	    SBP_M1 = M1_S;
	    idx = find(M1_C);
	    SBP_M1(idx) = SBP_M1(idx)+patternMappingriu2.num;

            SBP_M2 = M2_S;
	    idx = find(M2_C);
	    SBP_M2(idx) = SBP_M2(idx)+patternMappingriu2.num;

            SBP_M3 = M3_S;
	    idx = find(M3_C);
	    SBP_M3(idx) = SBP_M3(idx)+patternMappingriu2.num;

            SBP_M4 = M4_S;
	    idx = find(M4_C);
	    SBP_M4(idx) = SBP_M4(idx)+patternMappingriu2.num;

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

	    SBP_M1M3=[SBP_M1(:) SBP_M3(:)];
	    SBP_M1M3=double(SBP_M1M3);
	    Hist3D=hist3(SBP_M1M3,[patternMappingriu2.num*2,patternMappingriu2.num*2]);
	    M1M3H(ii,:)=reshape(Hist3D(:),1,numel(Hist3D));

	    SBP_M2M3=[SBP_M2(:) SBP_M3(:)];
	    SBP_M2M3=double(SBP_M2M3);
	    Hist3D=hist3(SBP_M2M3,[patternMappingriu2.num*2,patternMappingriu2.num*2]);
	    M2M3H(ii,:)=reshape(Hist3D(:),1,numel(Hist3D));

	    SBP_M1M4=[SBP_M1(:) SBP_M4(:)];
	    SBP_M1M4=double(SBP_M1M4);
	    Hist3D=hist3(SBP_M1M4,[patternMappingriu2.num*2,patternMappingriu2.num*2]);
	    M1M4H(ii,:)=reshape(Hist3D(:),1,numel(Hist3D));

	    SBP_M2M4=[SBP_M2(:) SBP_M4(:)];
	    SBP_M2M4=double(SBP_M2M4);
	    Hist3D=hist3(SBP_M2M4,[patternMappingriu2.num*2,patternMappingriu2.num*2]);
	    M2M4H(ii,:)=reshape(Hist3D(:),1,numel(Hist3D));

	    SBP_M3M4=[SBP_M3(:) SBP_M4(:)];
	    SBP_M3M4=double(SBP_M3M4);
	    Hist3D=hist3(SBP_M3M4,[patternMappingriu2.num*2,patternMappingriu2.num*2]);
	    M3M4H(ii,:)=reshape(Hist3D(:),1,numel(Hist3D));


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
	

	    M1M2_M1M3_M2M3H(ii,:)=[M1M2H(ii,:) M1M3H(ii,:) M2M3H(ii,:)];
	    M1M2_M1B3_M2B3H(ii,:)=[M1M2H(ii,:) M1B3H(ii,:) M2B3H(ii,:)];

	    M1M2_M1M3_M2M3_M3M4_M1M4_M2M4H(ii,:)=[M1M2_M1M3_M2M3H(ii,:) M1M4H(ii,:) M2M4H(ii,:) M3M4H(ii,:)];
            M1M2_M1B3_M2B3_B3B4_M1B4_M2B4H(ii,:)=[M1M2_M1B3_M2B3H(ii,:) M1B4H(ii,:) M2B4H(ii,:) B3B4H(ii,:)];
	



	end %endparfor

 	for j=1:length(trainIDs)
		[CP, fd]=ClassificationkNN(M1M2H,trainIDs{j},trainClassIDs{j},testIDs{j},testClassIDs{j},fid,[label 'M1M2_test_' num2str(j)]);
		[CP, fd]=ClassificationkNN(M1M2_M1M3_M2M3H,trainIDs{j},trainClassIDs{j},testIDs{j},testClassIDs{j},fid,[label 'M1M2_M1M3_M2M3_test_' num2str(j)]);
                [CP, fd]=ClassificationkNN(M1M2_M1B3_M2B3H,trainIDs{j},trainClassIDs{j},testIDs{j},testClassIDs{j},fid,[label 'M1M2_M1B3_M2B3_test_' num2str(j)]);

		[CP, fd]=ClassificationkNN(M1M2_M1M3_M2M3_M3M4_M1M4_M2M4H,trainIDs{j},trainClassIDs{j},testIDs{j},testClassIDs{j},fid,[label 'M1M2_M1M3_M2M3_M3M4_M1M4_M2M4H_test_' num2str(j)]);
		[CP, fd]=ClassificationkNN(M1M2_M1B3_M2B3_B3B4_M1B4_M2B4H,trainIDs{j},trainClassIDs{j},testIDs{j},testClassIDs{j},fid,[label 'M1M2_M1B3_M2B3_B3B4_M1B4_M2B4H_test_' num2str(j)]);


	end

        clear M1M2H;
        clear M2M3H;        
        clear M1M3H;
        clear M1M4H;
        clear M2M4H;        
        clear M3M4H;


        clear M1B3H;
        clear M2B3H;
        clear M1B4H;
        clear M2B4H;
        clear B3B4H;


        clear M1M2_M1M3_M2M3H;
        clear M1M2_M1B3_M2B3H;

	clear M1M2_M1M3_M2M3_M3M4_M1M4_M2M4H;
	clear M1M2_M1B3_M2B3_B3B4_M1B4_M2B4H;

    end %endrs
end %endR
	fclose(fid);
end %endtc

