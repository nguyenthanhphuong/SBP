% Author: Thanh Phuong NGUYEN, U2IS-ENSTA Paristech
%


function [CP fdetection]=ClassificationkNN(LBPH,trainIDs,trainClassIDs,testIDs,testClassIDs,fid,label)

trains = LBPH(trainIDs,:);
tests = LBPH(testIDs,:);
trainNum = size(trains,1);
testNum = size(tests,1);
DM = zeros(testNum,trainNum);
parfor i=1:testNum;
    test = tests(i,:);        
    DM(i,:) = distMATChiSquare(trains,test)';
end
[CP fdetection]=ClassifyOnNN(DM,trainClassIDs,testClassIDs,trainIDs);

fprintf(fid,'======================================\n');
fprintf(fid,'%s\n',label);
fprintf(fid,'CP: %f\n',CP);

fprintf('%s\n',label);
fprintf('CP: %f\n',CP);

