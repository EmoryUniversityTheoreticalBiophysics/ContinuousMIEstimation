function [] = finchDataExample()
%Runs through some different analyses of the Bengalese Finch data - comparing k and N
%dependence of mutual information estimates, and showing the effects of
%reparameterization.

load('finchData.mat')
%finchData.mat contains two variables, x and y, which 


%First, let's generate some basic plots, showing the dependence of the
%mutual information estimate on N and k
listOfKs = [1 2 3 4 5 7 10 15 20];
listSplitSizes = [1:10];
findMI_KSG_bias_kN(x,y,listOfKs, listSplitSizes, 1, 0, 1);





%Instead, we could make a single plot similar to figure 7 in Holmes &
%Nemenman.
ks = [2 5 20]; %to put it all on one plot, we need a shorter list of ks.
%we frequently find that in higher dimensions k = 1 performs poorly and unusually, and
%so here we'll compare using k = 2 for a low k.

%for each k, we need a list of values of the mutual information and error
%bars for those values.
means = zeros(length(ks), length(listSplitSizes)); 
stds = zeros(length(ks), length(listSplitSizes));

%instead of calling findMI_KSG_bias_kN.m, we can directly call
%findMI_KSG_subsampling.
for i = 1:length(ks)
    [MIs] = findMI_KSG_subsampling(x,y, ks(i), listSplitSizes,0);
    for j = 1:length(listSplitSizes)
        means(i,j) = mean(MIs{j,2});
        stds(i,j) = std(MIs{j,2});
    end
    [errorEstimate] = findMI_KSG_stddev(MIs, length(x), 1);
    stds(i,1) = errorEstimate;
end

figure
hold on
errorbar(listSplitSizes + 0.1,means(1,:),stds(1,:),'b') %the 0.1 offset for the x-axis is just for visualization purposes
errorbar(listSplitSizes,means(2,:),stds(2,:),'g')
errorbar(listSplitSizes - 0.1,means(3,:),stds(3,:),'r') 
xlabel('Number of divisions of the data')
ylabel('Estimated Mutual Information')
legend(horzcat('k = ', num2str(ks(1))), horzcat('k = ', num2str(ks(2))),horzcat('k = ', num2str(ks(3))))
title('N-dependence, unreparameterized data')



%Reparameterize the data to compare the effects of the reparameterization.

x2 = reparameterize_data(x);
y2 = reparameterize_data(y);



listOfKs = [1 2 3 4 5 7 10 15 20];
listSplitSizes = [1:10];
findMI_KSG_bias_kN(x2,y2,listOfKs, listSplitSizes, 1, 0, 1);

%and now we can generate a similar plot, using the reparameterized data
means = zeros(length(ks), length(listSplitSizes));
stds = zeros(length(ks), length(listSplitSizes));

for i = 1:length(ks)
    [MIs] = findMI_KSG_subsampling(x2,y2, ks(i), listSplitSizes,0);
    for j = 1:length(listSplitSizes);
        means(i,j) = mean(MIs{j,2});
        stds(i,j) = std(MIs{j,2});
    end
    [errorEstimate] = findMI_KSG_stddev(MIs, length(x), 1);
    stds(i,1) = errorEstimate;
end

figure
hold on
errorbar(listSplitSizes + 0.1,means(1,:),stds(1,:),'b')
errorbar(listSplitSizes,means(2,:),stds(2,:),'g')
errorbar(listSplitSizes - 0.1,means(3,:),stds(3,:),'r')
xlabel('Number of divisions of the data')
ylabel('Estimated Mutual Information')
legend(horzcat('k = ', num2str(ks(1))), horzcat('k = ', num2str(ks(2))),horzcat('k = ', num2str(ks(3))))
title('N-dependence, reparameterized data')
end
