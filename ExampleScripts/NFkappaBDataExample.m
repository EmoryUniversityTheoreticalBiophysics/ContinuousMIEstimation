function [] = NfkBDataAnalysis()
%Runs through some different analyses of the NF-kB data - comparing k and N
%dependence of mutual information estimates, and showing the effects of
%reparameterization.

load('NfkBData.mat')

%[Insert brief description of what NF-kB data actually is]



%First, let's generate some basic plots, showing the dependence of the
%mutual information estimate on N and k
listOfKs = [1 2 3 4 5 7 10 15 20];
listSplitSizes = [1:10];
findMI_KSG_bias_kN(X,Y,listOfKs, listSplitSizes, 1, 0, 1);





%Instead, we could make a single plot similar to figure 7 in Holmes &
%Nemenman.
ks = [1 3 20]; %to put it all on one plot, we need a shorter list of ks.

%for each k, we need a list of values of the mutual information and error
%bars for those values.
means = zeros(length(ks), length(listSplitSizes)); 
stds = zeros(length(ks), length(listSplitSizes));

%instead of calling findMI_KSG_bias_kN.m, we can directly call
%findMI_KSG_subsampling.
for i = 1:length(ks)
    [MIs] = findMI_KSG_subsampling(X,Y, ks(i), listSplitSizes, 0);
    for j = 1:length(listSplitSizes)
        means(i,j) = mean(MIs{j,2});
        stds(i,j) = std(MIs{j,2});
    end
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
ylim([0 1])


%Reparameterize the data to compare the effects of the reparameterization.
[X2] = reparameterize_data(X);
[Y2] = reparameterize_data(Y);



%and now we can generate a similar plot, using the reparameterized data
for i = 1:length(ks)
    [MIs] = findMI_KSG_subsampling(X2,Y2, ks(i), listSplitSizes, 0);
    for j = 1:length(listSplitSizes)
        means(i,j) = mean(MIs{j,2});
        stds(i,j) = std(MIs{j,2});
    end
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
ylim([0 1])