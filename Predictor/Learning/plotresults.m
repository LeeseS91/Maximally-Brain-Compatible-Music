load('resultsforchorusn1to68.mat')

for i=2:length(predscore)
    meanscore(i-1)=mean(predscore{i}(1,:));
    equalmean(i-1)=mean(predscore{i}(2,:));
end

figure
plot(2:length(predscore),meanscore,'b.')
xlabel('Max Length of Pattern Used', 'fontsize', 14)
ylabel('Score', 'fontsize', 14)
