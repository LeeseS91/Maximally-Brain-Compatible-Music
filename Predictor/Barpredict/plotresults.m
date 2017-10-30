load('barscores.mat')

for i=1:length(predscore)
    meanscore(i)=predscore{i}(1);
    equalmean(i)=predscore{i}(2);
end

figure
plot(1:length(predscore),meanscore,'k.','Markersize',15)
xlabel('Bar length used', 'fontsize', 14)
ylabel('Score', 'fontsize', 14)
