function [cell_pattern,patcount,probsarray]=datacollectlearn(song,j,n,patcount,cell_pattern,NN,probsarray)



temppattern=song(j-(n-1):j);
if patcount(n)==0 % for intial pattern of length found
    patcount(n)=patcount(n)+1;
    cell_pattern{n,1}(1,:)=temppattern;
    cell_pattern{n,2}(1)=1;
    cell_pattern{n,3}(1,:)=zeros(1,NN);
    cell_pattern{n,4}(1,:)=ones(1,NN)*1/NN;
    latestentry=1;
    patcount(n)=patcount(n)+1;
    
    % now if the pattern already exists in the cell
elseif sum(ismember(cell_pattern{n,1}(:,1:length(temppattern)), temppattern,'rows'))==1
    index=find(ismember(cell_pattern{n,1}(:,1:length(temppattern)), temppattern,'rows')==1);
    cell_pattern{n,2}(index)=cell_pattern{n,2}(index)+1;
    latestentry=index;
    
else % otherwise create a new pattern entry
    cell_pattern{n,1}(patcount(n),:)=temppattern;
    cell_pattern{n,2}(patcount(n))=1;
    cell_pattern{n,3}(patcount(n),:)=zeros(1,NN);
    cell_pattern{n,4}(patcount(n),:)=ones(1,NN)*1/NN;
    cell_pattern{n,5}(patcount(n),:)=1;
    latestentry=patcount(n);
    
    patcount(n)=patcount(n)+1;
end

if n>1
    findindex=find(ismember(cell_pattern{n-1,1}(:,1:length(temppattern)-1), temppattern(1:end-1),'rows')==1);
    cell_pattern{n-1,3}(findindex,cell_pattern{n,1}(latestentry,end))=...
        cell_pattern{n-1,3}(findindex,cell_pattern{n,1}(latestentry,end))+1;
    
    
    sumcount=sum(cell_pattern{n-1,3}(findindex,1:NN));
    if sumcount~=0
        cell_pattern{n-1,4}(findindex,1:NN)=cell_pattern{n-1,3}(findindex,1:NN)/sumcount;
    end
    clear sumcount
    
    %                 cell_pattern{n-1}=calcentropy2(cell_pattern{n,4},N,NN);
    
    cell_pattern{n-1,5}(findindex)=0;
    for iii=1:NN
        cell_pattern{n-1,5}(findindex)=cell_pattern{n-1,5}(findindex)...
            -(((cell_pattern{n-1,4}(findindex,iii)).*log2(cell_pattern{n-1,4}(findindex,iii))/log2(NN)));
    end
    if isnan(cell_pattern{n-1,5}(findindex))==1
        cell_pattern{n-1,5}(findindex)=0;
    end
    
    probindex=find(ismember(cell_pattern{n-1,1}(:,1:length(temppattern)-1), temppattern(end-(n-2):end),'rows')==1);
    probsarray(n,:)=(NN^(n-1))*(1-cell_pattern{n-1,5}(probindex))*cell_pattern{n-1,4}(probindex, :);
    
    
elseif n==1
    prob=cell_pattern{1,2}'/sum(cell_pattern{1,2});
    entrop=0;
    for iii=1:NN
        if prob(iii)~=0
            entrop=entrop-(((prob(iii)).*log2(prob(iii))/log2(NN)));
        end
    end
    probsarray(n,:)=(NN^(n-1))*(1-entrop)*prob;
    clear prob
    clear entrop
end
clear temppattern
clear index
clear findindex

end

