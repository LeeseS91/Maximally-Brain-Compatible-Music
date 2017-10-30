% function songreader(file)
clear all
file='chorus';
filename=[file '.txt'];
song=load(filename);
% songstrct=load('binaryopus1.mat');
% song=songstrct.textformat;
NN=2;
showsong=0;

for N=2:68
    for count=1:5
        error=0;
        % N=63;
        ptime=0.1; % speed of playback
        equalcount=0;
        equalindex=0;
        score=zeros(size(song,1),N);
        
        cell_pattern=cell(N,4);
        cell_pattern{1,1}=[];
        cell_pattern{1,2}=zeros(NN,1);
        cell_pattern{1,3}=zeros(NN,NN);
        cell_pattern{1,4}=ones(NN,NN)*1/NN;
        cell_pattern{1,5}=ones(NN,1);
        patcount=zeros(1,N);
        predict(1)=randi(NN)-1;
        equalcount=0;
        
        for j=1:(length(song)-1);
            for n=1:N
                if j>=n
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
                        cell_pattern{n-1,3}(findindex,cell_pattern{n,1}(latestentry,end)+1)=...
                            cell_pattern{n-1,3}(findindex,cell_pattern{n,1}(latestentry,end)+1)+1;
                        
                        
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
                        
                        probsarray(n,:)=(2^(n-1))*(1-cell_pattern{n-1,5}(probindex))*cell_pattern{n-1,4}(probindex, :);
                    elseif n==1
                        prob=cell_pattern{1,2}'/sum(cell_pattern{1,2});
                        entrop=0;
                        for iii=1:NN
                            if prob(iii)~=0
                                entrop=entrop-(((prob(iii)).*log2(prob(iii))/log2(NN)));
                            end
                        end
                        probsarray(n,:)=(2^(n-1))*(1-entrop)*prob;
                        clear prob
                        clear entrop
                    end
                    
                    clear temppattern
                    clear index
                    clear findindex
                    
                end
            end
            
            if j==1
                probsarray=ones(1,NN)*1/NN;
            end
            p=sort(probsarray,2);
            score(j,1:size(p,1))=p(:,end)';
            if size(probsarray,1)==1
                maxprob=find(probsarray==max(probsarray));
            else
                newprobs=sum(probsarray);
                maxprob=find(newprobs==max(newprobs));
            end
            
            if size(maxprob,2)>1
                choice=randi(size(maxprob,2));
                predict(j+1)=maxprob(choice)-1;
                equalcount=equalcount+1;
                equalindex(equalcount)=j+1;
            else
                predict(j+1)=maxprob-1;
            end
            
            if predict(j+1)~=song(j+1)
                error=error+1;
                errorindex(error)=j+1;
            end
            clear probsindex
            clear probsarray
            clear maxindex
            clear index
            clear maxprob
            clear probmin
            clear probmax
            
        end
        
        
        
        % % predict;
        conf=confusionmat(predict, song);
        predscore{N}(1,count)=1-((conf(1,2)+conf(2,1))/sum(sum(conf)));
        predscore{N}(2,count)=equalcount;
        
        if showsong==1
            comparmat{N}(1,:)=strrep(num2str(song),' ','');
            comparmat{N}(2,:)=strrep(num2str(predict),' ','');
            
            customshowalign(comparmat);
        end
        
%         y{N}(1,:)=strrep(num2str(song),' ','');
%         if isempty(errorindex)==0
%             y{N}(2,errorindex)='|';
%             altpeice=predict;
%             if isempty(equalindex)==0
%                 y{N}(3,equalindex)=':';
%                 y{N}(4,:)=strrep(num2str(altpeice),' ','');
%             else
%                 y{N}(3,:)=strrep(num2str(altpeice),' ','');
%             end
%         end
        clear conf
        clear error*
        clear equal*
        clear choice
        clear cell*
        clear iii
        clear j
        clear lates*
        clear n
        clear predict
        clear pa*
        clear p
        clear score
        clear alt*
        clear new*
        clear prob*
        count
    end
    N
end
