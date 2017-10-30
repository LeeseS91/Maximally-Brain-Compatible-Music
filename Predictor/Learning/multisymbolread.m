% function songreader(file)
clear all
file='symbols';
filename=[file '.txt'];
song=load(filename);
showsong=1;
for N=5 % max length of string examined
    clear y
    NN=4;
    error=0;
    playsong=0;
    
    ptime=0.1; % speed of playback
    nmark=0;
    jjmark=0;
    jjval=0;
    % import sound
    
    % play song
    if playsong==0;
        playsong4(song,ptime)
    end
    
    cell_count = cell(N-1,1);
    cell_pattern = cell(N,2);
    cell_probs = cell(N-1,1);
    for nn=1:N
        %         if nn==1
        %             n2m=[0:NN]';
        %         else
        cell_pattern{nn,1}=de2bi(0:((NN^nn)-1),nn,NN);
        cell_pattern{nn,2}= zeros(length(cell_pattern{nn,1}),1);
        %             for r=1:4^nn
        %                 for c=1:nn
        %                     n2m(r,c)=str2num(pattern(r,c));
        %                 end
        %             end
        %         end
        %         cell_pattern{nn,1}= n2m;
        %         cell_pattern{nn,2}= zeros(length(n2m),1);
        if nn==1
            cell_count{nn}=zeros(1,NN);
            cell_probs{nn}=ones(1,NN)*(1/NN);
        else
            cell_count{nn}=zeros(NN^(nn-1),NN);
            cell_probs{nn}=ones(NN^(nn-1),NN)*(1/NN);
        end
        %         clear n2m
        %         clear pattern
    end
    
    predict(1)=randi(NN)-1;
    equalcount=0;
    
    for j=1:(length(song)-1);
        for n=1:N
            if j>=n
                for jj=1:size(cell_pattern{n,1},1)
                    if song(j-(n-1):j)==cell_pattern{n,1}(jj,:)
                        cell_pattern{n,2}(jj)=cell_pattern{n,2}(jj)+1;
                        jjmark=jj;
                        if rem(jjmark,4^(n-1))==0
                            jjval=4^(n-1);
                            nmark=jjmark/(4^(n-1));
                        else
                            jjval=rem(jjmark,4^(n-1));
                            for ncount=1:4
                                if jjmark-(ncount*(4^(n-1)))<=0
                                    nmark=ncount;
                                    break
                                end
                            end
                            
                        end
                        if n==1
                            cell_count{n}(1,jjmark)=cell_count{n}(1,jjmark)+1;
                        else
                            cell_count{n}(jjval,nmark)=cell_count{n}(jjval,nmark)+1;
                        end
                    end
                end
                
                
                for f=1:size(cell_count{n},1)
                    sumcount=sum(cell_count{n}(f,1:NN));
                    if sumcount~=0
                        cell_probs{n}(f,1:NN)=cell_count{n}(f,1:NN)/sumcount;
                    end
                    clear sumcount
                end
     
                cell_entropy=calcentropy2(cell_probs,N,NN);
                if nmark~=0
                    if n<N
                        %                 if n==1
                        %                     jjmark=1;
                        %                 end
                        probsarray(n,:)=(2^(n-1))*(1-cell_entropy{n+1}(jjmark))*cell_probs{n+1}(jjmark, :);
                        %               (2^(n-1))*
                    end
                    jjmark=0;
                    nmark=0;
                    
                end
            end
        end
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
    
    predict;
    conf=confusionmat(predict, song);
    score(N,1)=1-((conf(1,2)+conf(2,1))/sum(sum(conf)));
    score(N,2)=equalcount;
    
    if showsong==1
        comparmat(1,:)=strrep(num2str(song),' ','');
        comparmat(2,:)=strrep(num2str(predict),' ','');
        
        customshowalign(comparmat);
    end
    
   y(1,:)=strrep(num2str(song),' ','');
if isempty(errorindex)==0
    y(2,errorindex)='|';
    altpeice=predict;
    if isempty(equalindex)==0
        y(3,equalindex)=':';
        y(4,:)=strrep(num2str(altpeice),' ','');
    else
        y(3,:)=strrep(num2str(altpeice),' ','');
    end
end

    %     clear cell*
    %     clear p*
    %     clear n*
    %     clear j*
    %     clear e*
    %     clear c*
end

