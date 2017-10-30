% function songreader(file)
clear all
% file='chorus';
% filename=[file '.txt'];
% song=load(filename);
songstrct=load('binaryopus1.mat');
song=songstrct.textformat;
NN=length(unique(song));
for barsize=4:4:32
    for N=4:4:32
        for M=4:16
            % NN=2;
            % barsize=8;
            
            error=0;
            % N=10;
            % M=10;
            ptime=0.1; % speed of playback
            % nmark=0;
            % jjmark=0;
            equalcount=0;
            equalindex=0;
            score=zeros(size(song,1),N);
            
            cell_pattern=cell(N,4);
            cell_pattern{1,1}=(0:NN-1)';
            cell_pattern{1,2}=zeros(NN,1);
            cell_pattern{1,3}=zeros(NN,NN);
            cell_pattern{1,4}=ones(NN,NN)*0.5;
            cell_pattern{1,5}=ones(NN,1);
            patcount=zeros(1,N);
            predict(1)=randi(NN)-1;
            barpredict(1)=1;
            equalcount=0;
            equalcountb=0;
            barerror=0;
            
            barcount=1;
            barpatcount=zeros(1,M);
            ii=0;
            add=0;
            barcell_pattern=cell(M,4);
            % barcell_pattern{1,1}=(0:NN-1)';
            % barcell_pattern{1,2}=zeros(NN,1);
            
            for j=1:(length(song)-1);
                
                if j>=barsize
                    if rem(j,barsize)==0;
                        ii=ii+1;
                        tempbarpattern=song(j-(barsize-1):j);
                        if barcount==1
                            bar_pattern(1,:)=tempbarpattern;
                            bar_patterncount(1)=1;
                            barindex=1;
                            barcount=barcount+1;
                        else
                            if sum(ismember(bar_pattern(:,1:barsize), tempbarpattern,'rows'))==1
                                barindex=find(ismember(bar_pattern(:,1:barsize), tempbarpattern,'rows')==1);
                                bar_patterncount(barindex)=bar_patterncount(barindex)+1;
                                add=0;
                            else
                                bar_pattern(barcount,:)=tempbarpattern;
                                bar_patterncount(barcount)=1;
                                barindex=barcount;
                                barcount=barcount+1;
                                add=1;
                            end
                        end
                        barseq(ii)=barindex;
                        
                        clear barindex
                        clear tempbarpattern
                        MM=length(bar_patterncount);
                        if add==1
                            barcell_pattern{1,1}=[];
                            barcell_pattern{1,1}=[1:length(bar_patterncount)]';
                            barcell_pattern{1,2}(MM,1)=0;
                            barcell_pattern{1,5}(MM,1)=1;
                            if isempty(barcell_pattern{1,3})==1
                                barcell_pattern{1,3}=zeros(MM,MM);
                                barcell_pattern{1,4}=ones(MM,MM)*1/MM;
                            else
                                barcell_pattern{1,3}=[barcell_pattern{1,3},zeros(MM-1,1);zeros(1,MM)];
                                barcell_pattern{1,4}=[barcell_pattern{1,4},ones(MM-1,1)*1/MM;ones(1,MM)*1/MM];
                            end
                            for pp=2:size(barcell_pattern,1)
                                barcell_pattern{pp,3}(:,MM)=0;
                                barcell_pattern{pp,4}(:,MM)=0;
                            end
                            add=0;
                        end
                        %             barcell_pattern{1,2}=[];
                        %             barcell_pattern{1,2}(:,1)=bar_patterncount;
                        %             barcell_pattern{1,5}(:,1)=ones
                        
                        for m=1:M
                            if ii>=m
                                tempbarseq=barseq(ii-(m-1):ii);
                                if barpatcount(m)==0 % for intial pattern of length found
                                    barpatcount(m)=barpatcount(m)+1;
                                    barcell_pattern{m,1}(1,:)=tempbarseq;
                                    barcell_pattern{m,2}(1,1)=1;
                                    barcell_pattern{m,3}(1,:)=zeros(1,MM);
                                    barcell_pattern{m,4}(1,:)=ones(1,MM)*(1/MM);
                                    barcell_pattern{m,5}(1,:)=1;
                                    latestbar=1;
                                    barpatcount(m)=barpatcount(m)+1;
                                    
                                    % now if the pattern already exists in the cell
                                elseif sum(ismember(barcell_pattern{m,1}(:,1:length(tempbarseq)), tempbarseq,'rows'))==1
                                    indexb=find(ismember(barcell_pattern{m,1}(:,1:length(tempbarseq)), tempbarseq,'rows')==1);
                                    barcell_pattern{m,2}(indexb,1)=barcell_pattern{m,2}(indexb,1)+1;
                                    latestbar=indexb;
                                    
                                else % otherwise create a new pattern entry
                                    barcell_pattern{m,1}(barpatcount(m),:)=tempbarseq;
                                    barcell_pattern{m,2}(barpatcount(m),1)=1;
                                    barcell_pattern{m,3}(barpatcount(m),:)=zeros(1,MM);
                                    barcell_pattern{m,4}(barpatcount(m),:)=ones(1,MM)*(1/MM);
                                    barcell_pattern{m,5}(barpatcount(m),:)=1;
                                    
                                    
                                    
                                    latestbar=barpatcount(m);
                                    
                                    barpatcount(m)=barpatcount(m)+1;
                                end
                                
                                
                                
                                if m>1
                                    findindexb=find(ismember(barcell_pattern{m-1,1}(:,1:length(tempbarseq)-1), tempbarseq(1:end-1),'rows')==1);
                                    barcell_pattern{m-1,3}(findindexb,barcell_pattern{m,1}(latestbar,end))=...
                                        barcell_pattern{m-1,3}(findindexb,barcell_pattern{m,1}(latestbar,end))+1;
                                    
                                    
                                    sumcountb=sum(barcell_pattern{m-1,3}(findindexb,1:MM));
                                    if sumcountb~=0
                                        barcell_pattern{m-1,4}(findindexb,1:MM)=barcell_pattern{m-1,3}(findindexb,1:MM)/sumcountb;
                                    end
                                    clear sumcountb
                                    
                                    %                 cell_pattern{n-1}=calcentropy2(cell_pattern{n,4},N,NN);
                                    
                                    barcell_pattern{m-1,5}(findindexb,1)=0;
                                    for iii=1:MM
                                        if barcell_pattern{m-1,4}(findindexb,iii)~=0
                                            barcell_pattern{m-1,5}(findindexb,1)=barcell_pattern{m-1,5}(findindexb,1)...
                                                -(((barcell_pattern{m-1,4}(findindexb,iii)).*log2(barcell_pattern{m-1,4}(findindexb,iii))/log2(MM)));
                                        end
                                    end
                                    if isnan(barcell_pattern{m-1,5}(findindexb))==1
                                        barcell_pattern{m-1,5}(findindexb)=0;
                                    end
                                    probbarindex=find(ismember(barcell_pattern{m-1,1}(:,1:length(tempbarseq)-1), tempbarseq(end-(m-2):end),'rows')==1);
                                    
                                    bararray(m,:)=(2^(m-1))*(1-barcell_pattern{m-1,5}(probbarindex))*barcell_pattern{m-1,4}(probbarindex, :);
                                    
                                elseif m==1
                                    clear bararray
                                    barprob=barcell_pattern{1,2}'/sum(barcell_pattern{1,2});
                                    barentrop=0;
                                    for iii=1:MM
                                        if barprob(iii)~=0
                                            barentrop=barentrop-(((barprob(iii)).*log2(barprob(iii))/log2(MM)));
                                        end
                                    end
                                    if isnan(barentrop)==1
                                        barentrop=0;
                                    end
                                    bararray(m,:)=(2^(m-1))*(1-barentrop)*barprob;
                                    clear barprob
                                    clear barentrop
                                end
                                
                                
                                clear tempbarseq
                                clear indexb
                                clear findindexb
                            end
                        end
                        if ii<floor(length(song)/barsize)
                            %                             if ii==1
                            %                                 bararray=ones(1,MM)*1/MM;
                            %                             end
                            b=sort(bararray,2);
                            scoreb(ii,1:size(b,1))=b(:,end)';
                            if size(bararray,1)==1
                                maxprobb=find(bararray==max(bararray));
                            else
                                newprobb=sum(bararray);
                                maxprobb=find(newprobb==max(newprobb));
                            end
                            
                            if size(maxprobb,2)>1
                                choiceb=randi(size(maxprobb,2));
                                barpredict(ii+1)=maxprobb(choiceb);
                                equalcountb=equalcountb+1;
                                equalindexb(equalcountb)=ii+1;
                            else
                                barpredict(ii+1)=maxprobb;
                            end
                            
                            if barpredict(ii)~=barseq(ii)
                                barerror=barerror+1;
                                errorindexb(barerror)=ii;
                            end
                            clear probsindex
                            clear bararray
                            clear maxindex
                            clear index
                            clear maxprobb
                            clear latestbar
                            clear newprobb
                        end
                    end
                end
                              
                cell_pattern{1,2}(song(j)+1)=cell_pattern{1,2}(song(j)+1)+1;
                for n=2:N
                    if j>=n
                        temppattern=song(j-(n-1):j);
                        if patcount(n)==0 % for intial pattern of length found
                            patcount(n)=patcount(n)+1;
                            cell_pattern{n,1}(1,:)=temppattern;
                            cell_pattern{n,2}(1)=1;
                            cell_pattern{n,3}(1,:)=zeros(1,NN);
                            cell_pattern{n,4}(1,:)=ones(1,NN)*0.5;
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
                            cell_pattern{n,4}(patcount(n),:)=ones(1,NN)*0.5;
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
                            
                            probsarray(n,:)=(2^(n-1))*(1-cell_pattern{n-1,5}(findindex))*cell_pattern{n-1,4}(findindex, :);
                        elseif n==1
                            probsarray(n,:)=(2^(n-1))*(1-cell_pattern{n-1,5}(findindex))*cell_pattern{1,2}'/sum(cell_pattern{1,2});
                            
                        end
                        
                        clear temppattern
                        clear index
                        clear findindex
                        clear latestentry
                        
                    end
                end
                
                if j==1
                    probsarray=ones(1,NN)*0.5;
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
                clear newprobs
                
            end
            
            
            
            barconf=confusionmat(barpredict, barseq);
            barpredscore{barsize/4,N/4}(M-3,1)=1-((barconf(1,2)+barconf(2,1))/sum(sum(barconf)));
            barpredscore{barsize/4,N/4}(M-3,2)=equalcountb;
            
            
            %
            % % predict;
            conf=confusionmat(predict, song);
            predscore{barsize/4,N/4}(M-3,1)=1-((conf(1,2)+conf(2,1))/sum(sum(conf)));
            predscore{barsize/4,N/4}(M-3,2)=equalcount;
            
            
            clear a*
            clear b
            clear bar_pattern*
            clear bararray
            clear barcount
            clear choice*
            clear equalcount*
            clear i*
            clear j*
            clear latest*
            clear m
            clear n
            clear neq
            clear p
            clear pp
            clear patcount
            clear probbarindex
            %% loop clear
            %             clear M*
            %             clear N*
            %             clear barsize
            clear conf
            clear barc*
            clear barer*
            clear barpat*
            clear barpredict
            clear barseq
            clear cell*
            clear e*
            clear predict
            clear score*
        M
        end
        N
    end
    barsize
end