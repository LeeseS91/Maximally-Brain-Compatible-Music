function [peice, cell_pattern, cell_count, cell_probs, jjmark,nmark,jjval]...
    =datacollectbar(j,n,NN,N, cell_pattern, peice, cell_count,...
    cell_probs)

jjmark=0;
nmark=0;
for jj=1:size(cell_pattern{n,1},1)
    if peice(j-(n-1):j)==cell_pattern{n,1}(jj,:)
        cell_pattern{n,2}(jj)=cell_pattern{n,2}(jj)+1;
        jjmark=jj;
        if rem(jjmark,NN^(n-1))==0
            jjval=NN^(n-1);
            nmark=jjmark/(NN^(n-1));
        else
            jjval=rem(jjmark,NN^(n-1));
            for ncount=1:N
                if jjmark-(ncount*(NN^(n-1)))<=0
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
        break;
    end
end
% update probabilities with new cell counts
for f=1:size(cell_count{n},1)
    sumcount=sum(cell_count{n}(f,1:NN));
    if sumcount~=0
        cell_probs{n}(f,1:NN)=cell_count{n}(f,1:NN)/sumcount;
    end
    clear sumcount
end

end