function cell_ent=calcentropylearn(cell_probs,N,NN)

entrop=cell(N,1);


for nn=1:N
    for ii=1:size(cell_probs{nn},1)
        sumprob=sum(cell_probs{nn}(ii,:));
        if sumprob==0
            entrop{nn}(ii)=1;
        elseif sum(isnan(cell_probs{nn}(ii,:)))>=1
            entrop{nn}(ii)=1;
        else
            for iii=1:NN
                if isnan((cell_probs{nn,1}(ii,iii).*log2(cell_probs{nn,1}(ii,iii)))/log2(NN))==0
            entrop{nn}(ii)=entrop{nn}(ii)-((cell_probs{nn,1}(ii,iii).*log2(cell_probs{nn,1}(ii,iii)))/log2(NN));
                end
            end
            
            if isnan(entrop{nn}(ii))==1
                entrop{nn}(ii)=0;
            end
        end
    end
end
cell_ent=entrop;