clear all
%% initialisation

% establish users preferences as to how model is initialised and as to
% type of threshold used for the predictability.
random=str2double(inputdlg('0=chosen file as intro, 1=random intro sequence', 'Intro type',1,{'1'}));
if random==0
    introfile=inputdlg('Input file name (DONT INCLUDE FILE EXTENSION)', 'Filename',1,{'barseqintro'});
else
    random=1;
end

prompt={'What type of threshold would you like to use? 1=meanlim, 2=normalised, 3=manual'};
defaultans={'2'};
input=str2double(inputdlg(prompt, 'Type of threshold',1,defaultans));
if input==1
    meanlim=1; normalise=0;
elseif input==2
    meanlim=0; normalise=1;
    predictability=str2double(inputdlg('What level of predictability would you like? (between 1 and 10)','Predictability',1,{'8'}));
elseif input==3
    meanlim=0; normalise=0;
    thresh=str2double(inputdlg('What threshold would you like to use?','Threshold',1,{'16'}));
else
    meanlim=1; normalise=0;
end
% random=1; meanlim=0; normalise=1; predictability=9;
% initialise data
N=32; % max length of pattern analysed
NN=4; % number of different sounds
MAX=100; %length of prediction
l=6; % length of random generated sequence before starting the song creation
error=0;
ptime=0.1; % speed of playback
nmark=0;
jjmark=0;
equalcount=0;

count=0;
change=0;
changeindex=[];
equalindex=[];
maxpred=0;
if normalise==0
    threshold=thresh;
elseif normalise==1
    for pp=1:N
        maxpred=maxpred+NN^(pp-1);
    end
    threshold=(predictability/10)*maxpred;
end

% initialise the song and probability matrices M=10;
if random==1;
    M=l;
    intro=randi(NN,[1 M]);
else
    introname=[introfile{1} '.txt'];
    intro=load(introname);
    M=length(intro);
end
% 
cell_pattern=cell(N,4);
cell_pattern{1,1}=[];
cell_pattern{1,2}=zeros(NN,1);
cell_pattern{1,3}=zeros(NN,NN);
cell_pattern{1,4}=ones(NN,NN)*1/NN;
cell_pattern{1,5}=ones(NN,1);
patcount=zeros(1,N);

piece(1:length(intro))=intro;

for j=1:M+MAX-1;
    probsarray=zeros(N,NN);
    for n=1:N
        if j>=n
            [cell_pattern,patcount,probsarray]=datacollectlearnbar(piece,j,n,patcount,cell_pattern,NN,probsarray);
            
        end
    end
    
    if j>=M
        if size(probsarray,2)>1
        probsarrayscore{j}=max(sum(probsarray));
        scoreent(j-(M-1))=max(sum(probsarray));
        count=count+1;
        else
            probsarrayscore{j}=max(probsarray);
        scoreent(j-(M-1))=max(probsarray);
        count=count+1;
        end
        if size(probsarray,1)==1
            maxprob{j+1}=find(probsarray==max(probsarray));
            minprob=find(probsarray==min(probsarray));
        else
            format long
            newprobs=sum(probsarray);
            maxprob{j+1}=find(newprobs==max(newprobs));
            minprob=find(newprobs==min(newprobs));
        end
        
        % make predictions as to the next note in the piece
        
        if size(maxprob{j+1},2)>1
            choice=randi(size(maxprob{j+1},2));
            piece(j+1)=maxprob{j+1}(choice);
            equalcount=equalcount+1;
            equalindex(equalcount)=j+1;
            
        else
            
            if scoreent(j-(M-1))>threshold %suprise(count)<=initsuprisefact-1
                change=change+1;
                changeindex(change)=j+1;
                piece(j+1)=minprob(randi(size(minprob,2)));
            else
                piece(j+1)=maxprob{j+1};
            end
        end
        clear probsindex
        clear probsarray
        clear maxindex
        clear index
        clear minprob
        clear probmin
        clear probmax
        
    end
    
end


finalpiece=barconversion(piece,NN-1);
[p, newsong]=createsong(finalpiece,32000/8,0,0,0,'binopus10.wav');

finalpieceN3=barconversionN3(piece,NN-1);
[pN3, newsongN3]=createsong(finalpieceN3,32000/8,0,0,0,'binopus10.wav');

finalpieceN4=barconversionN4(piece,NN-1);
[pN4, newsongN4]=createsong(finalpieceN4,32000/8,0,0,0,'binopus10.wav');


% y(1,:)=strrep(num2str(piece),' ','');
% if isempty(changeindex)==0
%     y(2,changeindex)='|';
%     altpeice=piece;
%     for ii=1:length(changeindex)
%         altpeice(changeindex(ii))=maxprob{changeindex(ii)}-1;
%         
%     end
%     if isempty(equalindex)==0
%         y(3,equalindex)=':';
%         y(4,:)=strrep(num2str(altpeice),' ','');
%     else
%         y(3,:)=strrep(num2str(altpeice),' ','');
%     end
% end

% [p,newsong]=createsong(piece,1000,0,0,0);