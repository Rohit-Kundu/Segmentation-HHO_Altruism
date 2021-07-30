function [Fbest, alfa, Xbest,metrics, time, Iout] = gsa_seg(SearchAgents_no,dim,Max_iter,Lb,Ub,Im)
format long;
format compact; 
rand('seed', sum(100 * clock));

[ih, ~] = imhist(Im(:,:,1));% histogram  Check Normalizaftion
max_nfes = SearchAgents_no * Max_iter;
[sz1,sz2] = size(Im(:,:,1));% im size
fa=0; % flag

format short e
NGen=10000;        % Maximum Number of Iterations
NP = 100;          % Population Size
G0 = 100;         % GSA G0
alfa = 20;        % GSA alfa
Rpower = 1;       % GSA power of 'R'
ElitistCheck=1; 
X = Lb+(Ub-Lb).*rand(NP,dim);
V = zeros(NP, dim);
for i = 1:NP
    FU=X(i,:)>Ub;FL=X(i,:)<Lb;X(i,:)=(X(i,:).*(~(FU+FL)))+Ub.*FU+Lb.*FL;
        % fitness of locations
    indivi = sort(fix(X(i,:)));
    fitness(i,1) = CE(indivi,ih);
    fitness = fitness';
%     fitness(i,1)=fitnessobj(func_num,X(i,:));
        
%     fitness(j,1) = feval(fhd,X(j,:),funcNum);
end
[best , bestX]=min(fitness);                  % minimization
Fbest=best;Xbest=X(bestX,:);
BestChart = [];
BestChart=[BestChart; Fbest];
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for it=2:NGen
    [M] = massCalculation(fitness);  % Calculation of M
    G = G0*exp(-alfa*it/NGen);       % Calculation of Gravitational constant
    a = Gfield(M,X,G,Rpower,ElitistCheck,it,NGen); %Calculation of accelaration
    V=rand(NP,dim).*V+a;
    X=X+V;
    for j = 1:NP
        Tp=X(j,:)>Ub;
        Tm=X(j,:)<Lb;
        X(j,:)=(X(j,:).*(~(Tp+Tm)))+((rand(1,dim).*(Ub-Lb)+Lb).*(Tp+Tm));  
        indivi = sort(fix(X(j,:)));
        fitness(j,1) = CE(indivi,ih);   
        fitness = fitness';
%         fitness(j,1)=fitnessobj(func_num,X(j,:));
        
%         fitness(j,1) = feval(fhd,X(j,1:dim),funcNum);     
    end
    [best , bestX]=min(fitness);                  % minimization
    if best<Fbest  %minimization.
        Fbest=best;Xbest=X(bestX,:);
    end     
    BestChart=[BestChart;Fbest];
end    %iteration
bsf_solution=Xbest;
bsf_solution = fix(bsf_solution);
BThresholds = sort(bsf_solution);

gBestR = BThresholds;   
Iout = imageGRAY(Im,gBestR); %Segmented Image

%% Metric calculation

window = 3;
Imd = double(Im);
Ioutd = double(Iout);
psnr = PSNR(Im, Iout); %1
SSIM = ssim(Im, Iout); %2
FSIM = FeatureSIM(Im, Iout); %3
UIQI = img_qi(Im, Iout, window); %4
QILV = qilv(Im, Iout, [window, window]); %5
HPSI = HaarPSI(Imd, Ioutd); %6
metrics = [psnr, SSIM, FSIM, UIQI, QILV, HPSI];

time = toc;
%converH = converH';
bsf_solution = sort(BThresholds);
Xbest=bsf_solution;


function [M] = massCalculation(fitness)
Fmax=max(fitness); Fmin=min(fitness); 
[N , ~]=size(fitness);
   if Fmax==Fmin
       M=ones(N,1);     
   else
      best=Fmin;worst=Fmax; 
      M=(fitness - worst)./(best - worst);
   end
M=M./sum(M);


function a = Gfield(M,X,G,Rpower,ElitistCheck,it,NGen)
 [N,dim]=size(X);
 final_per=2; %In the last iteration, only 2 percent of agents apply force to the others.
 if ElitistCheck==1
     kbest=final_per+(1-it/NGen)*(100-final_per); %kbest in eq. 21.
     kbest=round(N*kbest/100);
 else
     kbest=N; 
 end
 [~ , ds]=sort(M,'descend');
 for i=1:N
     E(i,:)=zeros(1,dim);
     for ii=1:kbest
         j=ds(ii);
         if j~=i
            R=norm(X(i,:)-X(j,:),2); %Euclidian distanse.
            for k=1:dim 
                E(i,k)=E(i,k)+rand*(M(j))*((X(j,k)-X(i,k))/(R^Rpower+eps)); %note that Mp(i)/Mi(i)=1              
            end
         end
     end
 end
a=E.*G; % acceleration ,note that Mp(i)/Mi(i)=1