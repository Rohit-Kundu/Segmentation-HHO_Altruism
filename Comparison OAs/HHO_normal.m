function [Rabbit_Energy,Rabbit_Location,metrics,time, Iout]=HHO_normal(N,T,lb,ub,dim,Im)
tic
format long;
format compact; 
rand('seed', sum(100 * clock));
%wnd = floor(1*Max_iter);
problem_size = dim;
[ih, ~] = imhist(Im(:,:,1));% histogram  Check Normalization
%max_nfes = SearchAgents_no * Max_iter;
[sz1,sz2] = size(Im(:,:,1));% im size
fa=0;

% initialize the location and Energy of the rabbit
Rabbit_Location=zeros(1,dim);
Rabbit_Energy=inf;

%Initialize the locations of Harris' hawks
X=initialization(N,dim,ub,lb);
size(Rabbit_Location);
[szf,szc] = size(X);

%szf is the size of the X vector 

for zz = 1:szf  % Evaluation in one shot
  indivi = sort(fix(X(zz,:)));
  fitness(zz) = CE(indivi,ih);
  fitness = fitness';
end
converH=[];
%CNVG=zeros(1,T);

t=0; % Loop counter

while t<T
    for i=1:size(X,1)
        % Check boundries
        FU=X(i,:)>ub;FL=X(i,:)<lb;X(i,:)=(X(i,:).*(~(FU+FL)))+ub.*FU+lb.*FL;
        % fitness of locations
        indivi = sort(fix(X(i,:)));
        fitness(i) = CE(indivi,ih);
        %indivi
        %ih
        %fitness=f_obj(X(i,:));
        % Update the location of Rabbit
        if fitness(i)<Rabbit_Energy
            Rabbit_Energy=fitness(i);
            Rabbit_Location=X(i,:);
        end
    end
    
    E1=2*(1-(t/T)); % factor to show the decreaing energy of rabbit
    % Update the location of Harris' hawks
    for i=1:size(X,1)
        E0=2*rand()-1; %-1<E0<1
        Escaping_Energy=E1*(E0);  % escaping energy of rabbit
        
        if abs(Escaping_Energy)>=1
            %% Exploration:
            % Harris' hawks perch randomly based on 2 strategy:
            
            q=rand();
            rand_Hawk_index = floor(N*rand()+1);
            X_rand = X(rand_Hawk_index, :);
            if q<0.5
                % perch based on other family members
                X(i,:)=X_rand-rand()*abs(X_rand-2*rand()*X(i,:));
            elseif q>=0.5
                % perch on a random tall tree (random site inside group's home range)
                X(i,:)=(Rabbit_Location(1,:)-mean(X))-rand()*((ub-lb)*rand+lb);
            end
            
        elseif abs(Escaping_Energy)<1
            %% Exploitation:
            % Attacking the rabbit using 4 strategies regarding the behavior of the rabbit
            
            %% phase 1: surprise pounce (seven kills)
            % surprise pounce (seven kills): multiple, short rapid dives by different hawks
            
            r=rand(); % probablity of each event
            temp_var=0;
            if r>=0.5 && abs(Escaping_Energy)<0.5 % Hard besiege
                X(i,:)=(Rabbit_Location)-Escaping_Energy*abs(Rabbit_Location-X(i,:));
            end
            
            if r>=0.5 && abs(Escaping_Energy)>=0.5  % Soft besiege
                Jump_strength=2*(1-rand()); % random jump strength of the rabbit
                X(i,:)=(Rabbit_Location-X(i,:))-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:));
            end
            
            %% phase 2: performing team rapid dives (leapfrog movements)
            if r<0.5 && abs(Escaping_Energy)>=0.5, % Soft besiege % rabbit try to escape by many zigzag deceptive motions
                
                Jump_strength=2*(1-rand());
                X1=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:));
                
                %a___=CE(X1',ih);
        %CE(sort(fix(X(i,:))),ih);
                X1=check_lb_n_ub(X1, lb, ub);
                temp_var=E0;
                    
                
                indivi2=CE(sort(fix(X(i,:))),ih);
                
                indivi = CE(sort(fix(X1)),ih);
                if indivi<indivi2 % improved move?
                    %[X1, X(i, :)]=altr(X1, X(i, :),temp_var , lb, ub, ih);
                
                    X(i,:)=X1;
                else % hawks perform levy-based short rapid dives around the rabbit
                    X2=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:))+rand(1,dim).*Levy(dim);
                    
                    X2=check_lb_n_ub(X2, lb, ub);
                    indivi2=sort(fix(X(i,:)));
                    indivi=sort(fix(X2));
                    if (CE(indivi,ih)<CE(indivi2,ih)), % improved move?
                        X(i,:)=X2;
                        %temp_var=E0;
                        %[X2, X(i, :)]=altr(X2, X(i, :),temp_var , lb, ub, ih);                    
                    end
                end
            end
            
            if r<0.5 && abs(Escaping_Energy)<0.5, % Hard besiege % rabbit try to escape by many zigzag deceptive motions
                % hawks try to decrease their average location with the rabbit
                Jump_strength=2*(1-rand());
                X1=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-mean(X));
                temp_var=E0*(1/szf);
                X1=check_lb_n_ub(X1, lb, ub);
                indivi2=sort(fix(X(i,:)));
                indivi=sort(fix(X1));
                if CE(indivi,ih)<CE(indivi2,ih) % improved move?
                    %[X1, X(i, :)]=altr(X1, X(i, :),temp_var , lb, ub, ih);
                    X(i,:)=X1;
                else % Perform levy-based short rapid dives around the rabbit
                    X2=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-mean(X))+rand(1,dim).*Levy(dim);
                    
                    X2=check_lb_n_ub(X2, lb, ub);
                    temp_var=E0*(1/szf);
                    indivi=sort(fix(X2));indivi2=sort(fix(X(i,:)));
                    if (CE(indivi,ih)<CE(indivi2,ih)), % improved move?
                        %[X2, X(i, :)]=altr(X2, X(i, :),temp_var , lb, ub, ih);
                    
                        X(i,:)=X2;
                    end
                end
            end
            %%
        end
    end
    t=t+1;
    %CNVG(t)=Rabbit_Energy;
%    Print the progress every 100 iterations
%    if mod(t,100)==0
%        display(['At iteration ', num2str(t), ' the best fitness is ', num2str(Rabbit_Energy)]);
%    end
end
bsf_solution=Rabbit_Location;
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
Rabbit_Location=bsf_solution;

end

% ___________________________________
function o=Levy(d)
beta=1.5;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;v=randn(1,d);step=u./abs(v).^(1/beta);
o=step;
end


function [X1]=check_lb_n_ub(X1, lb, ub)
    for i = 1:size(X1,2)
        if X1(i)<lb
            X1(i)=lb;
        end
        
        if X1(i)>ub
            X1(i)=ub;
        end
    end
end

function [X_exploit, X_original]=altr(X_exploit, X_original, r, lb, ub, ih)
    xtra=1+ (floor(rand()*4));% an xtra no that can be added or subtracted to the X_exploit between 0 and 3
    for i = 1:size(X_exploit)
        b=CE(sort(fix(X_exploit)),ih);
        c=CE(sort(fix(X_original)),ih);
        X_exploit(i)=X_exploit(i)+xtra;
        X_original(i)=X_original(i)-xtra;    
        X_exploit=check_lb_n_ub(X_exploit, lb, ub);
        X_original=check_lb_n_ub(X_original, lb, ub);
        b=b-CE(sort(fix(X_exploit)),ih);
        c=CE(sort(fix(X_original)),ih)-c;
        if r*b<=c || b<0
            X_exploit(i)=X_exploit(i)-xtra;
            X_original(i)=X_original(i)+xtra;    
            X_exploit=check_lb_n_ub(X_exploit, lb, ub);
            X_original=check_lb_n_ub(X_original, lb, ub);
        end
        
        b=CE(sort(fix(X_exploit)),ih);
        c=CE(sort(fix(X_original)),ih);
        X_exploit(i)=X_exploit(i)-xtra;
        X_original(i)=X_original(i)+xtra;    
        X_exploit=check_lb_n_ub(X_exploit, lb, ub);
        X_original=check_lb_n_ub(X_original, lb, ub);
        b=b-CE(sort(fix(X_exploit)),ih);
        c=CE(sort(fix(X_original)),ih)-c;
        if r*b<=c || b<0
            X_exploit(i)=X_exploit(i)+xtra;
            X_original(i)=X_original(i)-xtra;    
            X_exploit=check_lb_n_ub(X_exploit, lb, ub);
            X_original=check_lb_n_ub(X_original, lb, ub);
        end
        
    end
    %xtra
end

                    
% % no_search_agents, no_itns, lb, ub, dim, fobj
% %function [Rabbit_Energy,Rabbit_Location,CNVG]=HHO(N,T,lb,ub,dim,fobj)
% function [Rabbit_Energy,Rabbit_Location,metrics,time, Iout]=HHO(N,T,lb,ub,dim,Im)
% 
% %disp('HHO is now tackling your problem')
% tic
% format long;
% format compact; 
% rand('seed', sum(100 * clock));
% 
% %wnd = floor(1*Max_iter);
% problem_size = dim;
% [ih, ~] = imhist(Im(:,:,1));% histogram  Check Normalization
% %max_nfes = SearchAgents_no * Max_iter;
% [sz1,sz2] = size(Im(:,:,1));% im size
% fa=0;
% 
% % initialize the location and Energy of the rabbit
% Rabbit_Location=zeros(1,dim);
% Rabbit_Energy=inf;
% 
% %Initialize the locations of Harris' hawks
% X=initialization(N,dim,ub,lb);
% size(Rabbit_Location);
% [szf,szc] = size(X);
% 
% for zz = 1:szf  % Evaluation in one shot
%   indivi = sort(fix(X(zz,:)));
%   fitness(zz) = CE(indivi,ih);
%   fitness = fitness';
% end
% converH=[];
% %CNVG=zeros(1,T);
% 
% t=0; % Loop counter
% 
% while t<T
%     for i=1:size(X,1)
%         % Check boundries
%         FU=X(i,:)>ub;FL=X(i,:)<lb;X(i,:)=(X(i,:).*(~(FU+FL)))+ub.*FU+lb.*FL;
%         % fitness of locations
%         indivi = sort(fix(X(i,:)));
%         fitness(i) = CE(indivi,ih);
%         %fitness=f_obj(X(i,:));
%         % Update the location of Rabbit
%         if fitness(i)<Rabbit_Energy
%             Rabbit_Energy=fitness(i);
%             Rabbit_Location=X(i,:);
%         end
%     end
%     
%     E1=2*(1-(t/T)); % factor to show the decreaing energy of rabbit
%     % Update the location of Harris' hawks
%     for i=1:size(X,1)
%         E0=2*rand()-1; %-1<E0<1
%         Escaping_Energy=E1*(E0);  % escaping energy of rabbit
%         
%         if abs(Escaping_Energy)>=1
%             %% Exploration:
%             % Harris' hawks perch randomly based on 2 strategy:
%             
%             q=rand();
%             rand_Hawk_index = floor(N*rand()+1);
%             X_rand = X(rand_Hawk_index, :);
%             if q<0.5
%                 % perch based on other family members
%                 X(i,:)=X_rand-rand()*abs(X_rand-2*rand()*X(i,:));
%             elseif q>=0.5
%                 % perch on a random tall tree (random site inside group's home range)
%                 X(i,:)=(Rabbit_Location(1,:)-mean(X))-rand()*((ub-lb)*rand+lb);
%             end
%             
%         elseif abs(Escaping_Energy)<1
%             %% Exploitation:
%             % Attacking the rabbit using 4 strategies regarding the behavior of the rabbit
%             
%             %% phase 1: surprise pounce (seven kills)
%             % surprise pounce (seven kills): multiple, short rapid dives by different hawks
%             
%             r=rand(); % probablity of each event
%             
%             if r>=0.5 && abs(Escaping_Energy)<0.5 % Hard besiege
%                 X(i,:)=(Rabbit_Location)-Escaping_Energy*abs(Rabbit_Location-X(i,:));
%             end
%             
%             if r>=0.5 && abs(Escaping_Energy)>=0.5  % Soft besiege
%                 Jump_strength=2*(1-rand()); % random jump strength of the rabbit
%                 X(i,:)=(Rabbit_Location-X(i,:))-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:));
%             end
%             
%             %% phase 2: performing team rapid dives (leapfrog movements)
%             if r<0.5 && abs(Escaping_Energy)>=0.5, % Soft besiege % rabbit try to escape by many zigzag deceptive motions
%                 
%                 Jump_strength=2*(1-rand());
%                 X1=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:));
%                 
%                 %a___=CE(X1',ih);
%         %CE(sort(fix(X(i,:))),ih);
%                 indivi2=CE(sort(fix(X(i,:))),ih);
%                 
%                 %indivi = CE(sort(fix(X1)),ih);
%                 if indivi2<indivi2 % improved move?
%                     X(i,:)=X1;
%                 else % hawks perform levy-based short rapid dives around the rabbit
%                     X2=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:))+rand(1,dim).*Levy(dim);
%                     indivi2=sort(fix(X(i,:)));
%                     %indivi=sort(fix(X2));
%                     if (CE(indivi2,ih)<CE(indivi2,ih)), % improved move?
%                         X(i,:)=X2;
%                     end
%                 end
%             end
%             
%             if r<0.5 && abs(Escaping_Energy)<0.5, % Hard besiege % rabbit try to escape by many zigzag deceptive motions
%                 % hawks try to decrease their average location with the rabbit
%                 Jump_strength=2*(1-rand());
%                 X1=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-mean(X));
%                 
%                 indivi2=sort(fix(X(i,:)));
%                 %indivi=sort(fix(X1));
%                 if CE(indivi2,ih)<CE(indivi2,ih) % improved move?
%                     X(i,:)=X1;
%                 else % Perform levy-based short rapid dives around the rabbit
%                     X2=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-mean(X))+rand(1,dim).*Levy(dim);
%                     %indivi=sort(fix(X2));
%                     indivi2=sort(fix(X(i,:)));
%                     if (CE(indivi2,ih)<CE(indivi2,ih)), % improved move?
%                         X(i,:)=X2;
%                     end
%                 end
%             end
%             %%
%         end
%     end
%     t=t+1;
%     %CNVG(t)=Rabbit_Energy;
% %    Print the progress every 100 iterations
% %    if mod(t,100)==0
% %        display(['At iteration ', num2str(t), ' the best fitness is ', num2str(Rabbit_Energy)]);
% %    end
% end
% bsf_solution=Rabbit_Location;
% bsf_solution = fix(bsf_solution);
% BThresholds = sort(bsf_solution);
% 
% gBestR = BThresholds;   
% Iout = imageGRAY(Im,gBestR); %Segmented Image
% 
% %% Metric calculation
% 
% window = 3;
% Imd = double(Im);
% Ioutd = double(Iout);
% psnr = PSNR(Im, Iout); %1
% SSIM = ssim(Im, Iout); %2
% FSIM = FeatureSIM(Im, Iout); %3
% UIQI = img_qi(Im, Iout, window); %4
% QILV = qilv(Im, Iout, [window, window]); %5
% HPSI = HaarPSI(Imd, Ioutd); %6
% metrics = [psnr, SSIM, FSIM, UIQI, QILV, HPSI];
% 
% time = toc;
% %converH = converH';
% bsf_solution = sort(BThresholds);
% Rabbit_Location=bsf_solution;
% 
% end
% 
% % ___________________________________
% function o=Levy(d)
% beta=1.5;
% sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
% u=randn(1,d)*sigma;v=randn(1,d);step=u./abs(v).^(1/beta);
% o=step;
% end
% 
