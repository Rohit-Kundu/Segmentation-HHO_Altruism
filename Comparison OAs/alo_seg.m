% SearchAgents_no, dims, Max_iter, lb, ub, Im
function [Elite_antlion_fitness,Convergence_curve,Elite_antlion_position, metrics, time, Iout]=alo_seg(N,dim,Max_iter,lb,ub,Im)

tic
format long;
format compact; 
rand('seed', sum(100 * clock));

% wnd = floor(1*Max_iter);
% problem_size = dims;
[ih, ~] = imhist(Im(:,:,1));% histogram  Check Normalizaftion
% max_nfes = SearchAgents_no * Max_iter;
[sz1,sz2] = size(Im(:,:,1));% im size
fa=0; % flag

% Initialize the positions of antlions and ants
antlion_position=initializationALO(N,dim,ub,lb);
ant_position=initializationALO(N,dim,ub,lb);

% Initialize variables to save the position of elite, sorted antlions, 
% convergence curve, antlions fitness, and ants fitness
Sorted_antlions=zeros(N,dim);
Elite_antlion_position=zeros(1,dim);
Elite_antlion_fitness=inf;
Convergence_curve=zeros(1,Max_iter);
antlions_fitness=zeros(1,N);
ants_fitness=zeros(1,N);

% Calculate the fitness of initial antlions and sort them
for i=1:size(antlion_position,1)
%     antlions_fitness(1,i)=fobj(antlion_position(i,:)); 
    indivi = sort(fix(antlion_position(i,:)));
    antlions_fitness(1,i) = CE(indivi,ih);
%     fitness = fitness';
end

[sorted_antlion_fitness,sorted_indexes]=sort(antlions_fitness);
    
for newindex=1:N
     Sorted_antlions(newindex,:)=antlion_position(sorted_indexes(newindex),:);
end
    
Elite_antlion_position=Sorted_antlions(1,:);
Elite_antlion_fitness=sorted_antlion_fitness(1);

% Main loop start from the second iteration since the first iteration 
% was dedicated to calculating the fitness of antlions
Current_iter=2; 
while Current_iter<Max_iter+1
    
    % This for loop simulate random walks
    for i=1:size(ant_position,1)
        % Select ant lions based on their fitness (the better anlion the higher chance of catching ant)
        Rolette_index=RouletteWheelSelection(1./sorted_antlion_fitness);
        if Rolette_index==-1  
            Rolette_index=1;
        end
      
        % RA is the random walk around the selected antlion by rolette wheel
        RA=Random_walk_around_antlion(dim,Max_iter,lb,ub, Sorted_antlions(Rolette_index,:),Current_iter);
        
        % RA is the random walk around the elite (best antlion so far)
        [RE]=Random_walk_around_antlion(dim,Max_iter,lb,ub, Elite_antlion_position(1,:),Current_iter);
        
        ant_position(i,:)= (RA(Current_iter,:)+RE(Current_iter,:))/2; % Equation (2.13) in the paper          
    end
    
    for i=1:size(ant_position,1)  
        
        % Boundar checking (bring back the antlions of ants inside search
        % space if they go beyoud the boundaries
        Flag4ub=ant_position(i,:)>ub;
        Flag4lb=ant_position(i,:)<lb;
        ant_position(i,:)=(ant_position(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;  
        
        indivi = sort(fix(ant_position(i,:)));
        ants_fitness(1,i) = CE(indivi,ih);
%         ants_fitness(1,i)=fobj(ant_position(i,:));        
       
    end
    
    % Update antlion positions and fitnesses based of the ants (if an ant 
    % becomes fitter than an antlion we assume it was cought by the antlion  
    % and the antlion update goes to its position to build the trap)
    double_population=[Sorted_antlions;ant_position];
    double_fitness=[sorted_antlion_fitness ants_fitness];
        
    [double_fitness_sorted I]=sort(double_fitness);
    double_sorted_population=double_population(I,:);
        
    antlions_fitness=double_fitness_sorted(1:N);
    Sorted_antlions=double_sorted_population(1:N,:);
        
    % Update the position of elite if any antlinons becomes fitter than it
    if antlions_fitness(1)<Elite_antlion_fitness 
        Elite_antlion_position=Sorted_antlions(1,:);
        Elite_antlion_fitness=antlions_fitness(1);
    end
      
    % Keep the elite in the population
    Sorted_antlions(1,:)=Elite_antlion_position;
    antlions_fitness(1)=Elite_antlion_fitness;
  
    % Update the convergence curve
    Convergence_curve(Current_iter)=Elite_antlion_fitness;

    % Display the iteration and best optimum obtained so far
%     if mod(Current_iter,50)==0
%         display(['At iteration ', num2str(Current_iter), ' the elite fitness is ', num2str(Elite_antlion_fitness)]);
%     end

    Current_iter=Current_iter+1; 
end

bsf_solution=Elite_antlion_position;
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
Elite_antlion_position=bsf_solution;




function [RWs]=Random_walk_around_antlion(Dim,max_iter,lb, ub,antlion,current_iter)
if size(lb,1) ==1 && size(lb,2)==1 %Check if the bounds are scalar
    lb=ones(1,Dim)*lb;
    ub=ones(1,Dim)*ub;
end

if size(lb,1) > size(lb,2) %Check if boundary vectors are horizontal or vertical
    lb=lb';
    ub=ub';
end

I=1; % I is the ratio in Equations (2.10) and (2.11)

if current_iter>max_iter/10
    I=1+100*(current_iter/max_iter);
end

if current_iter>max_iter/2
    I=1+1000*(current_iter/max_iter);
end

if current_iter>max_iter*(3/4)
    I=1+10000*(current_iter/max_iter);
end

if current_iter>max_iter*(0.9)
    I=1+100000*(current_iter/max_iter);
end

if current_iter>max_iter*(0.95)
    I=1+1000000*(current_iter/max_iter);
end


% Dicrease boundaries to converge towards antlion
lb=lb/(I); % Equation (2.10) in the paper
ub=ub/(I); % Equation (2.11) in the paper

% Move the interval of [lb ub] around the antlion [lb+anlion ub+antlion]
if rand<0.5
    lb=lb+antlion; % Equation (2.8) in the paper
else
    lb=-lb+antlion;
end

if rand>=0.5
    ub=ub+antlion; % Equation (2.9) in the paper
else
    ub=-ub+antlion;
end

% This function creates n random walks and normalize accroding to lb and ub
% vectors 
for i=1:Dim
    X = [0 cumsum(2*(rand(max_iter,1)>0.5)-1)']; % Equation (2.1) in the paper
    %[a b]--->[c d]
    a=min(X);
    b=max(X);
    c=lb(i);
    d=ub(i);      
    X_norm=((X-a).*(d-c))./(b-a)+c; % Equation (2.7) in the paper
    RWs(:,i)=X_norm;
end




function choice = RouletteWheelSelection(weights)
  accumulation = cumsum(weights);
  p = rand() * accumulation(end);
  chosen_index = -1;
  for index = 1 : length(accumulation)
    if (accumulation(index) > p)
      chosen_index = index;
      break;
    end
  end
  choice = chosen_index;
  
  
  
  
function X=initializationALO(SearchAgents_no,dim,ub,lb)

Boundary_no= size(ub,2); % numnber of boundaries

% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1
    X=rand(SearchAgents_no,dim).*(ub-lb)+lb;
end

% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        X(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
    end
end