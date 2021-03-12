function [bsf_fit_var, converH, bsf_solution, metrics, time, Iout] = lshade(SearchAgents_no, dims, Max_iter, lb, ub, Im)

    tic
    format long;
    format compact; 
    rand('seed', sum(100 * clock));

    wnd = floor(1*Max_iter);
    problem_size = dims;
    [ih, ~] = imhist(Im(:,:,1));% histogram  Check Normalizaftion
    max_nfes = SearchAgents_no * Max_iter;
    [sz1,sz2] = size(Im(:,:,1));% im size
    fa=0; % flag


    %% Parameter settings for L-SHADE
    p_best_rate = 0.11; 
    arc_rate = 1.5;
    memory_size = 5;    
    pop_size = 100;
    max_pop_size = pop_size;
    min_pop_size = SearchAgents_no; 

    lu = [lb * ones(1, problem_size); ub * ones(1, problem_size)];

    %fprintf('\n-------------------------------------------------------\n')
    %fprintf(' Threshold = %d\n',  problem_size-1) 
    converH=[]; 

    %% Initialize the main population
    popold = repmat(lu(1, :), pop_size, 1) + rand(pop_size, problem_size) .* (repmat(lu(2, :) - lu(1, :), pop_size, 1));
    pop = popold; % the old population becomes the current population

    %--------------------------------------------------------------------
    [szf,szc] = size(pop);
    for zz = 1:szf  % Evaluation in one shot
      indivi = sort(fix(pop(zz,:)));
      fitness(zz) = CE(indivi,ih);
      fitness = fitness';
    end
    %--------------------------------------------------------------------

    nfes = 0;
    bsf_fit_var = 1e+30;
    bsf_solution = zeros(1, problem_size);

    %%%%%%%%%%%%%%%%%%%%%%%% for out
    [psize1, psize2] = size(pop);
      for i = 1 : psize1
        nfes = nfes + 1;

        if fitness(i) < bsf_fit_var
          bsf_fit_var = fitness(i);
          bsf_solution = pop(i, :);
        end   	 
    	
        if nfes > max_nfes 
          break;
        end

      end
    %%%%%%%%%%%%%%%%%%%%%%%% for out

    memory_sf = 0.5 .* ones(memory_size, 1);
    memory_cr = 0.5 .* ones(memory_size, 1);
    memory_pos = 1;

    archive.NP = arc_rate * pop_size; % Maximum size of the archive
    archive.pop = zeros(0, problem_size); % Solutions stored in te archive
    archive.funvalues = zeros(0, 1); % Function value of the archived solutions

    %% main loop
    while nfes < max_nfes

      pop = popold; % Old population becomes the current population
      [temp_fit, sorted_index] = sort(fitness, 'ascend');

      mem_rand_index = ceil(memory_size * rand(pop_size, 1));
      mu_sf = memory_sf(mem_rand_index);
      mu_cr = memory_cr(mem_rand_index);

      %% Generating crossover rate
      cr = normrnd(mu_cr, 0.1);
      term_pos = find(mu_cr == -1);
      cr(term_pos) = 0;
      cr = min(cr, 1);
      cr = max(cr, 0);

      %% Generating scaling factor
      sf = mu_sf + 0.1 * tan(pi * (rand(pop_size, 1) - 0.5));
      pos = find(sf <= 0);

      while ~ isempty(pos)
        	sf(pos) = mu_sf(pos) + 0.1 * tan(pi * (rand(length(pos), 1) - 0.5));
        	pos = find(sf <= 0);
      end

      sf = min(sf, 1); 
      
      r0 = [1 : pop_size];
      popAll = [pop; archive.pop];
      [r1, r2] = gnR1R2(pop_size, size(popAll, 1), r0);
      
      pNP = max(round(p_best_rate * pop_size), 2); % choose at least two best solutions
      randindex = ceil(rand(1, pop_size) .* pNP); % select from [1, 2, 3, ..., pNP]
      randindex = max(1, randindex); % Avoid the problem that rand = 0 and thus ceil(rand) = 0
      pbest = pop(sorted_index(randindex), :); % Randomly choose one of the top 100p solutions

      vi = pop + sf(:, ones(1, problem_size)) .* (pbest - pop + pop(r1, :) - popAll(r2, :));
      vi = boundConstraint(vi, pop, lu);

      mask = rand(pop_size, problem_size) > cr(:, ones(1, problem_size)); % Mask is used to indicate which elements of ui comes from the parent
      rows = (1 : pop_size)'; cols = floor(rand(pop_size, 1) * problem_size)+1; % Choose one position where the element of ui doesn't come from the parent
      jrand = sub2ind([pop_size problem_size], rows, cols); mask(jrand) = false;
      ui = vi; ui(mask) = pop(mask);

      %-------------------------------------------------------- FEVAL
      [szff,szcc]=size(ui);
      for ez=1:szff  % Evaluation in one shot
        child_target=sort(fix(ui(ez,:)));  
        children_fitness(ez)=CE(child_target,ih);  

        if fa == 0
          children_fitness=children_fitness;
          fa=fa+1;
        elseif fa==1
          children_fitness=children_fitness';
          fa=fa+1;
        else
          children_fitness=children_fitness;
        end

      end
      %-------------------------------------------------------FEVAL

      %%%%%%%%%%%%%%%%%%%%%%%% for out
      for i = 1 : pop_size
        nfes = nfes + 1;

        if children_fitness(i) < bsf_fit_var
          bsf_fit_var = children_fitness(i);
          bsf_solution = ui(i, :);
        end	  
	     
        converH=[converH; bsf_fit_var];% Historic

	      %%	if nfes > max_nfes; exit(1); end
        if nfes > max_nfes
          break;
        end

      end
      %%%%%%%%%%%%%%%%%%%%%%%% for out        
      
      q=0;      
      fitn=fitness(1:szff);
      children_fitness=children_fitness(1:szff);
      dif = abs( fitn - children_fitness );

      %% I == 1: the parent is better; I == 2: the offspring is better
      I = (fitness(1:szff) > children_fitness(1:szff) );
      goodCR = cr(I == 1);  
      goodF = sf(I == 1);
      dif_val = dif(I == 1);
      archive = updateArchive(archive, popold(I == 1, :), fitness(I == 1));
      [fitness, I] = min([fitness, children_fitness], [], 2);
      popold = pop;
      popold(I == 2, :) = ui(I == 2, :);
      num_success_params = numel(goodCR);

      if num_success_params > 0 
        sum_dif = sum(dif_val);
        dif_val = dif_val / sum_dif;

        %% for updating the memory of scaling factor 
        memory_sf(memory_pos) = (dif_val' * (goodF .^ 2)) / (dif_val' * goodF);

        %% for updating the memory of crossover rate
        if max(goodCR) == 0 || memory_cr(memory_pos) == -1
          memory_cr(memory_pos)  = -1;
        else
          memory_cr(memory_pos) = (dif_val' * (goodCR .^ 2)) / (dif_val' * goodCR);
        end

        memory_pos = memory_pos + 1;

        if memory_pos > memory_size  
          memory_pos = 1;
        end
      end

      %% for resizing the population size
      plan_pop_size = round((((min_pop_size - max_pop_size) / max_nfes) * nfes) + max_pop_size);

      if pop_size > plan_pop_size
        reduction_ind_num = pop_size - plan_pop_size;
        
        if pop_size - reduction_ind_num <  min_pop_size 
          reduction_ind_num = pop_size - min_pop_size; 
        end

        pop_size = pop_size - reduction_ind_num;
          
        for r = 1 : reduction_ind_num
          [valBest indBest] = sort(fitness, 'ascend');
          worst_ind = indBest(end);
          popold(worst_ind,:) = [];
          pop(worst_ind,:) = [];
          fitness(worst_ind,:) = [];
        end
              	  
        archive.NP = round(arc_rate * pop_size); 

        if size(archive.pop, 1) > archive.NP 
          rndpos = randperm(size(archive.pop, 1));
          rndpos = rndpos(1 : archive.NP);
          archive.pop = archive.pop(rndpos, :);
        end

      end % end if 

    end  % end while
    size(bsf_solution);
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
    converH = converH';
    bsf_solution = sort(BThresholds);

  
%% Function to generate column vector
function [r1, r2] = gnR1R2(NP1, NP2, r0)

    % gnR1R2 generate two column vectors r1 and r2 of size NP1 & NP2, respectively
    %    r1's elements are choosen from {1, 2, ..., NP1} & r1(i) ~= r0(i)
    %    r2's elements are choosen from {1, 2, ..., NP2} & r2(i) ~= r1(i) & r2(i) ~= r0(i)

    NP0 = length(r0);
    r1 = floor(rand(1, NP0) * NP1) + 1;
    %for i = 1 : inf
    for i = 1 : 99999999
      pos = (r1 == r0);
      if sum(pos) == 0
        break;
      else % regenerate r1 if it is equal to r0
        r1(pos) = floor(rand(1, sum(pos)) * NP1) + 1;
      end
      if i > 1000 % this has never happened so far
        error('Can not genrate r1 in 1000 iterations');
      end
    end

    r2 = floor(rand(1, NP0) * NP2) + 1;
    %for i = 1 : inf
    for i = 1 : 99999999
      pos = ((r2 == r1) | (r2 == r0));
      if sum(pos)==0
        break;
      else % regenerate r2 if it is equal to r0 or r1
        r2(pos) = floor(rand(1, sum(pos)) * NP2) + 1;
      end
      if i > 1000 % this has never happened so far
        error('Can not genrate r2 in 1000 iterations');
      end
    end

    
%% Function for bound constraints
function vi = boundConstraint (vi, pop, lu)

    % if the boundary constraint is violated, set the value to be the middle
    % of the previous value and the bound

    [NP, D] = size(pop);  % the population size and the problem's dimension

    %% check the lower bound
    xl = repmat(lu(1, :), NP, 1);
    pos = vi < xl;
    vi(pos) = (pop(pos) + xl(pos)) / 2;

    %% check the upper bound
    xu = repmat(lu(2, :), NP, 1);
    pos = vi > xu;
    vi(pos) = (pop(pos) + xu(pos)) / 2;

    
%% Function to update archive
function archive = updateArchive(archive, pop, funvalue)
    % Update the archive with input solutions
    % Step 1: Add new solution to the archive
    % Step 2: Remove duplicate elements
    % Step 3: If necessary, randomly remove some solutions to maintain the archive size

    if archive.NP == 0, return; end

    if size(pop, 1) ~= size(funvalue,1), error('check it'); end

    %% Method 2: Remove duplicate elements
    popAll = [archive.pop; pop ];
    funvalues = [archive.funvalues; funvalue ];
    [dummy IX]= unique(popAll, 'rows');
    if length(IX) < size(popAll, 1) % There exist some duplicate solutions
      popAll = popAll(IX, :);
      funvalues = funvalues(IX, :);
    end

    if size(popAll, 1) <= archive.NP   % add all new individuals
      archive.pop = popAll;
      archive.funvalues = funvalues;
    else                % randomly remove some solutions
      rndpos = randperm(size(popAll, 1)); % equivelent to "randperm";
      rndpos = rndpos(1 : archive.NP);
      archive.pop = popAll  (rndpos, :);
      archive.funvalues = funvalues(rndpos, :);
    end
   