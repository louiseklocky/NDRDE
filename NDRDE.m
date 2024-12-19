function output=NDRDE(FES,D,lu,f,benchmark,run_id,pop_size_NDR)

ub=lu(2,:);
lb=lu(1,:);
f=input.f;
bsf_fit_var=e+10^30;
optimalChart=[];
optimalChart_interval=1000;

p_best_rate = 0.11;
arc_rate = 1.2;
memory_size = 5;
pop_size = round(0.5*D);
max_pop_size = pop_size;
min_pop_size = 10;
memory_sf = 0.5 .* ones(memory_size, 1);
memory_cr = 0.5 .* ones(memory_size, 1);
memory_pos = 1;
archive.NP = arc_rate * pop_size;
archive.popu = zeros(0, D);
archive.funvalues = zeros(0, 1);

[popu, ~, ~, ~, ~, ~, ~, ~, ~] = Initialization(pop_size_NDR, D, benchmark, f);
fit=Evaluation(popu, benchmark, f);
fit = fit';
nFES=pop_size_NDR;

ibbb=1;
deit=1;
ibb=0;
bb=0;
popuOld = popu;
fitOld = fit;
tier=1;
TB=30000;
TD=70000;
bt=bsf_fit_var;

while nFES < FES

b = 1-nFES/FES;

if tier==1
if nFES>=TD
      if min(fit)<0.9*bt
         TD=TD+TB;
      end
      bt=bsf_fit_var;
end
end

%NDR%

if nFES<TD
    [~,I] = sort(fit);
    SF = zeros(1,popuSize);
    DSF = zeros(1,popuSize);
    for i = 1:popuSize
        SF(I(i)) = i / popuSize+ 0.1*randn;
        if SF(I(i)) < 0
            SF(I(i)) = i / popuSize;
        end
            if nFES<=2*popuSize
               DSF(i) = 0 + randperm(2,1);
            else
                if II(i) == 2
                   DSF(i)=dsf(i);
                else
                   if rand*dsf(i)<=0.2
                      DSF(i)=dsf(i) + randperm(2,1);
                   else
                      DSF(i)=dsf(i) - randperm(2,1);
                   end
                   if DSF(i)<=1
                      DSF(i)=0 + randperm(2,1);
                   end
                end
            end
    end
    dsf=DSF;
            
        for i = 1:popuSize
        randPopuList = randperm(popuSize);
        randPopuList = setdiff(randPopuList,1,'stable');
        indiR1 = popu(randPopuList(1),:);
        indiR2 = popu(randPopuList(2),:);
        indiR5 = popu(randPopuList(5),:);
        randDimList = randperm(D);
        randDSFList = randperm(DSF(i));
        seleDim = sort(randDimList(1:randDSFList(1)));
        seleDimLen = length(seleDim);
        if seleDimLen == 1
            s1= SF(i)*HyperSphereTransform_1D(indiR1,indiR2,seleDim);
        elseif seleDimLen == 2
            s1= SF(i)*HyperSphereTransform_2D(indiR1,indiR2,seleDim);
        else
            s1= SF(i)*HyperSphereTransform(indiR1,indiR2,seleDim);
        end
        if s1==0
           s1=0.05*(rand*(ub(1)-lb(1))+lb(1));
        end
           popu(i,seleDim) = indiR5(seleDim) + s1;
        end
            
            [popu]= BoundaryDetection(popu,lu);
            [fit]=Evaluation(popu, benchmark, f);
            fit=fit';

            [fitOld, II] = min([fitOld,fit], [], 2);
            popuOld(II == 2, :) = popu(II == 2, :);
            popu = popuOld;
            popold=popu;
            fit=fitOld;

            % get bsf_fit_var
            for i=1:popuSize
                nFES = nFES + 1;
                if min(fit) < bsf_fit_var
                   bsf_fit_var = min(fit);
                end
                if mod(nFES, optimalChart_interval) == 0
                   optimalChart = [optimalChart;bsf_fit_var];
                else
                    if nFES == FES
                       optimalChart = [optimalChart;bsf_fit_var];
                    end
                end
                if nFES > FES; break; end
            end

else

    if tier==1
       tier=2;
       ppsize=pop_size-popuSize;
       pop=popu;
       [~, indtest] = sort(fit, 'ascend');

       dsx=abs(repmat(pop(indtest(1),:),popuSize,1)-pop);
       bs=max(dsx,[],2);
       bb=max(bs,[],1);

       if bb>0.5*ub(1,1)
          bb=0.5*ub(1,1);
       elseif bb<0.01*ub(1,1)
          bb=0.01*ub(1,1);
       end
       bba=(rand(ppsize,D)*(bb+bb)-bb);
       pops=repmat(pop(indtest(1),:),ppsize,1)+bba;

       pops=BoundaryDetection(pops,lu);
       fitp=Evaluation(pops, benchmark, f);
       fitp=fitp';

    for i=1:ppsize
        nFES = nFES + 1;
        if min(fitp) < bsf_fit_var
           bsf_fit_var = min(fitp);
        end
        if mod(nFES, optimalChart_interval) == 0
           optimalChart = [optimalChart;bsf_fit_var];
        else
            if nFES == FES
               optimalChart = [optimalChart;bsf_fit_var];
            end
        end
        if nFES > FES; break; end
    end

    pop=[pop;pops];
    popold=pop;
    fit=[fit;fitp];

    end

if     ibbb>50
       ibbb=1;
       ibb=ibb+1;
       ppsize=pop_size-5;
       pop=popu;
       [~, indtest] = sort(fit, 'ascend');

       dsx=abs(repmat(pop(indtest(1),:),pop_size,1)-pop);
       bs=max(dsx,[],2);
       [bb,idxx]=max(bs,[],1);
       if bb>0.5*ub(1,1)*b
          bb=0.5*ub(1,1)*b;
       elseif bb<0.01*ub(1,1)*b
          bb=0.01*ub(1,1)*b;
       else
          bb=bb*b;
       end
       bba=(rand(ppsize,D)*(bb+bb)-bb);
       asize=round(0.3*ppsize);
       bsize=ppsize-2*asize;
       popsa=repmat(pop(indtest(1:2),:),asize,1)+0.25*bba(1:2*asize,D);
       popss=repmat(pop(idxx,:),bsize,1)+0.75*bba(2*asize+1:ppsize,D);
       pops=[popsa;popss];
       pops=BoundaryDetection(pops,lu);
       fitp=Evaluation(pops, benchmark, f);
       fitp=fitp';
       archive = updateArchive(archive,pop(indtest(6:pop_size),:),fit(indtest(6:pop_size)));

    for i=1:ppsize
        nFES = nFES + 1;
        if min(fitp) < bsf_fit_var
           bsf_fit_var = min(fitp);
        end
        if mod(nFES, optimalChart_interval) == 0
           optimalChart = [optimalChart;bsf_fit_var];
        else
            if nFES == FES
               optimalChart = [optimalChart;bsf_fit_var];
            end
        end
        if nFES > FES; break; end
    end

      [~,IX]= unique(pop, 'rows');
      if length(IX) < size(pop, 1) && length(IX)>4
         pop = pop(IX, :);
         fit = fit(IX, :);
      end
      [~, indtest] = sort(fit, 'ascend');

    pop=[pop(indtest(1:5),:);pops];
    popold=pop;
    fit=[fit(indtest(1:5),:);fitp];

end

%LSHADE%

      popu = popold;
      [~, sorted_index] = sort(fit, 'ascend');

      mem_rand_index = ceil(memory_size * rand(pop_size, 1));
      mu_sf = memory_sf(mem_rand_index);
      mu_cr = memory_cr(mem_rand_index);

      %% for generating crossover rate
      cr = normrnd(mu_cr, 0.1);
      term_pos = mu_cr == -1;
      cr(term_pos) = 0;
      cr = min(cr, 1);
      cr = max(cr, 0);

      %% for generating scaling factor
      sf = mu_sf + 0.1 * tan(pi * (rand(pop_size, 1) - 0.5));
      pos = find(sf <= 0);

      while ~ isempty(pos)
        	sf(pos) = mu_sf(pos) + 0.1 * tan(pi * (rand(length(pos), 1) - 0.5));
	        pos = find(sf <= 0);
      end

      sf = min(sf, 1); 
      
      r0 = 1 : pop_size;
      popAll = [popu; archive.popu];
      [r1, r2] = gnR1R2(pop_size, size(popAll, 1), r0);    
      pNP = max(round(p_best_rate * pop_size), 2);
      randindex = ceil(rand(1, pop_size) .* pNP);
      randindex = max(1, randindex);
      pbest = popu(sorted_index(randindex), :);
      vi = popu + sf(:, ones(1, D)) .* (pbest - popu + popu(r1, :) - popAll(r2, :));
      vi = BoundaryDetection(vi, lu);
      mask = rand(pop_size, D) > cr(:, ones(1, D));
      rows = (1 : pop_size)'; cols = floor(rand(pop_size, 1) * D)+1;
      jrand = sub2ind([pop_size D], rows, cols); mask(jrand) = false;
      ui = vi; ui(mask) = popu(mask);

      children_fitness = Evaluation(ui, benchmark, f);
      children_fitness = children_fitness';

            % get bsf_fit_var
            for i=1:pop_size
            nFES = nFES + 1;
            if min(children_fitness) < bsf_fit_var
               bsf_fit_var = min(children_fitness);
            end
            if mod(nFES, optimalChart_interval) == 0
               optimalChart = [optimalChart;bsf_fit_var];
            else
                if nFES == FES
                   optimalChart = [optimalChart;bsf_fit_var];
                end
            end
            if nFES > FES; break; end
            end

            dif = abs(fit - children_fitness);

      I1 = (fit > children_fitness);
      goodCR = cr(I1 == 1);  
      goodF = sf(I1 == 1);
      dif_val = dif(I1 == 1);

      archive = updateArchive(archive, popold(I1 == 1, :), fit(I1 == 1));

      [fit, I1] = min([fit, children_fitness], [], 2);
      
      popold = popu;
      popold(I1 == 2, :) = ui(I1 == 2, :);
      num_success_params = numel(goodCR);

      if num_success_params > 0 
	sum_dif = sum(dif_val);
	dif_val = dif_val / sum_dif;
	memory_sf(memory_pos) = (dif_val' * (goodF .^ 2)) / (dif_val' * goodF);
if max(goodCR) == 0 || memory_cr(memory_pos)  == -1
	memory_cr(memory_pos)  = -1;
else
	memory_cr(memory_pos) = (dif_val' * (goodCR .^ 2)) / (dif_val' * goodCR);
end
	memory_pos = memory_pos + 1;
	if memory_pos > memory_size;  memory_pos = 1;
    end
      end
   plan_pop_size = round((((min_pop_size - max_pop_size) / FES) * nFES) + max_pop_size);
if pop_size > plan_pop_size
	reduction_ind_num = pop_size - plan_pop_size;
	if pop_size - reduction_ind_num <  min_pop_size; reduction_ind_num = pop_size - min_pop_size;
    end
	pop_size = pop_size - reduction_ind_num;
for r = 1 : reduction_ind_num
	[~,indBest] = sort(fit, 'ascend');
	worst_ind = indBest(end);
	popold(worst_ind,:) = [];
	popu(worst_ind,:) = [];
	fit(worst_ind,:) = [];
end
archive.NP = round(arc_rate * pop_size);
if size(archive.popu, 1) > archive.NP 
   rndpos = randperm(size(archive.popu, 1));
   rndpos = rndpos(1 : archive.NP);
   archive.popu = archive.popu(rndpos, :);
end
end

      if deit==1
         fitde=min(fit);
      end
      if deit>=round(40*nFES/FES)+10
         deit=0;
         if (0.95+0.04*(1-b))*fitde<min(fit)
            ibbb=ibbb+1;
         end
      end
      deit=deit+1;
end

            
            %% Survival
            optimal = bsf_fit_var;
            fprintf('%3.8f %3.0f problem %5.0f time %5.0f |%5.0f -----> %9.30f\n',bb,ibb,f,run_id,nFES,optimal);
end


            while length(optimalChart(:,1))>FES/optimalChart_interval
                  optimalChart=optimalChart(2:popuSize*D+1,:);
            end
            while length(optimalChart(:,1))<FES/optimalChart_interval
                  optimalChart=[optimalChart;bsf_fit_var];
            end

            output.popu=popu;
            output.bsf_fit_var=bsf_fit_var;
            output.optimalChart=optimalChart;

function [r1, r2] = gnR1R2(NP1, NP2, r0)
NP0 = length(r0);
r1 = floor(rand(1, NP0) * NP1) + 1;
for i = 1 : 99999999
    pos = (r1 == r0);
    if sum(pos) == 0
        break;
    else
        r1(pos) = floor(rand(1, sum(pos)) * NP1) + 1;
    end
    if i > 1000
        error('Can not genrate r1 in 1000 iterations');
    end
end
r2 = floor(rand(1, NP0) * NP2) + 1;
for i = 1 : 99999999
    pos = ((r2 == r1) | (r2 == r0));
    if sum(pos)==0
        break;
    else
        r2(pos) = floor(rand(1, sum(pos)) * NP2) + 1;
    end
    if i > 1000
        error('Can not genrate r2 in 1000 iterations');
    end
end

function archive = updateArchive(archive, popu, funvalue)
if archive.NP == 0, return;
end
if size(popu, 1) ~= size(funvalue,1), error('check it');
end
popAll = [archive.popu; popu];
funvalues = [archive.funvalues; funvalue ];
[~,IX]= unique(popAll, 'rows');
if length(IX) < size(popAll, 1)
  popAll = popAll(IX, :);
  funvalues = funvalues(IX, :);
end

if size(popAll, 1) <= archive.NP
  archive.popu = popAll;
  archive.funvalues = funvalues;
else
  rndpos = randperm(size(popAll, 1));
  rndpos = rndpos(1 : archive.NP);
  archive.popu = popAll  (rndpos, :);
  archive.funvalues = funvalues(rndpos, :);
end

function ss = HyperSphereTransform(c,d,pp)
D=length(pp);
A=c(pp)-d(pp);
R=norm(A,2);
O(D-1)=  2*pi*rand ;
for i=1:D-2
    O(i)=  rand*pi  ;
end

% 转换直角坐眮E
C(1)=R*prod(sin(O));
for i=2:D-1
    C(i)= R*cos(O(i-1)) *prod(sin(O(i:D-1)));
end
C(D)=R*cos(O(D-1));
ss=C;

function ss = HyperSphereTransform_1D(c,d,pp)
R=  abs(c(pp)-d(pp));
C = R*cos(2*pi*rand);
ss=C;

function ss = HyperSphereTransform_2D(c,d,pp)
A=c(pp)-d(pp);
R=norm(A,2);
o1=2*pi*rand;
C=zeros(1,2);
C(1) =R*sin(o1);
C(2) =R*cos(o1);
ss=C;