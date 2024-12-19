function increm_eval(func_num, fitness)
  global initial_flag
  global bench_output_name
  global bench_bestfit
  persistent numevals
  persistent maxevals
  persistent bench_milestones
  persistent next_milestone
  persistent next_milestone_i
  persistent fid

  % initial_flag == 1 if benchmark is initialized, == 2 when this code is done
  if (initial_flag < 2)
    % Default output name
    if isempty(bench_output_name)
      bench_output_name = sprintf('results_%d.csv', func_num);
    end

    fid = func_num;
    maxevals = 3e6;
    numevals = 0;
    increm = maxevals/10;
    values = colon(increm, increm, maxevals);
    values = cat(2, values, [1.2e5, 6e5, 3e6]);
    bench_milestones = unique(sort(values));
    next_milestone = bench_milestones(1);
    next_milestone_i = 1;
    initial_flag = 2;
  end

  current_best = min(fitness);
  [dim, pop] = size(fitness);
  numevals = numevals+pop*dim;

  if numevals <= pop*dim
    bench_bestfit = current_best;
  elseif current_best < bench_bestfit
    bench_bestfit = current_best;
  end

  % Check the numevals is allowed
  if numevals <= maxevals
    if numevals == next_milestone
      file = fopen(bench_output_name, 'at');
      fprintf(file, '%d,%d,%e\n', numevals, fid, bench_bestfit);
      fclose(file);
      next_milestone_i = next_milestone_i+1;
      next_milestone = bench_milestones(next_milestone_i);
    end
  end
 
end
