function [xbest, fbest, hist] = PSO(n_particle, max_iter,lb,ub, func, options)
%
% Implements Particle swarm optimization. The following steps follows
% Algorithm 7.6 (page 318) of the Martins and Ning (2021).
%
% Joaquim R. R. A. Martins and Andrew Ning. Engineering Design Optimization.
% Cambridge University Press, 2021. ISBN: 9781108833417.
% 

% unpack
beta_max = options.beta_max;
gamma_max = options.gamma_max;
max_vel = options.max_vel;

params = func.params;
calc_obj = func.fobj;

% inertia parameters
alfa_upper = max(options.alfa); 
alfa_lower = min(options.alfa);
dalfa = alfa_upper - alfa_lower;

x_size = length(ub);         % number of design variables
x = zeros(n_particle,x_size);
vel = zeros(n_particle,x_size);
hist = struct;

k = 0;
% initialize all particles
for i=1:n_particle
    % generate postion
    x(i,:) = rand(1,x_size).*(ub - lb) + lb;
    % initialize velocity
    vel(i,:) = rand(1,x_size).*max_vel;

    % compute objective
    obj = calc_obj(x(i,:),params);
    % apply penalty method for constraint handling
    fitness(i) = obj;

    % update personal best/pbest
    pbest(i,:) = x(i,:);
    pbest_fitness(i,:) = fitness(i);
end

while true
    
    % select best individual
    gbest_fitness = min(pbest_fitness);
    
    for i=1:n_particle
        % compute pbest or best individual points
        if fitness(i)<= pbest_fitness(i)
            pbest(i,:) = x(i,:);
            pbest_fitness(i) = fitness(i);
        end
        % compute gbest or best swarm point
        if fitness(i) <= gbest_fitness
            xbest = pbest(i,:);
            gbest_fitness = fitness(i);
        end

    end   
    
    for i=1:n_particle
        
        alfa = alfa_upper - (k/max_iter)*dalfa;
        beta = beta_max.*rand(1,x_size);
        gamma = gamma_max.*rand(1,x_size);
        
        % calculate velocity
        vel(i,:) = alfa.*vel(i,:) + beta.*(pbest(i,:) - x(i,:)) + ...
                                   gamma.*(xbest - x(i,:));

        % limit velocity
        vel(i,:) = max(min(vel(i,:),max_vel ), -max_vel);

        % update position
        x(i,:) = x(i,:) + vel(i,:);

        % enforce bounds
        x(i,:) = max(min(x(i,:),ub),lb);

        % compute objective
        obj = calc_obj(x(i,:),params);

        % apply penalty method for constraint handling
        fitness(i) = obj;

    end

    k=k+1;

    % history for post-process
    hist(k).x = x;
    hist(k).xbest = xbest;
    hist(k).pbest = pbest;

    if k>max_iter, break; end

end

fbest = gbest_fitness;
