function [X,FVAL,EXITFLAG] = pso(EvalFun, X0, LB, UB, problem_data)

% X=pso(@FUN,X0,LB,UB) returns the best solution X obtained by the particle
% swarm optimization of the function FUN with the default parameters
% assigned by psoset. XO is an n X 1 or 1 X n array provided to make an
% initial guess.  Empty brackets,[], may be used for XO if an initial guess
% is not desired. LB and UB define a set of lower and upper bounds on the
% design variables, X, so that a solution is found in the range LB <= X
% <=UB. The dimensions of LB and UB must be the same as those of X and
% therefore must be the same as the dimensonality of the function FUN and
% they may be either n X 1 or 1 X n arrays. 

% X=pso(@FUN,X0,LB,UB,OPTIONS) minimizes with the default
% optimization parameters replaced by values in the structure OPTIONS, an
% argument created with the PSOSET function.  

% [X,FVAL]=pso(@FUN,X0...) additionally returns FVAL, the value of the
% objective function FUN at the solution X. 

% [X,FVAL,EXITFLAG]=pso(@FUN,X0,...) returns an EXITFLAG that describes the
% exit condition of the PSO algorithm. Possible values of EXITFLAG and the
% corresponding exit conditions are:

%       0  - Maximum number of function evaluations or iterations reached. 
%       1  - Distance between current and previous particle location in
%            design space less than the specified tolerance for every
%            particle.
%       2  - Absolute difference between current and previous objective
%            function value less than the specified tolerance for every
%            particle. 
%       3  - Distance between centroid of all particles and each particle
%            less than the specified tolerance for every particle.
%       4  - The maximum time specified in OPTIONS was reached.

% Options for PSO:
% psoset Create/alter PSO OPTIONS structure.
% 
% OPTIONS = psoset('PARAM1',VALUE1,'PARAM2',VALUE2,...) creates a PSO
% optimization options structure OPTIONS in which the named parameters have
% the specified values.  Any unspecified parameters are set to their
% default settings. Case is ignored. 

% Type psoset to see the options structure with the {default} options
% encloded in brackets.

% Type help psoset to obtain an explanation of each of the options.

startTime = clock;

% error check for the correct number of arguments
if nargin < 4,
    error('There are not enough arguments for the pso function. The minimum required arguments are pso(@EvalFun, X0, LB, UB)');
elseif nargin == 4,
    problem_data = psoset;
elseif nargin > 5,
    error('There are too many input arguments for the pso function.');
else %do nothing    
end


%ensure that the boundary conditions will work regardless of orientation
%(n X 1 or 1 X n).  Also sets the problems dimensions from the lower
%boundary condition.
[m,n] = size(LB);
[o,p] = size(UB);
if (p > o),
    UB=UB';
end
if (m > n) & (n ==1),
    problem_data.problem_dimensions = m;
elseif (n > m) & (m == 1),
    problem_data.problem_dimensions = n;
    LB=LB';
else
    error('The boundary conditions are not entered correctly');
end


%after reorienting the boundary conditions check to ensure the same
%dimensionality.
[m,n] = size(LB);
[o,p] = size(UB);
if ~(m==o),
    error('The boundaries conditions do not have the same dimensions');
end

%error check to make sure that the lower boundaries are indeed lower than
%the upper boundaries.
if LB < UB,
    %do nothing
else
    error('The Lower Bounds are larger than the Upper Bounds');
    %doing it this way causes the program to throw an error if any of the
    %lower bounds are greater than the upper bounds.
end



X = zeros(problem_data.problem_dimensions,1); %preallocate memory for the solution.
FVAL = 0; %preallocate memory for the solution.

inertreduce=0; %a troubleshooting parameter.  Tracks the number of times the inertia is reduced.

rand('state',sum(100*clock)) %randomly seeds the RN generator with the clock time.
problem_data.position_upperbound = UB;
problem_data.position_lowerbound = LB;
varphi = problem_data.SocialWeight + problem_data.CognitiveWeight; %a variable for constriction
%various flags and counters
solution_found_flag = 0; %
plotcounter=1;
no_improvement_counter = 0;
iter_counter = 1;

waitbarOn = strcmp(problem_data.WaitBar,'on');
if waitbarOn,
    wb = waitbar(0,'Percentage of max function evaluations completed');
end

% boolean calculations to determine the options selected
displayOnIter = strcmp(problem_data.Display,'iter');
displayOnEnd = strcmp(problem_data.Display,'on');
plothistory = strcmp(problem_data.Plot,'history');
plotoutput = strcmp(problem_data.Plot,'output');
plotdata = strcmp(problem_data.Plot,'designspace');
async = strcmp(problem_data.AlgorithmType,'asynchronous');
dynamic = strcmp(problem_data.ModificationMethod,'dynamic');
lhs = strcmp(problem_data.InitializationMethod,'lhs');

if problem_data.problem_dimensions == 2,
    plot2d=1;
else
    plot3d=1;
end


% histogram calculation not implemented in final copy, but left if desired
% later
% binEl = 100;
% binLength = (UB(1)-LB(1))/binEl;
% histoMatrix = zeros(binEl,binEl);

%preallocate memory for centroid calculation and tolfun and tolX
%calculations
tempdist2cen=zeros(problem_data.problem_dimensions,1);
dist2cen=zeros(problem_data.NumParticles,1);
Xnew = zeros(problem_data.problem_dimensions,problem_data.NumParticles);
Xold = Xnew;
Fnew = zeros(problem_data.NumParticles,1);
Fold = Fnew;


%a structure to contain iteration history
iteration_history = repmat(struct('FunEvals',0,'Objective',0,'MaxFun',0,'MaxX',0,'Centroid',...
    zeros(problem_data.problem_dimensions,1),'dist2cen',0),1,...
    ceil((problem_data.MaxFunEvals/problem_data.NumParticles)));


%Structures
%particle_data structure
particle = repmat(struct('coordinates', zeros(problem_data.problem_dimensions,1),'velocity',...
    zeros(problem_data.problem_dimensions,1),'inertia',problem_data.InitialInertia,'current_value',0,'best_value',0,...
    'normalized_infeasability',0,'best_value_coordinates',zeros(problem_data.problem_dimensions,1),...
    'type_of',0,'constraints_violated',0,'previous_value',0),1,problem_data.NumParticles);

%communication_data structure (contains global information)
comm_data.swarm_best_value = 0;
comm_data.swarm_best_coordinates = zeros(problem_data.problem_dimensions,1);
comm_data.swarm_repulse_value = 0;
comm_data.swarm_initial_value = 0;
comm_data.swarm_max_velocity = zeros(problem_data.problem_dimensions,1);
comm_data.swarm_best_normalized_infeasability = 0;
comm_data.function_evaluations = double(0);
comm_data.times_converged = double(0);
comm_data.tolerance = 0;

%function_result structure
function_result.value = 0;			%/* function value at position */
function_result.normalized_infeasability = 0;



%initialize particle positions
%/* this function initializes all the particle positions */

randomfloat = 0;
range = zeros(problem_data.problem_dimensions,1);

for ipp_j=1:problem_data.problem_dimensions, % sets range
    range(ipp_j) = problem_data.position_upperbound(ipp_j) - problem_data.position_lowerbound(ipp_j);
end

if ~isempty(X0), %used to make X0 particle 1 if X0 is entered.
    ipp_i = 2;
    j = 1;
    particle(1).coordinates = X0;
else
    ipp_i = 1;
end
if lhs, %generates Latin Hypercube Sample
    lhsRand =lhsdesign(problem_data.NumParticles,problem_data.problem_dimensions,'criterion','maximin');
end
while ipp_i <= problem_data.NumParticles,  %creates initial particle positions.
    ipp_j = 1;
    while ipp_j <= problem_data.problem_dimensions,
        if lhs,
            particle(ipp_i).coordinates(ipp_j) = lhsRand(ipp_i,ipp_j)*range(ipp_j) + problem_data.position_lowerbound(ipp_j);
        else
            randomfloat = rand;
		    particle(ipp_i).coordinates(ipp_j) = randomfloat*range(ipp_j) + problem_data.position_lowerbound(ipp_j);
        end
        ipp_j = ipp_j + 1;
    end %for ipp_j
    ipp_i = ipp_i + 1;
end %for ipp_i

%initialize particle velocities
%/* this function initializes all the particle velocities */

velocity_range = zeros(problem_data.problem_dimensions,1);

	%/* Init velocity bounds */

for j=1:problem_data.problem_dimensions, %creates velocity range
    velocity_range(j) = problem_data.VelocityFraction * (problem_data.position_upperbound(j)-problem_data.position_lowerbound(j));

end

for i=1:problem_data.NumParticles, %randomly creates initial velocities
    for j=1:problem_data.problem_dimensions,
        randomfloat = rand;
        particle(i).velocity(j) = randomfloat * velocity_range(j) - (velocity_range(j) * 0.5);
    end %for j
end %for i

%initialize particle inertia values
%/* this function sets up all the particle inertia values */

for i=1:problem_data.NumParticles, %sets all particle inertia values to the InitialInertia
    particle(i).inertia = problem_data.InitialInertia;
end

%initialize particle function values
%/* This function initializes all the particle function values */
infeasable_particles = double(0);

for i=1:problem_data.NumParticles, %solve for function values
    [function_result.value, function_result.normalized_infeasability] = feval(EvalFun,particle(i).coordinates);
	comm_data.function_evaluations = comm_data.function_evaluations + 1;
	particle(i).current_value = function_result.value;
	particle(i).normalized_infeasability = function_result.normalized_infeasability;
	particle(i).best_value = particle(i).current_value;
end

for i=1:problem_data.NumParticles, %set each initialized value to best value
    for j=1:problem_data.problem_dimensions,
        particle(i).best_value_coordinates(j) = particle(i).coordinates(j);
    end
    
% 	/* initialize swarm best if NI < NI_tol, and search for better one */

    if particle(i).constraints_violated == 0,
        comm_data.swarm_best_value = particle(i).best_value;  %/* for temporary comparison */
		comm_data.swarm_best_normalized_infeasability = particle(i).normalized_infeasability;
    end
end %for particle i        

	%/* search for best particle fitness value */
for i=1:problem_data.NumParticles,
    if ( (particle(i).best_value < comm_data.swarm_best_value) & (particle(i).constraints_violated == 0) ),
        comm_data.swarm_best_value = particle(i).best_value;
	    comm_data.swarm_best_normalized_infeasability = particle(i).normalized_infeasability;
        for j=1:problem_data.problem_dimensions,
            comm_data.swarm_best_coordinates(j) = particle(i).best_value_coordinates(j);
        end % for j
    end % if

end % for i

% end initialize particle function values

% initialize swarm max velocity

for j=1:problem_data.problem_dimensions,
    comm_data.swarm_max_velocity(j) = problem_data.VelocityFraction *...
        (problem_data.position_upperbound(j)-problem_data.position_lowerbound(j));
end %for j

% end initialize swarm max velocity

% begin swarm iterations
while (solution_found_flag == 0),

    % UPDATE FUNCTIONS
    % Velocity update
    improvement_flag=0;
    for i=1:problem_data.NumParticles, %loop through the number of particles to update each individually
        if (dynamic), %dynamic reduction
            for j=1:problem_data.problem_dimensions,
                    randomdouble1 = rand;
                    randomdouble2 = rand;
                    particle(i).velocity(j) = particle(i).inertia * particle(i).velocity(j) + problem_data.CognitiveWeight*randomdouble1*...
                        (particle(i).best_value_coordinates(j) - particle(i).coordinates(j)) +...
                        problem_data.SocialWeight*randomdouble2*(comm_data.swarm_best_coordinates(j)...
                        - particle(i).coordinates(j));
                    
            end %for j
        else %constriction
            randomdouble1 = rand;
            randomdouble2 = rand;
            for j=1:problem_data.problem_dimensions,
                constriction_factor = 2/abs(2 - varphi - sqrt(varphi*varphi - 4 * varphi));
                particle(i).velocity = constriction_factor*(particle(i).velocity) + problem_data.CognitiveWeight*randomdouble1*...
                    (particle(i).best_value_coordinates(j) - particle(i).coordinates(j)) +...
                    problem_data.SocialWeight*randomdouble2*(comm_data.swarm_best_coordinates(j)...
                    - particle(i).coordinates(j));                
            end %for j
        end %if
        
        if (problem_data.LimitMaxVelocity),

            for j=1:problem_data.problem_dimensions,
                if particle(i).velocity(j) > comm_data.swarm_max_velocity(j),
                   particle(i).velocity(j) = comm_data.swarm_max_velocity(j);
                end
                if particle(i).velocity(j) < -comm_data.swarm_max_velocity(j),
                   particle(i).velocity(j) = -comm_data.swarm_max_velocity(j);
                end
            end
        end
        %i = i + 1;
        %end %for i / end velocity update
    
    %position update

        for j =1:problem_data.problem_dimensions,
            particle(i).coordinates(j) = particle(i).coordinates(j) + particle(i).velocity(j);

            if (particle(i).coordinates(j) > problem_data.position_upperbound(j)) || (particle(i).coordinates(j) < problem_data.position_lowerbound(j)),
                
                switch problem_data.BoundsMethod,
                    case {'bounce'}
                        particle(i).velocity(j) = -particle(i).velocity(j);
                        particle(i).coordinates(j) = particle(i).coordinates(j) + particle(i).velocity(j);

                    case {'rand'}
                        for ipp_j=1:problem_data.problem_dimensions,
                            randomfloat = rand;
                            particle(i).coordinates(ipp_j) = randomfloat*range(ipp_j)...
                                + problem_data.position_lowerbound(ipp_j);
                        end %for ipp_j
                        break; 
                        %since the random respawn resets all dimensions,
                        %this break allows the position to be completely
                        %reset and to end the main 'for' loop from looking
                        %through the remaining dimensions to update position.
                    case {'drift'}
                        %do nothing and allow particles to drift out of
                        %bounds
                    otherwise
                        error('There is an error in the BoundsMethod setting.  The value must be: bounce|rand|drift');
                end % switch
            end % if
                
        end %for j

    
    % update particle function value
        [function_result.value, function_result.normalized_infeasability] = feval(EvalFun,particle(i).coordinates);
        comm_data.function_evaluations = comm_data.function_evaluations + 1;
        particle(i).previous_value = particle(i).current_value;
        particle(i).current_value = function_result.value;
        particle(i).normalized_infeasability = function_result.normalized_infeasability;
        
    %update best function value after each function evaluation if asynchronous
        if async,
            if (particle(i).current_value < particle(i).best_value) & (particle(i).constraints_violated == 0),
                particle(i).best_value = particle(i).current_value;
                for j=1:problem_data.problem_dimensions,
                    particle(i).best_value_coordinates(j) = particle(i).coordinates(j);
                end %for j        
            end % if
        
            if (particle(i).current_value < comm_data.swarm_best_value) & (particle(i).constraints_violated == 0),
                comm_data.swarm_best_value = particle(i).current_value;
                comm_data.swarm_normalized_infeasability = particle(i).normalized_infeasability;
                for j=1:problem_data.problem_dimensions,
                    comm_data.swarm_best_coordinates(j) = particle(i).coordinates(j);
                end %for j
                no_improvement_counter = 0;
                improvement_flag=1;
            elseif i==problem_data.NumParticles && ~improvement_flag,
                no_improvement_counter = no_improvement_counter + 1;
            else
                %do nothing
            end % if
        end %async search
        
        % plot the best found value of EvalFun if the options are selected
        % to do so
        if plothistory || plotoutput,
            figure(1);
            hold on
            plot(comm_data.function_evaluations,comm_data.swarm_best_value,'-');
            hold off
            if plothistory && rem(comm_data.function_evaluations,problem_data.IterationsNoImprovement*problem_data.NumParticles) == 0,
                drawnow;
            end
        end
        
        % plot the designspace if selected
        if (i == problem_data.NumParticles) & (plotdata),
	        figure(1);
            for j=1:problem_data.NumParticles,
                if (plot3d),
                    plot3(particle(j).coordinates(1),particle(j).coordinates(2),particle(j).coordinates(3),'bo');
                    hold on;
                    plot3(particle(j).best_value_coordinates(1),particle(j).best_value_coordinates(2),particle(j).best_value_coordinates(3),'y+');
                    axis([LB(1) UB(1) LB(2) UB(2) LB(3) UB(3)]);
                elseif (plot2d),
                    plot(particle(j).coordinates(1),particle(j).coordinates(2),'bo');
                    hold on;
                    plot(particle(j).best_value_coordinates(1),particle(j).best_value_coordinates(2),'y+');
                    axis([LB(1) UB(1) LB(2) UB(2)]);
                end
            end
             if (plot3d),
                 plot3(comm_data.swarm_best_coordinates(1),comm_data.swarm_best_coordinates(2),comm_data.swarm_best_coordinates(3),'rx');
             else
                 plot(comm_data.swarm_best_coordinates(1),comm_data.swarm_best_coordinates(2),'rx');
                 if iter_counter > 1,
                     plot(iteration_history(iter_counter-1).Centroid(1),iteration_history(iter_counter-1).Centroid(2),'k+');
                 end
             end
            hold off;
            drawnow;
        end %if
        
        % update waitbar if on.
        if waitbarOn,
            waitbar(comm_data.function_evaluations/problem_data.MaxFunEvals,wb);
        end
        
        % end if the function evaluatons have exceeded the maximum
        if (comm_data.function_evaluations >= problem_data.MaxFunEvals),
            solution_found_flag = 1;
            %EXITFLAG=comm_data.function_evaluations;
            EXITFLAG = 0;
            %break;
        end
       
	end %for particle loop

    %Update global best found after the swarm iteration if synchronous
	if ~async,
        for i=1:problem_data.NumParticles,
                if (particle(i).current_value < particle(i).best_value) & (particle(i).constraints_violated == 0),
                    particle(i).best_value = particle(i).current_value;
                    for j=1:problem_data.problem_dimensions,
                        particle(i).best_value_coordinates(j) = particle(i).coordinates(j);
                    end %for j        
                end % if
            
                if (particle(i).current_value < comm_data.swarm_best_value) & (particle(i).constraints_violated == 0),
                    comm_data.swarm_best_value = particle(i).current_value;
                    comm_data.swarm_normalized_infeasability = particle(i).normalized_infeasability;
                    for j=1:problem_data.problem_dimensions,
                        comm_data.swarm_best_coordinates(j) = particle(i).coordinates(j);
                    end %for j
                    no_improvement_counter = 0;
                    improvement_flag=1;
                elseif i==problem_data.NumParticles && ~improvement_flag,
                    no_improvement_counter = no_improvement_counter + 1;
                else
                    %do nothing
                end % if
                
        end %for particle loop
	end %sync search

        % update inertia values if no improvement has occured
	if no_improvement_counter >= problem_data.IterationsNoImprovement,
        inertreduce=inertreduce+1;
        for i=1:problem_data.NumParticles,
            particle(i).inertia = particle(i).inertia*(1-problem_data.InertiaReductionFraction);
        end
        for j=1:problem_data.problem_dimensions,
            comm_data.swarm_max_velocity(j)= comm_data.swarm_max_velocity(j) *(1-problem_data.VelocityReductionFraction);
        end
            
            %fprintf('Inertia Reduced %d times\n',inertreduce);
        no_improvement_counter = 0;
    end
	
    %calculations for tolerances
	Xold = Xnew;
	Fold = Fnew;
	for i=1:problem_data.NumParticles,
        Xnew(:,i)=particle(i).coordinates;
        Fnew(i)=particle(i).current_value;
	end

%2Dhistogram
% for i=1:problem_data.NumParticles,
% 
%     x1 = floor((particle(i).coordinates(1)-LB(1))/binLength)+1;
%     y1 = binEl-floor((particle(i).coordinates(2)-LB(2))/binLength);
%     histoMatrix(y1,x1) = histoMatrix(y1,x1) + 1;
% end
%2Dhistogram


%centroid calculation
    tempvar=zeros(problem_data.problem_dimensions,1);
    for j=1:problem_data.problem_dimensions,
        for i=1:problem_data.NumParticles
            tempvar(j)= tempvar(j) + particle(i).coordinates(j);
        end
        tempvar(j) = tempvar(j)/problem_data.NumParticles;
        iteration_history(iter_counter).Centroid(j)=tempvar(j);
    end

%iteration history update
    iteration_history(iter_counter).FunEvals = comm_data.function_evaluations;
    iteration_history(iter_counter).Objective = particle(1).current_value;
    iteration_history(iter_counter).MaxFun = abs(particle(1).current_value - particle(1).previous_value);
    iteration_history(iter_counter).MaxX = norm(particle(1).velocity);
    for i =1:problem_data.NumParticles,
        for j=1:problem_data.problem_dimensions,
            tempdist2cen(j) = iteration_history(iter_counter).Centroid(j) - particle(i).coordinates(j);
        end
        dist2cen(i) = norm(tempdist2cen);
    end
    iteration_history(iter_counter).dist2cen = dist2cen(1);

    for j=1:problem_data.NumParticles,
        if iteration_history(iter_counter).Objective >= particle(j).current_value & (particle(j).constraints_violated == 0),
            iteration_history(iter_counter).Objective = particle(j).current_value;
        end
        if iteration_history(iter_counter).MaxFun <= abs(particle(j).current_value - particle(j).previous_value) & (particle(j).constraints_violated == 0);
            iteration_history(iter_counter).MaxFun = abs(particle(j).current_value - particle(j).previous_value);
        end
        if iteration_history(iter_counter).MaxX <= norm(particle(j).velocity) & (particle(j).constraints_violated == 0);
            iteration_history(iter_counter).MaxX = norm(particle(j).velocity);
        end
        if iteration_history(iter_counter).dist2cen <= dist2cen(i),
            iteration_history(iter_counter).dist2cen = dist2cen(i);
        end
    end
    %end

    %check TolFun tolerance
    if~(abs(Fnew) < 2*eps), %this statement will not allow division by zero.  If zero is the minimum value the program will not end itself on the TolFun criteria.
        if (abs((Fnew-Fold)./Fnew) < problem_data.TolFun),
            solution_found_flag = 1;
            if displayOnIter,
                fprintf('A solution was reached based on the TolFun criteria in the following iteration:\n');
            end
            if waitbarOn,
               waitbar(1,wb);
            end
            EXITFLAG = 2;
        end
    end
    
    %check TolX tolerance
    if~(abs(Xnew)<2*eps),%this statement will not allow division by zero.  If zero is the location of the minimum value the program will evaluate the stopping condition on absolute tolerance only.
        if (abs(Xnew-Xold) < problem_data.TolX),
            solution_found_flag = 1;
            if displayOnIter,
                fprintf('A solution was reached based on the TolX criteria in the following iteration:\n');
            end
            if waitbarOn,
                waitbar(1,wb);
            end
            EXITFLAG = 1;
        end        
    else
        if (abs(Xnew-Xold) < problem_data.TolX | abs((Xnew-Xold)./Xnew) < problem_data.TolX),
            solution_found_flag = 1;
            if displayOnIter,
                fprintf('A solution was reached based on the TolX criteria in the following iteration:\n');
            end
            if waitbarOn,
                waitbar(1,wb);
            end
            EXITFLAG = 1;
        end
    end
    
    %check TolCen tolerance
    if (problem_data.TolCen > iteration_history(iter_counter).dist2cen),
        solution_found_flag = 1;
        if displayOnIter,
            fprintf('A solution was reached based on the TolCen criteria in the following iteration:\n');
        end
        if waitbarOn,
            waitbar(1,wb);
        end
        EXITFLAG = 3;
    end

    %%the below code is used to return early without doing all evaluations
    %%if the answer is previously known. (testing purposes
%     if(comm_data.swarm_best_value < 0 + 0.001),
%         solution_found_flag = 1;
%         EXITFLAG=comm_data.function_evaluations;
%     end
    
    %% end special code
    
    %check time limit
    if (etime(clock,startTime) > problem_data.MaxTime),
        solution_found_flag = 1;
        fprintf('A solution was reached based on the MaxTime (%dsec) criteria in the following iteration:\n',problem_data.MaxTime);
        if waitbarOn,
            waitbar(1,wb);
        end
    end
    
    %display Iteration History
    if (displayOnIter),
        if (iter_counter == 1),
            fprintf('Iteration      Function Evals      Objective            MaxFun                  MaxX         MaxCen\n');
        end
            fprintf('  %.5d            %.5d         %c      %d         %d        %d\n',...
                                iter_counter,iteration_history(iter_counter).FunEvals,comm_data.swarm_best_value,iteration_history(iter_counter).MaxFun,...
                iteration_history(iter_counter).MaxX,iteration_history(iter_counter).dist2cen);
                %iter_counter,iteration_history(iter_counter).FunEvals,iteration_history(iter_counter).Objective,iteration_history(iter_counter).MaxFun,...
                %iteration_history(iter_counter).MaxX,iteration_history(iter_counter).dist2cen);


    end
    iter_counter = iter_counter + 1;
end %while solution not found

    %display iteration history at end
if (displayOnEnd),
        fprintf('Iteration      Function Evals      Objective            MaxFun                  MaxX\n');
        totalIterations = iter_counter-1;

        for iter_counter=1:totalIterations,
            fprintf('  %.5d            %.5d         %c      %d         %d\n',iter_counter,iteration_history(iter_counter).FunEvals,iteration_history(iter_counter).Objective,iteration_history(iter_counter).MaxFun,iteration_history(iter_counter).MaxX);
        end
end
    
%assign values to output
FVAL = comm_data.swarm_best_value;
X = comm_data.swarm_best_coordinates;

%close the waitbar
if waitbarOn,
    close(wb);
end

%Troubleshooting aids:
%particleOneInert=particle(1).inertia
%InertReduced=inertreduce
% pcolor(histoMatrix);
% axis square;

 
