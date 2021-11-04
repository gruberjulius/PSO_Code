function options = psoset(varargin)

% OPTIONS = psoset('PARAM1,VALUE1,'PARAM2',VALUE2,...) creates a PSO
% options structure OPTIONS in which the named parameters have the
% specified values.  Any unspecified parameters are set to their default
% settings and case is ignored. 
%
% Type 'psoset' to see the options structure with the {default} options
% encloded in brackets.
%
% Options for pso:
%
%                Display: [ off | on | {iter} ]
%
% For 'on', the convergence history only displays at the very end.
% 'iter' displays the following information after each iteration (slows
% down the optimization significantly):
%
% Iteration   Function Evals   Objective   MaxFun    MaxX   MaxCen
%
% MaxFun displays the absolute difference between the current and previous
% function value for each particle, and reports the maximum.
%
% MaxX displays the magnitude of the largest velocity vector from the
% current iteration. This is the maximum distance between the current and
% previous particle position among all particles.
%
% MaxCen displays the magnitude of the largest distance between the
% centroid of all particles and the location of each particle.
%
%                 WaitBar: [ {off} | on ]
%
% For 'on', a waitbar is displayed that constantly updates the percentage
% of the maximum function evaluations that has been completed. 
%
%                   Plot: [ {off} | output | history | designspace ]
%
% 'output' displays a history plot of the best found value of FUN at the
% end of the optimization.
%
% 'history' displays a history plot of the best found value of FUN during
% the optimization (slows the calculation slightly).
%
% 'designspace' plots the first 2 or 3 design variables in the design
% space.  The global best found (x) and the best for for each particle (+)
% are also plotted in lighter colors.
% 
%           AlgorithmType: [ synchronous | {asynchronous} ]
%
% 'asynchronous' updates the global best found position after each function
% evaluation.
%
% 'synchronous' updates the global best found position after each swarm
% iteration.
%
%            NumParticles: [ positive scalar between 5 and 300] {20}
%
% NumParticles sets the number of particles in the swarm.
%
%    MaxFunEvals: [ positive scalar between NumParticles and 1e12]{10000}
%
% MaxFunEvals sets the limit of the maximum number of function evaluations.
%
%                  TolFun: [ positive scalar ] {1e-4}
%
% TolFun sets the stopping tolerance for the relative change in the value
% of FUN. TolFun is not calculated as the value of FUN approaches 0. Can be
% set to 0 to turn off this stopping criteria.
%
%                    TolX: [ positive scalar ] {1e-4}
%
% TolX sets the stopping tolerance of the relative and absolute change in
% X.  The pso will stop on whichever, relative or absolute, occurs first
% and only uses the absolute tolerance as X approaches 0. Can be set to 0
% to turn off this stopping criteria.
%
%                  TolCen: [ positive scalar ] {1e-4}
%
% TolCen sets the tolerance for the distance between each particle
% and the centroid of all the particles.  If all of the particles are
% closer to the centroid than TolCen, the program will stop execution. Can
% be set to 0 to turn off this stopping criteria.
%
%                 MaxTime: [ positive scalar ] {7200}
%
% MaxTime is the maximum time in seconds that the optimization will run.
% Can be set to inf to turn off this stopping criteria.
%
%    InitializationMethod: [ {rand} | lhs ]
%
% 'rand' randomly initializes the positions of the particles using the
% random number generator.
% 'lhs' initializes the positions of the particles based on a Latin
% Hypercube Sample (requires the Statistics Toolbox).
%
%            BoundsMethod: [ rand | {bounce} | drift ]
%
% 'rand' randomly respawns a particle if it moves beyond the boundaries.
% 'bounce' bounces the particles off of the boundaries.
% 'drift' allows particles to move beyond the boundaries.
%
%      ModificationMethod: [ constriction | {dynamic} ]
%
% 'dynamic' reduces the inertia of the particles dynamically as the
% optimization continues (superior method).
% 'constriction' reduces the inertia of the particles in a different manner
% than dynamic.
%
%          InitialPenalty: [ positive scalar between 0 and 1000] {100}
%
%             FinalPenalty: [ positive scalar between 0 and 1000] {1000}
%
%         PenaltyChangeEnd: [ positive scalar between 0 and MaxFunEvals] 
%
%           InitialInertia: [ positive scalar between 0 and 2] {1}
%
% This is the initial inertia of each of the particles.
%
%          CognitiveWeight: [ positive scalar between 0 and 5 ] {2}
%
% The CognitiveWeight is the weight given to the particles individual best
% found position.
%
%             SocialWeight: positive scalar between 0 and 5 ] {2}
%
% The SocialWeight is the weight given to the global best found position.
%
%  IterationsNoImprovement: [ positive scalar between 1 and MaxFunEvals]
%  {10}
%
% IterationsNoImprovement is the number of iterations without improvement
% at which the inertia and maximum velocity is constricted.
%
%         NoImprovementTol: [ positive scalar between 0 and inf] {0.01}
%
% NoImprovementTol is the tolerance on IterationsNoImprovement.
%
%         LimitMaxVelocity: [ binary 0 or 1 ] {1}
%
% 0 - the maximum velocity will not be limited.
% 1 - the maximum velocity will be limited and dynamically reduced.  It is
% initially calculated from the range between the UB and LB and multiplied
% by the VelocityFraction.
%
% InertiaReductionFraction: [ positive scalar between 0 and 1] {0.05}
%
% This is the amount by which the inertia is reduced when
% IterationsNoImprovement is reached.  0.05 represents 5% reduction.
%
% VelocityReductionFraction: [ positive scalar between 0 and 1] {0.05}
%
% This is the amount by which the maximum velocity is reduced when
% IterationsNoImprovement is reached.  0.05 represents 5% reduction.
%
%         VelocityFraction: [ positive scalar between 0 and 1] {0.50}
%
% VelocityFraction is multiplied by the range between the UB and LB in each
% dimension to obtain the maximum velocity in that dimension.

if (nargin == 0) & (nargout ==0)
    
    
 fprintf('                     Display: [ off | on | {iter} ]\n');   
 fprintf('                     WaitBar: [ {off} | on ]\n'); 
 fprintf('                        Plot: [ {off} | output | history | designspace ]\n'); 
 fprintf('               AlgorithmType: [ synchronous | {asynchronous} ]\n');  
 fprintf('                NumParticles: [ positive scalar between 5 and 500] {20}\n');  
 fprintf('                 MaxFunEvals: [ positive scalar between NumParticles and 1e12] {10000}\n'); 
 fprintf('                      TolFun: [ positive scalar ] {1e-4}\n'); 
 fprintf('                        TolX: [ positive scalar ] {1e-4}\n'); 
 fprintf('                      TolCen: [ positive scalar ] {1e-4}\n');
 fprintf('                     MaxTime: [ positive scalar | {7200} ]\n'); 
 fprintf('        InitializationMethod: [ {rand} | lhs ]\n'); 
 fprintf('                BoundsMethod: [ rand | {bounce} | drift ]\n'); 
 fprintf('          ModificationMethod: [ constriction | {dynamic} ]\n');
 fprintf('              InitialPenalty: [ positive scalar between 0 and 1000] {100}\n');
 fprintf('                FinalPenalty: [ positive scalar between 0 and 1000] {1000}\n');
 fprintf('            PenaltyChangeEnd: [ positive scalar between 0 and MaxFunEvals] {MaxFunEvals}\n');
 fprintf('              InitialInertia: [ positive scalar between 0 and 2 ] {1}\n');
 fprintf('             CognitiveWeight: [ positive scalar between 0 and 5 ] {2}\n');
 fprintf('                SocialWeight: [ positive scalar between 0 and 5 ] {2}\n');
 fprintf('     IterationsNoImprovement: [ positive scalar between 1 and MaxFunEvals/NumParticles] {10}\n'); 
 fprintf('            NoImprovementTol: [ positive scalar between 0 and inf] {0.01}\n'); 
 fprintf('            LimitMaxVelocity: [ binary 0 or 1] {1}\n');
 fprintf('    InertiaReductionFraction: [ positive scalar ] {0.0500}\n'); 
 fprintf('   VelocityReductionFraction: [ positive scalar ] {0.0500}\n'); 
 fprintf('            VelocityFraction: [ positive scalar ] {0.5000}\n');
 fprintf('\n');
 return;
end


options = struct('Display', ['iter'],...
    'WaitBar',['off'],...
    'Plot',['off'],...
    'AlgorithmType',['asynchronous'],...
    'NumParticles',[20],...
    'MaxFunEvals',[10000],...
    'TolFun',[1e-4],...
    'TolX',[1e-4],...
    'TolCen',[1e-4],...
    'MaxTime',[7200],...
    'InitializationMethod',['rand'],...
    'BoundsMethod',['bounce'],...
    'ModificationMethod',['dynamic'],...
    'InitialPenalty',[100],...
    'FinalPenalty',[1000],...
    'PenaltyChangeEnd',[10000],...
    'InitialInertia',[1],...
    'CognitiveWeight',[2],...
    'SocialWeight',[2],...
    'IterationsNoImprovement',[10],...
    'NoImprovementTol',[0.01],...
    'LimitMaxVelocity',[1],... 
    'InertiaReductionFraction',[0.05],...
    'VelocityReductionFraction',[0.05],...
    'VelocityFraction',[0.5]);   



optNames = fieldnames(options);
[m,n]=size(optNames);
optnames=lower(optNames);
PenaltyChgEndSetByUser = 0;

i=1;
while i<=nargin,
    arg = lower(varargin{i});
    arg2 = lower(varargin{i+1});
    if ~ischar(arg)
        error(sprintf('psoset was expecting a string in argument %d, but did not find one',i));
    end
    index = -1;
    for j=1:m,
        if (strcmp(arg,optnames(j))),
                index = j;
                validityCheck(arg,arg2);
                options.(optNames{j}) = arg2;
        end
    end
    if index == 16,
        PenaltyChgEndSetByUser = 1;
    end
    if index < 0,
        error(sprintf('%s is not a valid field name',arg));
    end
    
    i = i +2;
end

if ~PenaltyChgEndSetByUser,
    options.PenaltyChangeEnd = options.MaxFunEvals;
end

if options.MaxFunEvals < options.NumParticles,
    error(sprintf('"MaxFunEvals" (%d) is less than "NumParticles" (%d)',options.MaxFunEvals,options.NumParticles));
end
if options.PenaltyChangeEnd > options.MaxFunEvals,
    error(sprintf('"PenaltyChangeEnd" (%d) is greater than "MaxFunEvals" (%d)',options.PenaltyChangeEnd,options.MaxFunEvals));
end
if options.IterationsNoImprovement < 1 || options.IterationsNoImprovement > (options.MaxFunEvals/options.NumParticles),
    error(sprintf('The value for "IterationsnoImprovement" is not an accepted value.\n*%d* was entered and it must be: [positive scalar between 1 and MaxFunEvals/NumParticles (%d)]',options.IterationsNoImprovement,options.MaxFunEvals/options.NumParticles));
end
    
function bool = validityCheck(fieldname, value)
 bool = 0;
switch fieldname
    case {'display'}
        if (ischar(value))&&(strcmp(value,'off')||strcmp(value,'on')||strcmp(value,'iter')),
            bool = 1;
            return;
        else
            error(sprintf('The value for "Display" is not an accepted value.\n*%s* was entered and it must be either: off|on|iter',value));
        end
        
    case {'waitbar'}
        if (ischar(value))&&(strcmp(value,'off')||strcmp(value,'on')),
            bool = 1;
            return;
        else
            error(sprintf('The value for "WaitBar" is not an accepted value.\n*%s* was entered and it must be either: off|on',value));
        end
    case {'plot'}
        if (ischar(value))&&(strcmp(value,'off')||strcmp(value,'output')||strcmp(value,'history')||strcmp(value,'designspace')),
            bool = 1;
            return;
        else
            error(sprintf('The value for "Plot" is not an accepted value.\n*%s* was entered and it must be either: off|output|history|designspace',value));
        end
    case {'maxfunevals'}
        if (isnumeric(value))&& value > 0 && value <= 1e12,
            bool = 1;
            return;
        else
            error(sprintf('The value for "MaxFunEvals" is not an accepted value.\n*%d* was entered and it must be: [positive scalar between NumParticles and 1e12]',value));
        end
    case {'tolfun'}
        if (isnumeric(value))&& (value >= 0),
            bool = 1;
            return;
        else
             error(sprintf('The value for "TolFun" is not an accepted value.\n*%d* was entered and it must be: [positive scalar]',value));
        end
    case {'tolx'}
        if (isnumeric(value))&&(value >= 0),
            bool = 1;
            return;
        else
             error(sprintf('The value for "TolX" is not an accepted value.\n*%d* was entered and it must be: [positive scalar]',value));
        end
    case {'noimprovementtol'}
        if (isnumeric(value))&& value > 0,
            bool = 1;
            return;
        else
            error(sprintf('The value for "NoImprovementTol" is not an accepted value.\n*%d* was entered and it must be: [positive scalar]',value));
        end
    case {'maxtime'}
        if (isnumeric(value))&& value > 0,
            bool = 1;
            return;
        else
            error(sprintf('The value for "MaxTime" is not an accepted value.\n*%d* was entered and it must be: [positive scalar]',value));
        end
    case {'initializationmethod'}
        if (ischar(value))&&(strcmp(value,'rand')||strcmp(value,'lhs')),
            bool = 1;
            return;
        else
            error(sprintf('The value for "InitializationMethod" is not an accepted value.\n*%s* was entered and it must be: rand|lhs',value));
        end
    case {'boundsmethod'}
        if (ischar(value))&&(strcmp(value,'rand')||strcmp(value,'bounce')||strcmp(value,'drift')),
            bool = 1;
            return;
        else
            error(sprintf('The value for "BoundsMethod" is not an accepted value.\n*%s* was entered and it must be: rand|bounce|drift',value));
        end
    case {'modificationmethod'}
        if (ischar(value))&&(strcmp(value,'constriction')||strcmp(value,'dynamic')),
            bool = 1;
            return;
        else
            error(sprintf('The value for "ModificationMethod" is not an accepted value.\n*%s* was entered and it must be: constriction|dynamic',value));
        end
    case {'numparticles'}
        if (isnumeric(value))&& value >= 5 && value <=500,
            bool = 1;
            return;
        else
            error(sprintf('The value for "NumParticles" is not an accepted value.\n*%d* was entered and it must be: [positive integer between 5 and 500]',value));
        end
    case {'algorithmtype'}
        if (ischar(value))&&(strcmp(value,'synchronous')||strcmp(value,'asynchronous')),
            bool = 1;
            return;
        else
            error(sprintf('The value for "AlgorithmType" is not an accepted value.\n*%s* was entered and it must be: synchronous|asynchronous',value));
        end
    case {'initialpenalty'}
        if (isnumeric(value))&& ((value >= 0) & (value <= 1000)),
            bool = 1;
            return;
        else
            error(sprintf('The value for "InitialPenalty" is not an accepted value.\n*%d* was entered and it must be: [positive scalar between 0 and 1000]',value));
        end
    case {'finalpenalty'}
        if (isnumeric(value))&& ((value >= 0) & (value <= 1000)),
            bool = 1;
            return;
        else
            error(sprintf('The value for "FinalPenalty" is not an accepted value.\n*%d* was entered and it must be: [positive scalar between 0 and 1000]',value));
        end
    case {'penaltychangeend'}
        if (isnumeric(value))&& (value >= 0),
            bool = 1;
            return;
        else
            error(sprintf('The value for "PenaltyChangeEnd" is not an accepted value.\n*%d* was entered and it must be: [positive scalar between 0 and MaxFunEvals]',value));
        end  
    case {'initialinertia'}
        if (isnumeric(value))&& (value >= 0 & value <=2),
            bool = 1;
            return;
        else
            error(sprintf('The value for "InitialInertia" is not an accepted value.\n*%d* was entered and it must be: [positive scalar between 0 and 2]',value));
        end
    case {'cognitiveweight'}
        if (isnumeric(value))&& ((value >= 0) & (value <= 5))
            bool = 1;
            return;
        else
            error(sprintf('The value for "CognitiveWeight" is not an accepted value.\n*%d* was entered and it must be: [positive scalar between 0 and 5]',value));
        end
    case {'socialweight'}
        if (isnumeric(value))&& ((value >= 0) & (value <= 5))
            bool = 1;
            return;
        else
            error(sprintf('The value for "SocialWeight" is not an accepted value.\n*%d* was entered and it must be: [positive scalar between 0 and 5]',value));
        end
    case {'iterationsnoimprovement'}
        if (isnumeric(value))&& ((value >0))
            bool = 1;
            return;
        else
            error(sprintf('The value for "IterationsNoImprovement" is not an accepted value.\n*%d* was entered and it must be: [positive scalar between 1 and MaxFunEvals/NumParticles]',value));
        end
    case {'inertiareductionfraction'}
        if (isnumeric(value))&&(value >=0 & value <=1)
            bool = 1;
            return;
        else
            error(sprintf('The value for "InertiaReductionFraction" is not an accepted value.\n*%d* was entered and it must be: [positive scalar between 0 and 1]',value));
        end
    case {'velocityreductionfraction'}
        if (isnumeric(value))&&(value >=0 & value <=1)
            bool = 1;
            return;
        else
            error(sprintf('The value for "VeloctiyReductionFraction" is not an accepted value.\n*%d* was entered and it must be: [positive scalar between 0 and 1]',value));
        end
    case {'velocityfraction'}
        if (isnumeric(value))&&(value >=0 & value <=1)
            bool = 1;
            return;
        else
            error(sprintf('The value for "VeloctiyFraction" is not an accepted value.\n*%d* was entered and it must be: [positive scalar between 0 and 1]',value));
        end
    case {'noimprovementtol'}
        if (isnumeric(value))&&(value >=0)
            bool = 1;
            return;
        else
            error(sprintf('The value for "NoImprovementTol" is not an accepted value.\n*%d* was entered and it must be: [positive scalar between 0 and inf]',value));
        end   
    case {'limitmaxvelocity'}
        if (isnumeric(value))&&(value >=0 & value <2)
            bool = 1;
            return;
        else
            error(sprintf('The value for "LimitMaxVelocity" is not an accepted value.\n*%d* was entered and it must be: [boolean 0 or 1]',value));
        end          
    otherwise
        disp('this case is not yet covered in code');
        bool =1;
        return
end