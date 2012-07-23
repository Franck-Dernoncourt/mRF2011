function [  ] = draw_objective( objective_number, file_name )
%DRAW_OBJECTIVE Draws an objective's best score during evolution.
%   draw_objective() draws the 1st objective's best score during evolution. 
%                    default data file is 'pareto.dat'
%
%   draw_objective(objective_number) draws the #objective_number objective's 
%                    best score during evolution. Default data file is 'pareto.dat'
%
%   draw_objective(objective_number, file_name) draws the #objective_number 
%                    objective's best score during evolution. Data file is
%                    file_name.
%
%   Expected format of the data file:
%   Generation# Individual# Score#1 Score#2 ... Score#n
% 
%   Values in columns must be real number (integer, float, etc)
%   Columns must be separated with a file delimiter (the character
%   between elements in each row) recognized by MATLAB's load(), i.e. it 
%   can be a blank, comma, semicolon, or tab character. The file can 
%   contain MATLAB comments (lines that begin with a percent sign, %).
%   See http://www.mathworks.com/help/techdoc/ref/load.html for more 
%   information on load().
% 
%   Example of data file's content:
%   0 0 -0.458084 2 
%   0 1 -0.470283 4 
%   0 2 -0.516762 5 
%   0 3 -0.52249 6 
%   1 0 -0.458084 2 
%   1 1 -0.46491 4 
%   1 2 -0.487641 5 
%   1 3 -0.519193 6 
%   2 0 -0.451223 5 
%   2 1 -0.512767 6 
%   3 0 -0.435107 4 
%   3 1 -0.44803 6 
%   4 0 -0.435107 4 
%   4 1 -0.445121 5 
%   4 2 -0.44803 6 
%   5 0 -0.435107 4 
%   5 1 -0.445121 5 
%   [...]
%
%   Note: Documents describing MATLAB coding standards or a good
%         programming guide can be found here:
%         http://www.mathworks.com/support/solutions/en/data/1-1AN6R/?solution=1-1AN6R

%% Argument parsing
if nargin < 1
    objective_number = 1;
end

if nargin < 2
    file_name = 'pareto.dat';
end

% Load the file output by sferes
scores = load(file_name);

% Check argument 'objective_number'
if ( not(isequal(fix(objective_number),objective_number))... % Integer
    || (objective_number > size(scores, 2) - 2)... % Not too big
    || (objective_number < 1)) % Strictly positive
    error('Sorry, the objective number is not valid')
end


%% Filter: only keep the best individual's score for each generation.
% (the best individual depends on objective_number)

% Since the 1st objective is numbered generation 0
if objective_number == 1
    scores_to_plot = scores(scores(:, 2)==0, [1 3]);
% Otherwise we'll have to loop to find which generation has the best score
else 
    scores_to_plot = []; 
    for generation_number = 0:max(scores(:, 1))
        generation_scores = scores(scores(:, 1) == generation_number, :);
        [generation_number max(generation_scores(:, 2+objective_number))];
        scores_to_plot = vertcat(scores_to_plot,[generation_number max(generation_scores(:, 2+objective_number))]);
    end
end

%% Draw the matrix
figure
[~,AX,~,~] = plotmatrix(scores_to_plot);
title(strcat('Evolution of the objective #', int2str(objective_number)), 'FontWeight','bold')
set(get(AX(1),'ylabel'),'string','Generation#')
set(get(AX(2),'xlabel'),'string','Generation#')
set(get(AX(2),'ylabel'),'string',strcat('Objective #', ...
    int2str(objective_number), ' score'))
set(get(AX(4),'xlabel'),'string',strcat('Objective #', ...
    int2str(objective_number), ' score'))

end