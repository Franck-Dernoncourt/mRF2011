function [  ] = draw_data_stats( nb_objectives, file_name, nb_ind_per_gen )
%DRAW_DATA_STATS Draws all the data stats.
%   draw_data_stats() draws all the data stats.
%                    Default nb_objectives is 2.
%                    Default data file is 'all_phenotypes_data_stats.dat'.
%                    Default nb_ind_per_gen is 500.
%
%   draw_data_stats(nb_objectives) draws all the data stats, knowing thaat
%                    the number of objectives is nb_objectives. 
%                    Default data file is 'all_phenotypes_data_stats.dat'.
%
%   draw_data_stats(nb_objectives, file_name) draws the #objective_number 
%                    objectives' best score during evolution. Data file is
%                    file_name.
%
%   draw_data_stats(nb_objectives, file_name, nb_ind_per_gen) draws the 
%                    objective_number objectives' best score of the
%                    nb_ind_per_gen first individuals during evolution. 
%                    nb_ind_per_gen can be useful since drawing a number of
%                    points can cause slowdown. Data file is file_name.
%
%   Expected format of the data file:
%   Generation# Individual# Score#1 Score#2 ... Score#m Stat#1 Stat#2 ... Stat#n
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
%   0,0,-0.616708,2.93416,0.230769,0,30,2,0.333333
%   0,1,-0.629622,4.00559,0.142857,0,41,2,0.208333
%   0,2,-0.631582,6.17913,0.363636,0,62,2,0.315789
%   0,3,-0.631582,6.17292,0.382353,0,62,2,0.35
%   0,4,-0.631582,5.77476,0.363636,0,58,2,0.333333
%   0,5,-0.644906,3.01381,0.214286,0,31,2,0.375
%   0,6,-0.655063,2.93416,0.230769,0,30,2,0.333333
%   0,7,-0.664914,2.04261,0.521739,0,21,2,0.428571
%   0,8,-0.664914,1.94261,0.521739,0,20,2,0.428571
%   0,9,-0.666667,7.0575,0.25,0,71,2,0.25
%   0,10,-0.666667,6.95601,0.25,0,70,2,0.263158
%   0,11,-0.666667,6.95601,0.25,0,70,2,0.263158
%   [...]
%

%% Argument parsing
if nargin < 1
    nb_objectives = 2;
end

if nargin < 2
    file_name = 'all_phenotypes_data_stats.dat';
end

if nargin < 3
    nb_ind_per_gen = 500;
end

% Load the file output by sferes
data_stats = load(file_name);


%% Filtering
scores_to_plot = data_stats;
scores_to_plot = scores_to_plot(scores_to_plot(:, 2)<nb_ind_per_gen, :);

%% Draw the matrix
figure
[~,AX,~,~] = plotmatrix(scores_to_plot);
title('Evolution of the data stats', 'FontWeight','bold')
set(get(AX(1, 1),'ylabel'),'string','Gen#')
set(get(AX(2, 1),'ylabel'),'string','Ind#')
set(get(AX(end, 1),'xlabel'),'string','Gen#')
set(get(AX(end, 2),'xlabel'),'string','Ind#')
for i=1:nb_objectives
    set(get(AX(2 + i, 1),'ylabel'),'string', strcat('Obj#',int2str(i)))
    set(get(AX(end, 2 + i),'xlabel'),'string', strcat('Obj#',int2str(i)))
end
stat_number = 1;
for i=(2+nb_objectives + 1):size(AX, 1)
    set(get(AX(i, 1),'ylabel'),'string', strcat('Stat#',int2str(stat_number)))
    set(get(AX(end, i),'xlabel'),'string', strcat('Stat#',int2str(stat_number)))
    stat_number = stat_number + 1;
end

end