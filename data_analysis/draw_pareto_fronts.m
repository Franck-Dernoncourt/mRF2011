function [  ] = draw_pareto_fronts( objective_numbers, file_name, video_name )
%DRAW_PARETO_FRONTS Draws pareto fronts (2D or 3D) of all generations. It 
%                   can also record a rotation of the plot.
% 
%   draw_pareto_fronts() draws the 1st objective's best score during evolution. 
%                       default data file is 'pareto.dat'
%
%   draw_pareto_fronts(objective_numbers) draws the objectives
%                       objective_numbers's (matrix expected) best score 
%                       during evolution. Default data file is 'pareto.dat'
%                       objective_numbers's size must be 1, 2 or 3.
%                       if objective_numbers's size is 1, then we call 
%                       draw_objective().
%
%   draw_pareto_fronts(objective_numbers, file_name) draws the objectives
%                       objective_numbers's best score during evolution.
%                       Data file is file_name.
% 
%   draw_pareto_fronts(objective_numbers, file_name, record_video) draws 
%                       the objectives objective_numbers's best score 
%                       during evolution. Data file is file_name.
%                       It will record a rotation of the plot in video_name
%                       (AVI format).
%
%   Examples: 
%      * draw_pareto_fronts()
%      * draw_pareto_fronts([1:3])
%      * draw_pareto_fronts([1:3], 'pareto5.dat')
%      * draw_pareto_fronts([1], 'pareto.dat')
%      * draw_pareto_fronts([1:3], 'pareto.dat', 'pareto_video.avi')
%
%   Dependencies: draw_pareto_fronts.m depends on draw_objective.m
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

%% Argument parsing
if nargin < 1
    objective_numbers = [1 2];
end

if nargin < 2
    file_name = 'pareto.dat';
end

record_video = (nargin > 2);

% Load the file output by sferes
scores = load(file_name);


%% Draw pareto fronts for each generation
scrsz = get(0,'ScreenSize');
h = figure('Position',[30 (scrsz(4)/2 - 100) scrsz(3)/2 scrsz(4)/2]);
% 1D
if (length(objective_numbers) == 1 || size(scores, 2) == 3)
    draw_objective(objective_numbers(1), file_name);
    return

% 2D
elseif (length(objective_numbers) == 2 || size(scores, 2) == 4)
    scores_to_plot = scores(:, [3 4 1]);
    scatter(scores_to_plot(:, 1), scores_to_plot(:, 2),5, ...
        scores_to_plot(:, 3));

% 3D
elseif (size(scores, 2) == 5)
    scores_to_plot = scores(:, [3 4 5 1]);
    scores_plot = scatter3(scores_to_plot(:, 1), scores_to_plot(:, 2), ...
        scores_to_plot(:, 3),5, scores_to_plot(:, 4));
    zlabel('Objective 3');


end
cbar_axes = colorbar;
set(get(cbar_axes, 'title'), 'string','Generation#');
xlabel('Objective 1');
ylabel('Objective 2');

%% Record video
if record_video
    % Wait bar
    h_wait = waitbar(0,'Please wait...', ...
        'Name','Please wait, the video is being recorded...', ...
        'Position', [800 500 275 50]);

    % Video parameters
    fps = 100; 
    rotation_step_x = 0.3;
    sec = 360 / rotation_step_x / fps ;
    assert(rotation_step_x * fps * sec == 360,...
        'ERROR: please make sure that rotation_step_x * fps * sec == 360, otherwise the video will not show a full rotation')
    vidObj = VideoWriter(video_name);
    vidObj.Quality = 100;
    vidObj.FrameRate = fps;

    % Video recording
    open(vidObj);
    for i=1:fps*sec
      if (mod(i, 5) == 0)
          waitbar(i/(fps*sec),h_wait,horzcat(int2str(i),...
              ' frames out of ', int2str(fps*sec)))
      end
      camorbit(rotation_step_x,0);
      drawnow;
      writeVideo(vidObj,getframe(gcf));
    end
    close(vidObj);
    delete(h_wait)
end