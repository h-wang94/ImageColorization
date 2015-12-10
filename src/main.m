% default to run through everything

run_colorization = true;

if run_colorization == true
    disp('Running Image Colorization');
    run('initialize.m');
else
    disp('NOT running run_colorization. Edit main.m to toggle');
end

disp('Done with script!');
