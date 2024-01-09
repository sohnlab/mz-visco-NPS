function EventSelectionGUI(t, signal, startTimes, saveDir, saveName, ...
                           plotBuffer, endBuffer)
% This function runs the selection GUI for manual selection of cell transit
% events identified in [signal].
%
%          [t] Signal time vector.
%     [signal] Signal vector, should be 'smoothV1'.
% [startTimes] A [nx1] vector contianing the sample indices at which a 
%              cell transit event was identified.
%    [saveDir] File directory to save the output [selectedEvents] .mat
%              file. Should be same directory as original raw data file.
%   [saveName] File name of saved output. Should reference the name of 
%              orginal raw data file.
% [plotBuffer] How many extra samples to plot before the identified start
%              time, and after the anticipated end time.
%  [endBuffer] The number of samples anticipated to be the end of the cell
%              transit event after the identified start time.

%% Initialize Variables
selectedEvents = [];
numEvents = length(startTimes);

%% GUI Component Set-up for Manual Filter of Auto Event Selection

% Create UIFigure and hide until all components are created
fig = uifigure('Visible', 'off');
fig.Position = [100 100 424 488];
fig.Name = 'ViscoNPS Event Selection Tool';

% Create UIAxes
figAxes = uiaxes(fig);
figAxes.Position = [10 58 404 382];
figAxes.Title.String = sprintf('Cell Transit Event 1 of %d', numEvents);
figAxes.XLabel.String = 'Time [s]';
figAxes.YLabel.String = 'Voltage [V]';

% Create UILabels
figLabel1 = uilabel(fig);
figLabel1.FontSize = 16;
figLabel1.FontWeight = 'bold';
figLabel1.Position = [90 467 405 22];
figLabel1.Text = '6-Zone ViscoNPS - Event Selection';

figLabel2 = uilabel(fig);
figLabel2.FontAngle = 'italic';
figLabel2.FontColor = [0.4667 0.6745 0.1882];
figLabel2.Position = [150 446 250 22];
figLabel2.Text = sprintf('Total: %d events identified.', numEvents);

% Create ACCEPT Button
figACCEPTButton = uibutton(fig, 'push');
figACCEPTButton.BackgroundColor = [0.4667 0.6745 0.1882];
figACCEPTButton.FontWeight = 'bold';
figACCEPTButton.FontColor = [1 1 1];
figACCEPTButton.Position = [200 21 100 23];
figACCEPTButton.Text = 'ACCEPT';

% ACCEPT Button Response
figACCEPTButton.ButtonPushedFcn = @(btn,event) acceptButtonPushed();

% Create REJECT Button
figREJECTButton = uibutton(fig, 'push');
figREJECTButton.BackgroundColor = [1 0 0];
figREJECTButton.FontWeight = 'bold';
figREJECTButton.FontColor = [1 1 1];
figREJECTButton.Position = [310 21 100 23];
figREJECTButton.Text = 'REJECT';

% ACCEPT Button Response
figREJECTButton.ButtonPushedFcn = @(btn,event) rejectButtonPushed();

% Add function for key press
fig.KeyPressFcn = @keyboardButtonPushed;

%% Initialize First Plot

plotIndex = 1;
startSample = startTimes(plotIndex)-plotBuffer;
if startSample < 1
    startSample = 1;
end
endSample = startTimes(plotIndex)+endBuffer+plotBuffer;

p = plot( figAxes,t(startSample:endSample),signal(startSample:endSample) );
figAxes.XLim = [t(startSample), t(endSample)];
drawnow;

%% Show GUI After All Components Initialized
fig.Visible = 'on';

%% Callback Functions and Helper Functions

    function keyboardButtonPushed(src, event)
        if strcmp(event.Key, 'return')
            acceptButtonPushed()
        end
        if strcmp(event.Key, 'space')
            rejectButtonPushed()
        end
    end

    function acceptButtonPushed()
        % GUI Button Callback for "Accept Button"
        % This function is executed when the user presses "Accept"
        
        if plotIndex < numEvents
            selectedEvents = [selectedEvents; startTimes(plotIndex)];
            newPlot()
        elseif plotIndex == numEvents
            selectedEvents = [selectedEvents; startTimes(plotIndex)];
            saveData()
        else
            msgbox('No more cell transit events to select!','Complete!');
        end
        
    end

    function rejectButtonPushed()
        % GUI Button Callback for "Reject Button"
        % This function is executed when the user presses "Reject"
        
        if plotIndex == numEvents
            saveData()
        elseif plotIndex > numEvents
            msgbox('No more cell transit events to select!','Complete!');
        else
            newPlot() % plot next cell transit event
        end
            
    end

    function newPlot()
        % This executes after one of the button callback functions are
        % executed to plot the next cell transit event
        
        %cla(figAxes); % clear the figure axes
        plotIndex = plotIndex + 1; % iterate plot index
        
        % create new plot only if there are more events
        if plotIndex <= numEvents
            startSample = startTimes(plotIndex)-plotBuffer;
            endSample = startTimes(plotIndex)+endBuffer+plotBuffer;
            if endSample > length(signal)
                endSample = length(signal);
            end
            
            % create new plot on the same figAxes
            %plot(figAxes,t(startSample:endSample),...
            %             signal(startSample:endSample));
            p.XData = t(startSample:endSample);
            p.YData = signal(startSample:endSample);
            figAxes.XLim = [t(startSample), t(endSample)];
            figAxes.Title.String = sprintf('Cell Transit Event %d of %d', ...
                            plotIndex,numEvents);
            drawnow;
        else
           msgbox('No more cell transit events to select!','Complete!');
        end
        
    end

    function saveData()
        plotIndex = plotIndex + 1;
        fullDir = strcat(saveDir, saveName, '-EventSelection.mat');
        save(fullDir, 'selectedEvents');
        message = sprintf(['Data Saved!\nTotal: %d events out of %d' ...
                           ' selected.'],...
                          length(selectedEvents), numEvents);
        waitfor(msgbox(message,'Complete!'));
        close(fig);
    end


end


