function ZoneQCGUI(EventSignals, EventSignalsParsed, numZones, saveDir, saveName)
% Andre Lai | June 25, 2022
%
% This function runs the Zone QC GUI for manual selection of cell transit
% events with parsed contraction zone signal in [EventSignals].
%
%       [EventSignals] Matrix containing selected events.
% [EventSignalsParsed] Matrix containing the identified start and end times
%                      of each zone signal for each event in [EventSignals].
%           [numZones] Number of zones measured (including sizing pore).
%            [saveDir] File directory to save the output [selectedEvents]
%                      .mat file. Should be same directory as original raw 
%                      data file.
%           [saveName] File name of saved output. Should reference the name 
%                      of orginal raw data file.

%% Initialize Variables
selectedEvents2 = false(size(EventSignals,3),1);
numEvents = size(EventSignals,3);

%% GUI Component Set-up for Manual Filter of Auto Event Selection

% Create UIFigure and hide until all components are created
fig = uifigure('Visible', 'off');
left = 400;
bottom = 150;
width = left + 100;
height = bottom + numZones*140;
fig.Position = [left bottom width height];
fig.Name = 'ViscoNPS Zone QC Tool';

% initialize and add [numZones] number of plots axes
t = tiledlayout(fig, numZones, 1, 'tilespacing', 'compact');
axesCell = cell(numZones,1);
for i = 1:numZones
    axesCell{i} = nexttile(t);
    if i ~= numZones
        axesCell{i}.XTickLabel = "";
    end
end
axesCell{1}.Title.String = sprintf('Cell Transit Event 1 of %d', numEvents);
axesCell{1}.Title.FontSize = 14;
t.YLabel.String = 'Resistance [\Omega]';
t.YLabel.FontSize = 14; 
t.XLabel.String = 'Time [s]';
t.XLabel.FontSize = 14; 

% Create UILabels
figLabel1 = uilabel(fig);
figLabel1.FontSize = 16;
figLabel1.FontWeight = 'bold';
figLabel1.Position = [1 0.96*height fig.Position(3) 22];
figLabel1.Text = sprintf('%d-Zone ViscoNPS - Zone Identification QC',numZones);
figLabel1.HorizontalAlignment = 'center';

% Create ACCEPT Button
figACCEPTButton = uibutton(fig, 'push');
figACCEPTButton.BackgroundColor = [0.4667 0.6745 0.1882];
figACCEPTButton.FontWeight = 'bold';
figACCEPTButton.FontColor = [1 1 1];
figACCEPTButton.Position = [250 10 100 30];
figACCEPTButton.Text = 'ACCEPT';

% ACCEPT Button Response
figACCEPTButton.ButtonPushedFcn = @(btn,event) acceptButtonPushed();

% Create REJECT Button
figREJECTButton = uibutton(fig, 'push');
figREJECTButton.BackgroundColor = [1 0 0];
figREJECTButton.FontWeight = 'bold';
figREJECTButton.FontColor = [1 1 1];
figREJECTButton.Position = [360 10 100 30];
figREJECTButton.Text = 'REJECT';

% ACCEPT Button Response
figREJECTButton.ButtonPushedFcn = @(btn,event) rejectButtonPushed();

% Add function for key press
fig.KeyPressFcn = @keyboardButtonPushed;

%% Initialize Plot of First Event

eventNum = 1;

% common parameters
timeSignal = EventSignals(:,numZones*2+1,eventNum);
startTime = timeSignal(1); 
endTime = timeSignal(end);
yStart = min(min(EventSignals(:,1:numZones,eventNum)));
yEnd = max(max(EventSignals(:,1:numZones,eventNum)));

colors = [[60/255,60/255,60/255];[0, 0.4470, 0.7410]; ...
          [0.8500, 0.3250, 0.0980] ;[0.9290, 0.6940, 0.1250]; ...
          [0.4940, 0.1840, 0.5560] ;[0.4660, 0.6740, 0.1880]];

% cell array to hold plot components
% there are 3 plots per axes: 
% 1-parsed zone signal, 2-full event signal, and 3-baseline
plotCell = cell(numZones,3); 
      
for z = 1:numZones
    
    hold( axesCell{z}, 'on' )
    % parsed zone signal
    plotCell{z,1} = plot(axesCell{z}, timeSignal(EventSignalsParsed(1,z,eventNum): ...
                                                 EventSignalsParsed(2,z,eventNum)), ...
                                      EventSignals(EventSignalsParsed(1,z,eventNum): ...
                                                   EventSignalsParsed(2,z,eventNum),...
                                                   z,eventNum), ...
                        'LineWidth', 3, 'Color',colors(z,:), ...
                        'ButtonDownFcn',@mouseClick);
    % full event signal                
    plotCell{z,2} = plot(axesCell{z},timeSignal, EventSignals(:,z,eventNum), ...
                         'Color',[160/255,160/255,160/255],'LineWidth', 1+z/100, ...
                         'ButtonDownFcn',@mouseClick);
    % linear baseline 
    plotCell{z,3} = plot(axesCell{z},timeSignal, EventSignals(:,z,eventNum) .* 0 + ...
                                                 median(EventSignals(:,z,eventNum)), ...
                         'LineWidth', 0.5, 'color', [0.9247, 0.3297, 0.2244]);
    % plot parameters
    axesCell{z}.XLim = [startTime, endTime];
    if z ~= 1
        axesCell{z}.YLim = [yStart, yEnd];
    end
    hold( axesCell{z}, 'off' )
    
end

drawnow;

%% Show GUI After All Components Initialized
fig.Visible = 'on';

%% Callback Functions and Helper Functions

    function keyboardButtonPushed(~, event)
        if strcmp(event.Key, 'return')
            acceptButtonPushed()
        elseif strcmp(event.Key, 'space')
            rejectButtonPushed()
        elseif strcmp(event.Key, 'leftarrow')
            if eventNum > 1
                eventNum = eventNum - 2;
                newPlot()
            else
                uiwait(warndlg(sprintf('Invalid Keystroke!\nThere is no previous plot.'),...
                           'Error', 'replace'))
                figure(fig) % call back figure window to the front
            end
        else
            uiwait(warndlg(sprintf(['Invalid Keystroke!\n[left click] on plot to update ',...
                                    'start and end point parsing.\nPress [enter] for Accept.',...
                                    '\nPress [space] for Reject.\nPress [left arrow] ',...
                                    'to go back to the previous plot.']),...
                           'Error', 'replace'))
            figure(fig) % call back figure window to the front
        end
    end

    function mouseClick(src, event)
        
        % get coordinates of click 
        xPoint = event.IntersectionPoint(1);
        
        % the zone clicked is coded in the 100th decimal of the LineWidth property
        zone = cast((src.LineWidth-1)*100, "uint8");
        
        % find the index of the clicked xPoint using the timeSignal vector
        tol = 0.0001;
        index_of_xPoint = find(abs(timeSignal-xPoint) < tol);
        if length(index_of_xPoint) > 1
            index_of_xPoint = index_of_xPoint(end);
        end
        
        % update the signal parsing
        if xPoint > timeSignal(EventSignalsParsed(1,zone,eventNum)) && ...
           xPoint > timeSignal(EventSignalsParsed(2,zone,eventNum))
            EventSignalsParsed(2,zone,eventNum) = index_of_xPoint;
            fprintf("Event %d, Zone %d, end-point updated.\n", eventNum, zone)
        elseif xPoint < timeSignal(EventSignalsParsed(1,zone,eventNum)) && ...
               xPoint < timeSignal(EventSignalsParsed(2,zone,eventNum))
            EventSignalsParsed(1,zone,eventNum) = index_of_xPoint;
            fprintf("Event %d, Zone %d, start-point updated.\n", eventNum, zone)
        elseif abs(xPoint - timeSignal(EventSignalsParsed(1,zone,eventNum))) < ...
               abs(xPoint - timeSignal(EventSignalsParsed(2,zone,eventNum)))
            EventSignalsParsed(1,zone,eventNum) = index_of_xPoint;
            fprintf("Event %d, Zone %d, start-point updated.\n", eventNum, zone)
        else
            EventSignalsParsed(2,zone,eventNum) = index_of_xPoint;
            fprintf("Event %d, Zone %d, end-point updated.\n", eventNum, zone)
        end
        
        % update the plot
        eventNum = eventNum - 1;
        newPlot()
        
    end

    function acceptButtonPushed()
        % GUI Button Callback for "Accept Button"
        % This function is executed when the user presses "Accept"
        
        if eventNum < numEvents
            selectedEvents2(eventNum) = 1;
            newPlot()
        elseif eventNum == numEvents
            selectedEvents2(eventNum) = 1;
            saveData()
        else
            msgbox('No more cell transit events to select!','Complete!');
        end
        
    end

    function rejectButtonPushed()
        % GUI Button Callback for "Reject Button"
        % This function is executed when the user presses "Reject"
        
        if eventNum == numEvents
            saveData()
        elseif eventNum > numEvents
            msgbox('No more cell transit events to select!','Complete!');
        else
            newPlot() % plot next cell transit event
        end
            
    end

    function newPlot()
        % This executes after one of the button callback functions are
        % executed to plot the next cell transit event
        
        %cla(figAxes); % clear the figure axes
        eventNum = eventNum + 1; % iterate event plot index
        
        % create new plot only if there are more events
        if eventNum <= numEvents
            
            % updated parameters
            timeSignal = EventSignals(:,numZones*2+1,eventNum);
            startTime = timeSignal(1); 
            endTime = timeSignal(end);
            yStart = min(min(EventSignals(:,1:numZones,eventNum)));
            yEnd = max(max(EventSignals(:,1:numZones,eventNum)));       
            
            % update plot data on the same UIAxes
            for zn = 1:numZones
                
                plotCell{zn,1}.XData = timeSignal(EventSignalsParsed(1,zn,eventNum): ...
                                                  EventSignalsParsed(2,zn,eventNum));
                plotCell{zn,1}.YData = EventSignals(EventSignalsParsed(1,zn,eventNum): ...
                                                    EventSignalsParsed(2,zn,eventNum),zn,eventNum);
                plotCell{zn,2}.XData = timeSignal;
                plotCell{zn,2}.YData = EventSignals(:,zn,eventNum);
                plotCell{zn,3}.XData = timeSignal;
                plotCell{zn,3}.YData = EventSignals(:,zn,eventNum) .* 0 + ...
                                      median(EventSignals(:,1,eventNum));
                axesCell{zn}.XLim = [startTime, endTime];
                if zn ~= 1
                    axesCell{zn}.YLim = [yStart, yEnd];
                end
                
            end
            
            axesCell{1}.Title.String = sprintf('Cell Transit Event %d of %d', ...
                                                eventNum,numEvents);
            drawnow;
            
        else
           msgbox('No more cell transit events!','Complete!');
        end
        
    end

    function saveData()
        eventNum = eventNum + 1;
        fullDir = strcat(saveDir, saveName, '-EventSignalsParsedManual.mat');
        save(fullDir, 'EventSignalsParsed');
        fullDir = strcat(saveDir, saveName, '-EventSelection2.mat');
        save(fullDir, 'selectedEvents2');
        message = sprintf(['Data Saved!\nTotal: %d events out of %d' ...
                           ' accepted.'],...
                          nnz(selectedEvents2), numEvents);
        waitfor(msgbox(message,'Complete!'));
        close(fig);
    end

end


