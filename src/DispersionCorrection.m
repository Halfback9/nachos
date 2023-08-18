classdef DispersionCorrection < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)

        % Components
        UIFigure  matlab.ui.Figure
        UIAxes    matlab.ui.control.UIAxes
        FinishButton matlab.ui.control.Button
        Instructions matlab.ui.control.Label
        ColourRangeEditField       matlab.ui.control.NumericEditField
        ColourRangeEditFieldLabel  matlab.ui.control.Label
        
        % Data
        TASurface (:, :) {mustBeNumeric}
        pixels (1, :) {mustBeNumeric}
        times (:, 1) {mustBeNumeric}

        % Plots
        surfPlot
        scatterPlot
        linePlot
        dispersionCurveXVals {mustBeNumeric}
        dispersionCurveYVals {mustBeNumeric}
        dispersionFitCoefficients (1, :) {mustBeNumeric}
        dispersionFit  (1, :) {mustBeNumeric}
        dispersionZVals {mustBeNumeric}


    end

    % Callbacks that handle component events
    methods (Access = private)

        % Button down function: UIAxes
        function onSurfaceClicked(app, event)

            if event.Button == 1

                app.dispersionCurveXVals = [app.dispersionCurveXVals, event.IntersectionPoint(1)];
                app.dispersionCurveYVals = [app.dispersionCurveYVals, event.IntersectionPoint(2)];
                app.updateDispersionPoints();
            elseif event.Button == 3
                app.dispersionCurveXVals(end) = [];
                app.dispersionCurveYVals(end) = [];
                app.updateDispersionPoints();
            else
            end

        end

        function onFinishButtonPushed(app, event)
            
            app.close_window(false);
        end

        function onKeyPressed(app, event)
            key = event.Key;
            switch key
                case 'return'
                    app.close_window(false);
                case 'escape'
                    app.close_window(true);
            end
        end

        function close_window(app, cancelled)
            if cancelled
                app.dispersionFitCoefficients = zeros(1, 3);
                app.dispersionFit = zeros(size(app.pixels));
            end

            close(app.UIFigure);
        end

        function updateDispersionPoints(app)
            Zs = ones(size(app.dispersionCurveXVals))*max(max(app.TASurface))*1.1;
            

            set(app.scatterPlot, 'XData', app.dispersionCurveXVals, 'YData', app.dispersionCurveYVals, 'ZData', Zs);

            if length(Zs) == 2
                %app.dispersionFitCoefficients = [0, polyfit(app.dispersionCurveXVals, app.dispersionCurveYVals, 1)];
                app.dispersionFitCoefficients = polyfit(app.dispersionCurveXVals, app.dispersionCurveYVals, 1);
                %app.dispersionFit = app.dispersionFitCoefficients(1).*app.pixels.^2 + app.dispersionFitCoefficients(2).*app.pixels;
                app.dispersionFit = app.dispersionFitCoefficients(1).*app.pixels + app.dispersionFitCoefficients(2);
                set(app.linePlot, 'XData', app.pixels, 'YData', app.dispersionFit, 'ZData', app.dispersionZVals);
            elseif length(Zs) == 3
                %app.dispersionFitCoefficients = polyfit(app.dispersionCurveXVals, app.dispersionCurveYVals, 2);
                app.dispersionFitCoefficients = polyfit(app.dispersionCurveXVals, app.dispersionCurveYVals, 2);
                app.dispersionFit = app.dispersionFitCoefficients(1).*app.pixels.^2 + app.dispersionFitCoefficients(2).*app.pixels + app.dispersionFitCoefficients(3);
                set(app.linePlot, 'XData', app.pixels, 'YData', app.dispersionFit, 'ZData', app.dispersionZVals);
            elseif length(Zs) > 3
                %app.dispersionFitCoefficients = [0, polyfit(app.dispersionCurveXVals, app.dispersionCurveYVals, 3)];
                app.dispersionFitCoefficients = polyfit(app.dispersionCurveXVals, app.dispersionCurveYVals, 3);
                app.dispersionFit = app.dispersionFitCoefficients(1).*app.pixels.^3 + app.dispersionFitCoefficients(2).*app.pixels.^2 + app.dispersionFitCoefficients(3).*app.pixels + app.dispersionFitCoefficients(4);
                set(app.linePlot, 'XData', app.pixels, 'YData', app.dispersionFit, 'ZData', app.dispersionZVals);
            else
                app.dispersionFit = zeros(size(app.pixels));
            end

        end

        function changeClim(app, event)
            val = event.Value;
            try
                clim(app.UIAxes, [-1*abs(val), abs(val)]);
            catch
                caxis(app.UIAxes, [-1*abs(val), abs(val)]);
            end
        end
    end

    methods (Access = private)
        function createComponents(app)
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 640 480];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.KeyPressFcn = createCallbackFcn(app, @onKeyPressed, true);
    
            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            title(app.UIAxes, 'Dispersion Correction');
            xlabel(app.UIAxes, 'Pixels');
            ylabel(app.UIAxes, 'Time');
            app.UIAxes.ButtonDownFcn = createCallbackFcn(app, @onSurfaceClicked, true);
            app.UIAxes.Position = [50 50 540 380];

            % Create FinishButton
            app.FinishButton = uibutton(app.UIFigure, 'push');
            app.FinishButton.ButtonPushedFcn = createCallbackFcn(app, @onFinishButtonPushed, true);
            app.FinishButton.Position = [367 8 98 34];
            app.FinishButton.Text = 'Finish';

            app.Instructions = uilabel(app.UIFigure);
            app.Instructions.Position = [98 1 205 70];
            app.Instructions.Text = {'Left Click: Add point'; 'Right Click: Delete Point'; 'Enter: Finish'; 'Esc: Cancel'};

            % Create ColourRangeEditFieldLabel
            app.ColourRangeEditFieldLabel = uilabel(app.UIFigure);
            app.ColourRangeEditFieldLabel.HorizontalAlignment = 'right';
            app.ColourRangeEditFieldLabel.Position = [50 444 54 22];
            app.ColourRangeEditFieldLabel.Text = 'Colour Range';

            % Create ColourRangeEditField
            app.ColourRangeEditField = uieditfield(app.UIFigure, 'numeric');
            app.ColourRangeEditField.ValueChangedFcn = createCallbackFcn(app, @changeClim, true);
            app.ColourRangeEditField.Position = [144 444 100 22];
    
            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end

    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = DispersionCorrection(pixels, times, data)
            
            app.pixels = pixels(:)';
            app.times = times(:);

            if size(app.times*app.pixels) ~= size(data)
                data = data';
            end

            app.TASurface = data;

            app.dispersionZVals = ones(size(app.pixels))*max(max(app.TASurface))*1.05;

            % Create UIFigure and components
            createComponents(app)

            [app.dispersionCurveXVals, app.dispersionCurveYVals] = deal([]);

            app.dispersionFit = zeros(size(app.pixels));
            
            app.surfPlot = surf(app.UIAxes, app.pixels, app.times, app.TASurface, 'EdgeColor','none', 'HitTest','off');
            %app.surfPlot = pcolor(app.UIAxes, app.pixels, app.times, app.TASurface, 'HitTest','off');
            hold(app.UIAxes, 'on');
            app.scatterPlot = scatter(app.UIAxes, app.dispersionCurveXVals, app.dispersionCurveYVals, 20, 'r', 'filled');
            
            app.linePlot = plot(app.UIAxes, app.pixels, app.dispersionFit, 'Color', 'r', 'LineWidth', 2, 'HitTest', 'off');
            set(app.linePlot, 'ZData',app.dispersionZVals);
            
            view(app.UIAxes, [0, 90]);
            xlim(app.UIAxes, [min(app.pixels), max(app.pixels)]);
            %ylim(app.UIAxes, [min(app.times), max(app.times)]);
            ylim(app.UIAxes, [-1, 1]);
            try
                clim(app.UIAxes, [-0.003, 0.003]);
            catch
                caxis(app.UIAxes, [-0.003, 0.003]);
            end

            if nargout == 0
                clear app
            end
        end

        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end

    end
end