classdef TAExperiment < handle
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here

    properties
       
        json;

        fileName (1,1) string

        nPixels (1, 1) {mustBeNumeric}
        pixels (1, :) {mustBeNonnegative}
        wavelengths(1, :) {mustBeNonnegative}

        Calibration (1, 2) {mustBeNumeric}
        CalibrationPoints (:, 1) CalibrationPoint

        pumpWavelength (1, 1) {mustBeNumeric}
        pumpBandwidth (1, 1) {mustBeNumeric}
        pumpEnergy (1, 1) {mustBeNumeric}

        analytes (1, :) string
        solvent (1, 1) string
        
        nTimes (1, 1) {mustBeNumeric}
        times (:, 1) {mustBeNumeric}

        nScans (1, 1) {mustBeNonnegative}

        scans %(1, :) TAScan

        TAMean (:, :) {mustBeNumeric}
        TAVariance (:, :) {mustBeNumeric}
        TANShots (:, :) {mustBeNumeric}

        pumpOnMean (:, :) {mustBeNumeric}
        pumpOnVariance (:, :) {mustBeNumeric}
        pumpOnNShots (:, :) {mustBeNumeric}

        pumpOffMean (:, :) {mustBeNumeric}
        pumpOffVariance (:, :) {mustBeNumeric}
        pumpOffNShots (:, :) {mustBeNumeric}

        dispersionFitCoefficients (1, :) {mustBeNumeric}
        dispersionFit (1, :) {mustBeNumeric}
        
    end

    methods (Access = public)
        
        function obj = TAExperiment(fileName, options)
            
            arguments
                fileName (1, 1) string
                options.Experiment = 'Burritos';
            end

            obj.fileName = fileName;
            
            if isfile(obj.fileName)
                switch lower(options.Experiment)
                    case 'burritos'
                        obj.loadFile();
                    case 'specter'
                        obj.loadFile716();
                end
            else
                error('File not found');
            end

        end

        function subtractTABackground(obj, upperLimit)

            index = obj.findTime(upperLimit);

            for i = 1:obj.nScans

                background = mean(obj.scans(i).TAMean(1:index, :), 1);

                for j = 1:obj.nTimes

                    obj.scans(i).TAMean(j, :) = obj.scans(i).TAMean(j, :) - background;

                end
                
            end

            obj.combine_scans();
        end

        function correctDispersion(obj)

            correctionApp = DispersionCorrection(obj.pixels, obj.times, obj.TAMean);
            uiwait(correctionApp.UIFigure);
            obj.dispersionFitCoefficients = correctionApp.dispersionFitCoefficients;
            obj.dispersionFit = correctionApp.dispersionFit;

            [~, timeGrid] = meshgrid(obj.pixels, obj.times);
            time_shifted_grid = timeGrid - obj.dispersionFit; % Care must be taken that the dispersion curve is a row vec while the timegrid has columns as pixels amnd the rows as times
            
            for i = 1:obj.nScans
                for j = 1:obj.nPixels
                    obj.scans(i).TAMean(:, j) = interp1(time_shifted_grid(:, j), obj.scans(i).TAMean(:, j), obj.times, 'linear', 'extrap');
                end
            end

            obj.combine_scans();


        end

        function [figTA, axTA] = peekTA(obj, options)
            arguments
                obj TAExperiment
                options.XValues = "Wavelengths"

            end

            figTA = figure();
            axTA = axes('Parent',figTA);
            if options.XValues == "Wavelengths"
                surf(axTA, obj.wavelengths, obj.times, obj.TAMean, 'EdgeColor', 'none');
                xlabel(axTA, 'Wavelengths / nm');
                xlim(axTA, [min(obj.wavelengths), max(obj.wavelengths)]);
            else
                surf(axTA, obj.pixels, obj.times, obj.TAMean, 'EdgeColor', 'none');
                xlabel(axTA, 'Pixels');
                xlim(axTA, [min(obj.pixels), max(obj.pixels)]);
            end
            ylabel(axTA, 'Times / ps');
            ylim(axTA, [min(obj.times), max(obj.times)]);
            view(axTA, [0, 90]);
            title(axTA, obj.makeTitle());
            subtitle(axTA, obj.makeSubtitle());

        end

        function [timeTraceMean, timeTraceVar, timeTraceNshots] = getTimeTrace(obj, wavelength, options)
            arguments
                obj TAExperiment
                wavelength (1,:)
                options.dataType = 'TA';
            end

            switch options.dataType
                % Use function for all data types default to TA
                case 'TA'
                    datamean = obj.TAMean;
                    datavar = obj.TAVariance;
                    dataNShots = obj.TANShots;
                case 'Pump On'
                    disp('Here');
                    datamean = obj.pumpOnMean;
                    datavar = obj.pumpOnVariance;
                    dataNShots = obj.pumpOnNShots;
                case 'Pump Off'
                    datamean = obj.pumpOffMean;
                    datavar = obj.pumpOffVariance;
                    dataNShots = obj.pumpOffNShots;
            end
            
            if length(wavelength) == 1
                % Return single wavelength trace
                index = obj.findWavelength(wavelength);
                timeTraceMean = datamean(:, index);
                timeTraceVar = datavar(:, index);
                timeTraceNshots = dataNShots(:, index);
                
            elseif length(wavelength) == 2
                indices = [obj.findWavelength(wavelength(1,1)), obj.findWavelength(wavelength(1,2))];

                minIndex = min(indices);
                maxIndex = max(indices);

                timeTraceNshots = sum(dataNShots(:, minIndex:maxIndex), 2);
                timeTraceMean = sum(dataNShots(:, minIndex:maxIndex).*datamean(:, minIndex:maxIndex), 2, 'omitnan')./timeTraceNshots;
                timeTraceVar = sum(dataNShots(:, minIndex:maxIndex).*datavar(:, minIndex:maxIndex), 2, 'omitnan')./timeTraceNshots;
                
            else
                error('Wavelength argument must be a single wavelength or an array of two represneing upper and lower limits');
            end


        end

        function [wavelengthTraceMean, wavelengthTraceVar, wavelengthTraceNshots] = getWavelengthTrace(obj, time, options)
            arguments
                obj TAExperiment
                time (1,:)
                options.dataType = 'TA';
            end

            switch options.dataType
                % Use function for all data types default to TA
                case 'TA'
                    datamean = obj.TAMean;
                    datavar = obj.TAVariance;
                    dataNShots = obj.TANShots;
                case 'Pump On'
                    disp('Here');
                    datamean = obj.pumpOnMean;
                    datavar = obj.pumpOnVariance;
                    dataNShots = obj.pumpOnNShots;
                case 'Pump Off'
                    datamean = obj.pumpOffMean;
                    datavar = obj.pumpOffVariance;
                    dataNShots = obj.pumpOffNShots;
            end
            
            if length(time) == 1
                % Return single wavelength trace
                index = obj.findTime(time);
                wavelengthTraceMean = datamean(index, :);
                wavelengthTraceVar = datavar(index, :);
                wavelengthTraceNshots = dataNShots(index, :);
                
            elseif length(time) == 2
                indices = [obj.findTime(time(1,1)), obj.findTime(time(1,2))];

                minIndex = min(indices);
                maxIndex = max(indices);

                wavelengthTraceNshots = sum(dataNShots(minIndex:maxIndex, :), 1);
                wavelengthTraceMean = sum(dataNShots(minIndex:maxIndex, :).*datamean(minIndex:maxIndex, :), 1, 'omitnan')./wavelengthTraceNshots;
                wavelengthTraceVar = sum(dataNShots(minIndex:maxIndex, :).*datavar(minIndex:maxIndex, :), 1, 'omitnan')./wavelengthTraceNshots;
                
            else
                error('Wavelength argument must be a single wavelength or an array of two represneing upper and lower limits');
            end


        end

    end

    methods (Access = private)

        function index = findWavelength(obj, wavelength)
            [~, index] = min(abs(obj.wavelengths - wavelength));
        end

        function index = findTime(obj, time)
            [~, index] = min(abs(obj.times - time));
        end

        function loadFile(obj)

            data = jsondecode(fileread(obj.fileName));
            obj.json = data;
            
            obj.nPixels = data.Spectrometer.num_pixels;
            
            obj.buildPixelsVector();
            
            obj.Calibration(1, 1) = data.Spectrometer.wavelength_calibration.slope;
            obj.Calibration(1, 2) = data.Spectrometer.wavelength_calibration.intercept;

            obj.buildWavelengthVector;

            if isempty(data.Spectrometer.calibration_points)
                warning('No calibration points in data file.');
            else
                for i = 1:length(data.Spectrometer.calibration_points)
                    obj.CalibrationPoints = [obj.CalibrationPoints; CalibrationPoint(data.Spectrometer.calibration_points(i).wavelength, data.Spectrometer.calibration_points(i).pixel, data.Spectrometer.calibration_points(i).uncertainty)];
                end
            end

            obj.pumpWavelength = data.TransientAbsorption.pump_wavelength;
            obj.pumpBandwidth = data.TransientAbsorption.pump_bandwidth;
            %obj.pumpEnergy = data.TransientAbsorption.pump_energy;
            obj.solvent = data.TransientAbsorption.solvent;

            obj.nTimes = length(data.TransientAbsorption.time_delays);
            obj.times = data.TransientAbsorption.time_delays;

            obj.nScans = length(data.TransientAbsorption.scans);
            obj.scans = TAScan.empty(obj.nScans, 0);

            for i = 1:length(data.TransientAbsorption.analytes)
                mol = data.TransientAbsorption.analytes(i).analyte;
                conc = sprintf('%.3f', data.TransientAbsorption.analytes(i).concentration);
                anal = {mol, conc};
                obj.analytes = [obj.analytes, anal];
            end

            for i = 1:obj.nScans
                obj.scans(i) = TAScan(data.TransientAbsorption.scans(i), obj.nTimes, obj.nPixels);
            end

            [obj.TAMean, obj.TAVariance, obj.TANShots, obj.pumpOnMean, obj.pumpOnVariance, obj.pumpOnNShots, obj.pumpOffMean, obj.pumpOffVariance, obj.pumpOffNShots] = deal(zeros(obj.nTimes, obj.nPixels));
            
            obj.combine_scans();

        end

        function loadFile716(obj)

            listTimes = [];

            lines = readlines(obj.fileName);
            
            obj.nPixels = 256;
            obj.buildPixelsVector();

            listNShots = [];
            listPumpOff = [];
            listPumpOn = [];
            listTAMean = [];
            listTAVar = [];

            for i = 1:length(lines)

                line = lines(i);

                switch line
                    case "# Spectrograph:  Slope and intercept"
                        slopeAndIntercept = split(lines(i+1));
                        obj.Calibration(1, 1) = str2double(slopeAndIntercept(1));
                        obj.Calibration(1, 2) = str2double(slopeAndIntercept(2));
                        obj.buildWavelengthVector();
                    case "# Pump-Probe Delay"
                        listTimes = [listTimes; str2double(lines(i+1))];
                        listNShots = [listNShots; str2double(split(lines(i+4)))'];
                        listPumpOff = [listPumpOff; str2double(split(lines(i+7)))'];
                        listPumpOn = [listPumpOn; str2double(split(lines(i+13)))'];
                        listTAMean = [listTAMean; 1E3*str2double(split(lines(i+13)))'];
                        listTAVar = [listTAVar; (1E3*str2double(split(lines(i+13)))).^2'];
                end

                
            end
            
            obj.times = unique(listTimes);
            obj.nTimes = length(obj.times);

            obj.nScans = ceil(length(listTimes) / length(obj.times));

            pad = zeros(obj.nScans*length(obj.times) - length(listTimes), obj.nPixels);

            listNShots = [listNShots; pad];
            listPumpOff = [listPumpOff; pad];
            listPumpOn = [listPumpOn; pad];
            listTAMean = [listTAMean; pad];
            listTAVar = [listTAVar; pad];

            cubeNShots = reshape(listNShots, obj.nScans, obj.nTimes, obj.nPixels);
            cubePumpOff = reshape(listPumpOff, obj.nScans, obj.nTimes, obj.nPixels);
            cubePumpOn = reshape(listPumpOn, obj.nScans, obj.nTimes, obj.nPixels);
            cubeTAMean = reshape(listTAMean, obj.nScans, obj.nTimes, obj.nPixels);
            cubeTAVar = reshape(listTAVar, obj.nScans, obj.nTimes, obj.nPixels);

           %%%%%% NEED TO REFACTOR HOW INITIALISING TASCAN WORKS
            

        end

        function buildPixelsVector(obj)
            obj.pixels = linspace(1, obj.nPixels, obj.nPixels);
        end

        function buildWavelengthVector(obj)
            obj.wavelengths = obj.pixels * obj.Calibration(1,1) + obj.Calibration(1, 2);
        end

        function combine_scans(obj)
            % Compute a weightd mean of the data from the different scans.
            % Check both the mean and variance in each scan for nan values
            % and exclude that point in that scan for the averaging
            % Pooled variance = weighted average of variance

            [TAMeanCube, TAVarCube, TANShotsCube, pumpOnMeanCube, pumpOnVarCube, pumpOnNShotsCube, pumpOffMeanCube, pumpOffVarCube, pumpOffNShotsCube] = deal(zeros(obj.nScans, obj.nTimes, obj.nPixels));

            TANShotsCube(isnan(TAMeanCube + TAVarCube)) = 0;
            pumpOnNShotsCube(isnan(pumpOnMeanCube + pumpOnVarCube)) = 0;
            pumpOffNShotsCube(isnan(pumpOffMeanCube + pumpOffVarCube)) = 0;

            for i = 1:obj.nScans
                TAMeanCube(i, :, :) = obj.scans(i).TAMean;
                TAVarCube(i, :, :) = obj.scans(i).TAVariance;
                TANShotsCube(i, :, :) = obj.scans(i).TANShots;
            end

            obj.TAMean = reshape(sum(TANShotsCube.*TAMeanCube, 1, 'omitnan')./sum(TANShotsCube, 1), obj.nTimes, obj.nPixels);
            obj.TAVariance = reshape(sum(TANShotsCube.*TAVarCube, 1, 'omitnan')./sum(TANShotsCube, 1), obj.nTimes, obj.nPixels);
            obj.TANShots = reshape(sum(TANShotsCube, 1), obj.nTimes, obj.nPixels);
            obj.pumpOnMean = reshape(sum(pumpOnNShotsCube.*pumpOnMeanCube, 1, 'omitnan')./sum(pumpOnNShotsCube, 1), obj.nTimes, obj.nPixels);
            obj.pumpOnVariance = reshape(sum(pumpOnNShotsCube.*pumpOnVarCube, 1, 'omitnan')./sum(pumpOnNShotsCube, 1), obj.nTimes, obj.nPixels);
            obj.pumpOnNShots = reshape(sum(pumpOnNShotsCube, 1), obj.nTimes, obj.nPixels);
            obj.pumpOffMean = reshape(sum(pumpOffNShotsCube.*pumpOffMeanCube, 1, 'omitnan')./sum(pumpOffNShotsCube, 1), obj.nTimes, obj.nPixels);
            obj.pumpOffVariance = reshape(sum(pumpOffNShotsCube.*pumpOffVarCube, 1, 'omitnan')./sum(pumpOffNShotsCube, 1), obj.nTimes, obj.nPixels);
            obj.pumpOffNShots = reshape(sum(pumpOffNShotsCube, 1), obj.nTimes, obj.nPixels);

            % Perform a final check for any points which were a nan across
            % al scans (typical of a portion of the spectrometer upon which
            % no white light falls and set them to zero for future
            % algorithms.
            %obj.TAMean(isnan(obj.TAMean)) = 0;
            obj.TAVariance(isnan(obj.TAVariance)) = 0;
            obj.pumpOnMean(isnan(obj.pumpOnMean)) = 0;
            obj.pumpOnVariance(isnan(obj.pumpOnVariance)) = 0;
            obj.pumpOffMean(isnan(obj.pumpOffMean)) = 0;
            obj.pumpOffVariance(isnan(obj.pumpOffVariance)) = 0;

        end

        function titlestring = makeTitle(obj)
            titlestring = "";
            for i = 1:size(obj.analytes, 1)
                analytestring = strcat(obj.analytes(i, 1), " (", obj.analytes(i, 2), " \muM)  ");
                titlestring = strcat(titlestring, analytestring);
            end

            titlestring = strcat(titlestring, '/ ', obj.solvent);
        end

        function subtitlestring = makeSubtitle(obj)
            wavelength = sprintf('%.0f', obj.pumpWavelength);
            bandwidth = sprintf('%.0f', obj.pumpBandwidth);
            %energy = sprintf('%.1f', obj.pumpEnergy);
            subtitlestring = strcat("\lambda_{exc} = ", wavelength, " nm \Delta_{\lambda} = ", bandwidth, " nm");
            %subtitlestring = strcat("\lambda_{exc} = ", wavelength, " nm \Delta_{\lambda} = ", bandwidth, " nm E = ", energy, " nJ");
        end

    end

end
