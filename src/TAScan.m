classdef TAScan < handle
    properties

        dateStarted (1, 1) string
        timeStarted (1, 1) string
        dateFinished (1, 1) string
        timeFinished (1, 1) string

        TAMean (:, :) {mustBeNumeric}
        TAVariance (:, :) {mustBeNumeric}
        TANShots (:, :) {mustBeNumeric}

        pumpOnMean (:, :) {mustBeNumeric}
        pumpOnVariance (:, :) {mustBeNumeric}
        pumpOnNShots (:, :) {mustBeNumeric}

        pumpOffMean (:, :) {mustBeNumeric}
        pumpOffVariance (:, :) {mustBeNumeric}
        pumpOffNShots (:, :) {mustBeNumeric}
        

    end
    
    methods
        function obj = TAScan(jsonscan, nTimes, nPixels)
            obj.dateStarted = jsonscan.date_started;
            obj.timeStarted = jsonscan.time_started;
            obj.dateFinished = jsonscan.date_finished;
            obj.timeFinished = jsonscan.time_finished;

            [obj.TAMean, obj.TAVariance, obj.TANShots, obj.pumpOnMean, obj.pumpOnVariance, obj.pumpOnNShots, obj.pumpOffMean, obj.pumpOffVariance, obj.pumpOffNShots] = deal(zeros(nTimes, nPixels));
            
            for i = 1:nTimes
                obj.TAMean(i, :) = jsonscan.scan(i).spectrum.transient_absorption.mean;
                obj.TAVariance(i, :) = jsonscan.scan(i).spectrum.transient_absorption.variance;
                obj.TANShots(i, :) = jsonscan.scan(i).spectrum.transient_absorption.num_shots;

                obj.pumpOnMean(i, :) = jsonscan.scan(i).spectrum.pump_on.mean;
                obj.pumpOnVariance(i, :) = jsonscan.scan(i).spectrum.pump_on.variance;
                obj.pumpOnNShots(i, :) = jsonscan.scan(i).spectrum.pump_on.num_shots;

                obj.pumpOffMean(i, :) = jsonscan.scan(i).spectrum.pump_off.mean;
                obj.pumpOffVariance(i, :) = jsonscan.scan(i).spectrum.pump_off.variance;
                obj.pumpOffNShots(i, :) = jsonscan.scan(i).spectrum.pump_off.num_shots;
            end
        end
    end
    
end