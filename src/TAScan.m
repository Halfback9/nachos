classdef TAScan < handle

    % TAScan - Data corresponding to a single transient absorption scan
    %   Upon initialisation (with number of pixels and number of time
    %   points) creates blank matriced which correspond to the mean,
    %   variance and number of shots for the transient absorption, pump off
    %   and pump on signals
    %
    %
    %
    % TAScan Properties:
    %
    %
    %
    %   dateStarted                         -   Date upon which the scan
    %                                           started.
    %
    %   timeStarted                         -   Time at which the scan was
    %                                           started.
    %
    %   dateFinished                        -   Date upon which the scan
    %                                           finished.
    %
    %   timeFinished                        -   Time at which the scan
    %                                           finished.
    %
    %
    %
    %   nTimes                              -   Number of time delays
    %                                           collected.
    %
    %   nPixels                             -   Number of pixels measured.
    %
    %
    %
    %   TAMean                              -   A matrix containing the
    %                                           measured transient
    %                                           absorption signals.
    %
    %   TAVariance                          -   A matrix containing the
    %                                           variance in the measured
    %                                           transient absorption
    %                                           signals.
    %
    %   TANShots                            -   The number of shots over
    %                                           which the mean and variance
    %                                           of the transient absorption
    %                                           signals was measured.
    %
    %
    %
    %   pumpOnMean                          -   A matrix containing the
    %                                           measured pump on white
    %                                           light probe signals.
    %
    %   pumpOnVariance                      -   A matrix containing the
    %                                           variance in the measured
    %                                           pump on white light probe
    %                                           signals.
    %
    %   pumpOnNShots                        -   The number of shots over
    %                                           which the mean and variance
    %                                           of the pump on white light
    %                                           probe signals were
    %                                           measured.
    %
    %
    %
    %   pumpOffMean                         -   A matrix containing the
    %                                           measured pump off white
    %                                           light probe signals.
    %
    %   pumpOffVariance                     -   A matrix containing the
    %                                           variance in the measured
    %                                           pump off white light probe
    %                                           signals.
    %
    %   pumpOffNShots                       -   The number of shots over
    %                                           which the mean and variance
    %                                           of the pump off white light
    %                                           probe signals were
    %                                           measured.
    %
    %
    %
    % TAScan Methods:
    %
    %
    %
    %   TAScan(nTimes, nPixels)             -   Class constructor.
    %
    %   populate716(jsonscan)               -   Parse data and fill 
    %                                           properties from a scan in
    %                                           716.
    %
    %   populate714(jsonscan)               -   Parse data and fill 
    %                                           properties from a scan in
    %                                           714.
            


    properties

        dateStarted (1, 1) string 
        timeStarted (1, 1) string
        dateFinished (1, 1) string
        timeFinished (1, 1) string

        nTimes (1, 1) {mustBeNumeric}
        nPixels (1, 1) {mustBeNumeric}

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
        function obj = TAScan(nTimes, nPixels)
            % TAScan(nTimes, nPixels) - Class constructor.
            %   Initalises the TAScan class. Requires two arguments, the
            %   number of time delays and number of pixels, in order to
            %   initialises the matrices.

            arguments
                nTimes (1, 1) {mustBeNumeric}
                nPixels (1, 1) {mustBeNumeric}
            end

            obj.nTimes = nTimes;
            obj.nPixels = nPixels;

            [obj.TAMean, obj.TAVariance, obj.TANShots, obj.pumpOnMean, obj.pumpOnVariance, obj.pumpOnNShots, obj.pumpOffMean, obj.pumpOffVariance, obj.pumpOffNShots] = deal(zeros(nTimes, nPixels));
     
        end

        function populate714(obj, jsonscan)
            % populate714(jsonscan) - Parse data and fill properties from a
            % scan in 714.
            %   Takes in a struct derived from decoding json formatted data
            %   produced by the HRRTA experiment in 714. 


            obj.dateStarted = jsonscan.date_started;
            obj.timeStarted = jsonscan.time_started;
            obj.dateFinished = jsonscan.date_finished;
            obj.timeFinished = jsonscan.time_finished;

            for i = 1:obj.nTimes
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

        function populate716(obj, TA, TAVar, TANShots, pumpOff, pumpOn)
            % populate716(jsonscan) - Parse data and fill properties from a
            % scan in 716.
            %   Takes in a set of matrices derived from the reading in the
            %   data produced in 716. That experiment doesn't note the
            %   variance in the white light signals and so the variance is
            %   set to zero and the number of shots is assumed to be the
            %   same as that in the TA.

            obj.TAMean = TA;
            obj.TAVariance = TAVar;
            obj.TANShots = TANShots;
            obj.pumpOffMean = pumpOff;
            obj.pumpOffNShots = TANShots;
            obj.pumpOnMean = pumpOn;
            obj.pumpOnNShots = TANShots;
        end

    end
    
end