classdef CalibrationPoint < handle

    properties
        wavelength (1,1) {mustBeNumeric}
        pixel (1,1) {mustBeNumeric}
        variance (1, 1) {mustBeNumeric}
    end
    
    methods
        function obj = CalibrationPoint(wavelength, pixel, variance)

            obj.wavelength = wavelength;
            obj.pixel = pixel;
            obj.variance = variance;
        end

    end
end

