classdef Analyte < handle
    % Analyte - Information about an analyte
    
    properties (Access = public)
        name string = "";
        IUPACName string = "";
        molecularFormula string = "";
        concentration (1, 1) {mustBeNumeric} % Concentration in micromoles
        electronicAbsorption (:, 1) {mustBeNumeric}
        idCAS string = "";
        idCID string = "";
        molecularWeight (1,1) {mustBeNumeric}

        % Temp
        jsonblah
       
    end

    methods

        function obj = Analyte(name, concentration)
            arguments
                name string = "";
                concentration (1, 1) {mustBeNumeric} = 1000;
            end
            
            obj.name = name;
            obj.concentration = concentration;

            

        end

        function findOnPubChem(obj)

            strHTTP = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/';
            strProperties = '/property/IUPACName,MolecularFormula,MolecularWeight/JSON';
            query = strcat(strHTTP, obj.name, strProperties);
            try
                response = webread(query);
    
                obj.IUPACName = response.PropertyTable.Properties.IUPACName;
                obj.molecularFormula = response.PropertyTable.Properties.MolecularFormula;
                obj.molecularWeight = str2double(response.PropertyTable.Properties.MolecularWeight);
                disp("Added information from PubChem database");
            catch
                disp("Molecule not found on PubChem database.");
            end

        end

    end

    methods % Setters

        function set.name(obj, name)
            obj.name = name;
            obj.findOnPubChem();
        end
        
    end
end