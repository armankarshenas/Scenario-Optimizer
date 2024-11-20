% BRCM_DEMOFILE Step by step introduction to the BRCM Toolbox. Heavily commented.
% ------------------------------------------------------------------------
% This file is part of the BRCM Toolbox v1.03.
%
% The BRCM Toolbox - Building Resistance-Capacitance Modeling for Model Predictive Control.
% Copyright (C) 2013  Automatic Control Laboratory, ETH Zurich.
%
% The BRCM Toolbox is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% The BRCM Toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with the BRCM Toolbox.  If not, see <http://www.gnu.org/licenses/>.
%
% For support check www.brcm.ethz.ch.
% ------------------------------------------------------------------------
function [A,Bu,Bv,Bvu,Bxu, B] = dynamic(Th)

global g_debugLvl
g_debugLvl = 1;
BRCMRootPath = getBRCMRootPath();
cd(BRCMRootPath);
% The name of the folder in which the building data is stored
Model = 'Firstexample';
% The address of the excel files
thermalModelDataDir =   [BRCMRootPath,filesep,'BuildingData',filesep,Model,filesep,'ThermalModel'];
EHFModelDataDir =       [BRCMRootPath,filesep,'BuildingData',filesep,Model,filesep,'EHFM'];


%% --------------------------------------------------------------------------------------
% 1) Create a building
% --------------------------------------------------------------------------------------
% Create an empty Building object with an optional identifier argument.
buildingIdentifier = Model;
B = Building(buildingIdentifier);
%% --------------------------------------------------------------------------------------
% 2) Load the thermal model data
% --------------------------------------------------------------------------------------

% Load the thermal model data.
B.loadThermalModelData(thermalModelDataDir);
% The thermal model data consists of zones, building elements, constructions,
% materials, windows and parameters. The data of each element group must
% be provided by a separate .xls files and all base files are required for
% loading the builing data. We require the file names and the file contents to follow a
% specific convention, see the Documentation.
B.checkThermalModelDataConsistency;
%% --------------------------------------------------------------------------------------
% 3) Declare external heat flux models that should be included
% --------------------------------------------------------------------------------------
% Heat exchange with ambient air and solar gains
EHFModelClassFile = 'BuildingHull.m';                                         % This is the m-file defining this EHF model's class.
EHFModelDataFile = [EHFModelDataDir,filesep,'buildinghull'];                  % This is the spreadsheet containing this EHF model's specification.
EHFModelIdentifier = 'BuildingHull';                                          % This string identifies the EHF model uniquely
B.declareEHFModel(EHFModelClassFile,EHFModelDataFile,EHFModelIdentifier);

% Occupancy
EHFModelClassFile = 'InternalGains.m';
EHFModelDataFile = [EHFModelDataDir,filesep,'internalgains'];
EHFModelIdentifier = 'IG';
B.declareEHFModel(EHFModelClassFile,EHFModelDataFile,EHFModelIdentifier);

% Cooling
EHFModelClassFile = 'BEHeatfluxes.m';
EHFModelDataFile = [EHFModelDataDir,filesep,'BEHeatfluxes'];
EHFModelIdentifier = 'TABS';
B.declareEHFModel(EHFModelClassFile,EHFModelDataFile,EHFModelIdentifier);

% Heating
EHFModelClassFile = 'Radiators.m';
EHFModelDataFile = [EHFModelDataDir,filesep,'radiators'];
EHFModelIdentifier = 'Rad';
B.declareEHFModel(EHFModelClassFile,EHFModelDataFile,EHFModelIdentifier);

%% --------------------------------------------------------------------------------------
% 4) Display thermal model data to Command Window and draw Building (optional)
% --------------------------------------------------------------------------------------
% 3-D plot of Building
B.drawBuilding;
% 2-D plot of Building
B.drawBuilding('Floorplan');

%% --------------------------------------------------------------------------------------
% 5) Generate thermal model and full model
% --------------------------------------------------------------------------------------
% Generate thermal model (optional)
B.generateThermalModel;

% Generate (full) building model (includes thermal model generation if not yet done)
B.generateBuildingModel;

% Disretization
B.building_model.setDiscretizationStep(Th);
B.building_model.discretize();

%% --------------------------------------------------------------------------------------
% 6) Retrieve Matrices and generate costs and constraints
% --------------------------------------------------------------------------------------
% Access of full model matrices
A = B.building_model.discrete_time_model.A;
Bu = B.building_model.discrete_time_model.Bu;
Bv = B.building_model.discrete_time_model.Bv;
Bvu = B.building_model.discrete_time_model.Bvu;
Bxu = B.building_model.discrete_time_model.Bxu;

end

