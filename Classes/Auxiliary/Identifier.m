classdef Identifier
   %IDENTIFIER This class stores all the identifiers of a model.
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
   

   
   
   
   properties
      x@cell = {};            % Identifiers for the x-state
      q@cell = {};            % Identifiers for the q-state
      u@cell = {};            % Identifiers for the u control inputs
      v@cell = {};            % Identifiers for the v disturbance inputs
      y@cell = {};            % Identifiers for the y outputs
      constraints@cell = {};  % Identifiers for the constraints
   end
   
   methods
      %constructor
      function obj = Identifier()
         
      end % Identifier
   end % methods
end %classdef
