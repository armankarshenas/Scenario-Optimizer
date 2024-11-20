function check_file_extension(ext,fileXLS,type) %#ok<INUSD>
   %CHECK_FILE_EXTENSION Checks whether the file has the appropriate/supported extension.
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
   
   
   % This assignement required since it does not work in sprintf below, but from command line sprintf works great. why?
   % remove * from supported extensions
   
   
   extensions = cellfun(@(x) x(2:end),Constants.supported_file_extensions,'UniformOutput',0);
   
   if isempty(intersect(ext,extensions))
      error('File:Extension','File ''%s'' has incompatible extension.\nSUPPORTED exensions:\t %s\n',...
         sprintf('%s',regexprep(fileXLS,'\','\\')),sprintf('''%s'' ',extensions{:}));
   end
end