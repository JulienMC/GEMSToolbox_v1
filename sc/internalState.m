%% This program allows for the computation of water and heat flow through a mine network
%%     Copyright (C) 2022  Durham University
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
%%

% JMC a class which is used to convert indices from GIS nodes to internal
% Matlab code. This is required to be able to be able to set any node as a
% fixed head node. And/or set the last node as a injection/abstraction
% node.
% The is a static class which stores the relvant data and methods to
% perform the operation described above.
classdef internalState

   methods (Static)
       % Method used to store the GIS indices as values of a vector
       % who's indices correspond to the Matlab indexing used in our model.
       function out = nodes(varargin)
         persistent Var;
         if length(varargin) == 0
             out = Var;
             return;
         end
         % initialisation case
         if varargin{1} == "init"
            Var = [1:varargin{2}+varargin{3}];
         % setting a GIS index to a matlab one
         elseif varargin{1} == "set" % arg{2} is the GIS node, arg{3} is the matlab node
             if varargin{2} > length(Var) || varargin{3} > length(Var)
                    error("Index passed as argument, exceeds the total amount of nodes nn+no. Hence cannot be swapped. Make sure both arguments are in nn+no range.")
             elseif varargin{3} <= internalState.nn()
                    error("You are swapping a node with another unknown head node. If you are trying to assign a GIS node to fixed head node make sure your third argument is > than nn.")
             end
            tempA = find(Var == varargin{3}); % location of the value to be fetched
            tempB = Var(varargin{2}); % value at the location to be replaced
            Var(varargin{2}) = varargin{3};
            Var(tempA) = tempB; % replacing the value at the fetch location with the value that was at the replaced location
            internalState.swaps(varargin{2},tempA); % ID of the swaps are memorised
         % fetching of a GIS index at the matlab location
         elseif varargin{1} == "get" && length(varargin) > 1
            out = internalState.nestedFunctionCall(Var, varargin{2:end});
         else
            out = Var;
         end
       end
        
       % stores the number of unknown nodes locally
       function out = nn(varargin)
         persistent Var;
         if nargin > 0
             Var = varargin{1};
         end
            out = Var;
       end

       % stores the number of known nodes locally
       function out = no(varargin)
         persistent Var;
         if nargin > 0
             Var = varargin{1};
         end
            out = Var;
       end

      % Stores the swaps done so they can be easily accessed and/or reversed for outputting
      % results, or any other reason. Assumes each swap is a fixed head
      % assignment !
      function out = swaps(varargin)
         persistent GISi;
         persistent MLi;
         if nargin > 0 && ~isnumeric(varargin{1}) && varargin{1} == "clear"
             GISi = {};
             MLi = {};
             return
         end

         if nargin > 0 % Add the swap
             newi = length(GISi)+1;
             GISi{newi} = varargin{1};
             MLi{newi} = varargin{2};
         end
            out = [GISi; MLi]; % return the array of the swaps
      end


        % Init method. To be called in the code before setting up fixed
        % heads.
        function Init(nn,no)
           internalState.swaps("clear");
           internalState.nn(nn);
           internalState.no(no);
           internalState.nodes("init",nn,no);
        end

        % Wrapper used to easily perform the head assignement without
        % needing to specify which internal fixed head node to assign the
        % head to. This is computed by the number of swaps (i.e.
        % assignments) already performed.
        function Ho = SetAsFixedHead(GISnode, head, Ho)
            if length(Ho) ~= internalState.no
                error("Argument Ho length mismatches the number of fixed head nodes, no, stored internally to internalState. Make sure internalState.Init(nn, no) has been called, or check the length of Ho.")
            end            
            s = size(internalState.swaps);
            nOfFixH = s(2);
            if nOfFixH + 1 > internalState.no
                error("Trying to assign a fixed head to more nodes than the number, no, of fixed head nodes specified. Increase no, or reduce the number of fixed head node assignments.")
            end
            internalState.nodes("set", GISnode, internalState.nn+nOfFixH+1); % sets the GIS node to the fixed head node internally
            Ho(nOfFixH+1) = head;
        end

        % Method used to convert the global incidence matrix M120 into the
        % fixed head node incidence matrix M10 and the unknown head nodes
        % incidence matrix M12. Not this can also be used for any other
        % node data (e.g. xtotal, or Tn); FIXED HEAD NODES ARE MOVED AT THE
        % END OF M120
       function [M12, M10] = MatSetup(M120, varagin)
           if nargin > 1
                M120 = M120.';                
           end
           % Note the following might missbehave if successive swaps affect the same positions
           % performing all the swaps
           temp = M120;
           M120(:, [internalState.swaps{2,:}]) = M120(:, [internalState.swaps{1,:}]);
           M120(:, [internalState.swaps{1,:}]) = temp(:, [internalState.swaps{2,:}]);
            % splitting the matrix into the unknown and known matrices
           M12 = M120(:, [1:internalState.nn]);
           M10 = M120(:,[internalState.nn+1:end]);

           if nargin > 1
                M12 = M12.'; 
                M10 = M10.'; 
           end
       end

      % Method who reverts the swaps back to the GIS coordinates.
      function [M12, M10] = MatToGIS(M12, M10, varagin)       
           if nargin > 2
                M12 = M12.';
                M10 = M10.';
           end
            M120 = [M12, M10];
           % Note the following might missbehave if successive swaps affect the same positions
           % performing all the swaps
           temp = M120;

           s = size(internalState.swaps);
           nOfFixH = s(2);
           for i = 1:nOfFixH
               GISi = internalState.swaps{1,i};
               MLi = internalState.nodes("get",GISi);
               M120(:, [GISi]) = M120(:, [MLi]);
               M120(:, [MLi]) = temp(:, [GISi]);
           end
            % splitting the matrix into the unknown and known matrices
           M12 = M120(:, [1:internalState.nn]);
           M10 = M120(:,[internalState.nn+1:end]);

           if nargin > 1
                M12 = M12.'; 
                M10 = M10.'; 
           end
      end
        
      % conveinence function to perform nested function calls, a la Python.
       function out = nestedFunctionCall(x, varargin)
           if sum(cell2mat(varargin) > length(x))
               error("Requested node(s) %s, exceeds the total number of nodes in the model: %d. Please check your model size, your inflow outflow points, and any fixed head value node.\n",string(strjoin(string(varargin))),length(x))
           end
            out = x(varargin{:});
       end
   end

end

%% OLD STUFF Potentially not needed unless we want to also use nn and no specific tracking of indices
%        function out = indicesNO(varargin)
%          persistent Var;
%          if length(varargin) == 0
%              out = Var;
%              return;
%          end
%          if varargin{1} == "init"
%             Var = [1:varargin{2}];
%          elseif varargin{1} == "set"  % arg{2} is the GIS node, arg{3} is the matlab node
%             if varargin{2} > length(Var) || varargin{3} > length(Var)
%                        error("Index passed as argument, exceeds the total amount of nodes nn+no. Hence cannot be swapped. Make sure both arguments are in nn+no range.")
%             end   
%             tempA = find(Var == varargin{3}); % location of the value to be fetched
%             tempB = Var(varargin{2}); % value at the location to be replaced
%             Var(varargin{2}) = varargin{3}; % replacing the value at the replace loc with the fectched value
%             Var(tempA) = tempB; % replacing the value at the fetch location with the value that was at the replaced location
%          elseif varargin{1} == "get" && length(varargin) > 1
%             out = internalState.nestedFunctionCall(Var, varargin{2:end});
%          else
%             out = Var;
%          end
%        end

%        function out = indicesNN(varargin)
%          persistent Var;
%          if length(varargin) == 0
%              out = Var;
%              return;
%          end
%          if varargin{1} == "init"
%             Var = [1:varargin{2}];
%          elseif varargin{1} == "set"  % arg{2} is the GIS node, arg{3} is the matlab node
%              if varargin{2} > length(Var) || varargin{3} > length(Var)
%                        error("Index passed as argument, exceeds the total amount of nodes nn+no. Hence cannot be swapped. Make sure both arguments are in nn+no range.")
%              end    
%             tempA = find(Var == varargin{3}); % location of the value to be fetched
%             tempB = Var(varargin{2}); % value at the location to be replaced
%             Var(varargin{2}) = varargin{3};
%             Var(tempA) = tempB; % replacing the value at the fetch location with the value that was at the replaced location
%          elseif varargin{1} == "get" && length(varargin) > 1
%             out = internalState.nestedFunctionCall(Var, varargin{2:end});
%          else
%             out = Var;
%          end
%        end