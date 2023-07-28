% This program allows for the computation of water and heat flow through
% a mine network
%     Copyright (C) 2022  Durham University
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

%% DESCRIPTION
% This function will assign the well flow rates

function q = setupwells(qset, q_in, q_out, nn)
    % Checks for invalid q_in or q_out
    for i = 1:length(q_in)
        if internalState.nodes("get",q_in(i)) > nn
            error('q_in and/or q_out locations are greater than nodal space of mine model. Choose different q_in and/or q_out.');
        end
    end

    for i = 1:length(q_out)
        if internalState.nodes("get",q_out(i)) > nn
            error('q_in and/or q_out locations are greater than nodal space of mine model. Choose different q_in and/or q_out.');
        end
    end

    % set any external in/outflow for each (non-fixed) node:
    q    = zeros(nn,1);

    % check whether qset is generic or well specific
    if length(qset) == 1
        fprintf("Assigning injection flow of %f m^3.s^-1 to all node(s): [%s]\nAssigning abstraction flow of  %f to all node(s): [%s]\n", -qset, strjoin(string(q_in)), qset, strjoin(string(q_out)))
        for i = 1:length(q_in)
            q(internalState.nodes("get",q_in(i)))  = -qset;
        end
        for i = 1:length(q_out)
            q(internalState.nodes("get",q_out(i))) = qset;
        end
    else % well specific flow rates summing to 0 are specified
        if abs(length(q_in) + length(q_out) - length(qset)) > 0
            error("Mismatch between specified qset and q_in and q_out inputs. Please make sure that the number of q_in and q_out nodes specified is equal to the number of qset values specified.")
        end
        for i = 1:length(q_in)
            fprintf("Assinging injection node %d flow of %f m^3.s^-1\n", q_in(i), -qset(1,i))
            q(internalState.nodes("get",q_in(i)))  = -qset(1,i);
        end
        for i = 1:length(q_out)
            fprintf("Assinging abstraction node %d flow of %f m^3.s^-1\n", q_out(i), qset(1,end+1-i))
            q(internalState.nodes("get",q_out(i))) = qset(1,length(q_in)+i);
        end
    end
end