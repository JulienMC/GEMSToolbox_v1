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

function Outputs = source_terms(case_st,varargin)
   % Assiging the source terms for the model
   %    currently only the temperature is assigned here
   %    the flow st q is still assigned in geometries
   %
    switch case_st
        %% Temperature Source Terms
        case 0
            % Rock temperature set to Tr everywhere
            [PhysicalProperties, Outputs] = varargin{1:end};
            nn = Outputs.nn;
            no = Outputs.no;
            xtotal = Outputs.xtotal;
            Tr = PhysicalProperties.Tr;
            if length(Tr) > 1
                error("Invalid Tr specified for source term case 0. A single value is expected.")
            end
            N = nn+no;
            % apply Tr to all nodes to start with:
            T_st    = ones(N,1)*Tr;

            Outputs.('T_st') = T_st;

        % Injection point set to Tf_ini, everything else to Tr
        case 1 %BROKEN
            [q_in, q_out, Tr, Tf_ini, no, nn, Ho] = varargin{1:end};

            % Inflow temperature for all the nodes
            T_st = zeros(no+nn,1) + Tr;
            % for injection node override with Tf_ini
            for i = 1:length(q_in)
                T_st(q_in{i} + no) = Tf_ini;
            end
            Outputs.('T_st') = T_st;

        % If head defined model Tf_ini is applied to the highest fixed
        % pressure node.
        case 2 %BROKEN
            [q_in, q_out, Tr, Tf_ini, no, nn, Ho] = varargin{1:end};
            % Inflow temperature for all the nodes
            T_st = zeros(no+nn,1) + Tr;
            % for max pressure node override with Tf_ini
            maxHi = find(Ho == max(Ho)); % index of max pressure node
            T_st(nn+maxHi) = Tf_ini;
            Outputs.('T_st') = T_st;

        % Injection points set to the values specified in the Tf_ini array
        % from the input file
        case 3
            [PhysicalProperties, Outputs] = varargin{1:end};
            Tf_ini = PhysicalProperties.Tf_ini;
            internal_nodes = internalState.nodes;
            T_st = Outputs.Tr_st;
            q_in = Outputs.q_in;
            % apply Tr to all nodes water temperature to start with:

            % check whether Tf_ini is generic or well specific
            if length(Tf_ini) == 1
                for i = 1:length(q_in)
                    T_st(internal_nodes(q_in(i)))  = Tf_ini;
                end
            else % well specific inflow temperatures specified
                for i = 1:length(q_in)
                    T_st(internal_nodes(q_in(i)))  = Tf_ini(i);
                end
            end
            Outputs.('T_st') = T_st;

        % non-injection ndoes set using geothermal gradient specified in the Tf_ini array
        % from the input file
        case 5
            [PhysicalProperties, Outputs] = varargin{1:end};
            Tr = PhysicalProperties.Tr;
            nn = Outputs.nn;
            no = Outputs.no;
            xtotal = Outputs.xtotal;
            N = nn+no;

            % check whether Tr is array or scalar
            if length(Tr) == 1
                warning("A single value is specified for the geothermal gradient under the input keyword Tr. It is assumed to be the gradient and a surface temperature of 5C will be assumed.");
                Tgrad = Tr;
                Tsurf = 5;
            else % well specific inflow temperatures specified
                Tgrad = Tr(1);
                Tsurf = Tr(2);
            end

            % check if xtotal is 2D or 3D and adds a blank 1 column if
            % required, so that the rock temperature will be set to the
            % Tgrad + Tsurf value
            if size(xtotal,2) < 3
                xtotal(:,3) = 0;
            end

            % apply the geothermal gradient to all the nodes
            Tr_st  = Tgrad * xtotal(:,3) + Tsurf;
            Outputs.('Tr_st') = Tr_st;

        %% Groundwater flow Source Terms
        case 4
            [PhysicalProperties, Outputs] = varargin{1:end};
            % Setting groundwater interaction
            V_gw = PhysicalProperties.v_gw;%1e-10; % groundwater darcy velocity m/s;
            nn = Outputs.nn;
            no = Outputs.no;
            q = Outputs.q;
            q_in = Outputs.q_in;
            L = Outputs.L;
            d = Outputs.d;
            A12 = Outputs.A12;
            x = Outputs.x;
            
            A = zeros(nn, 1);
            for i = 1:nn
                pidxs = find(A12(:,i)); % find all pipes connected to node i
                % divide by 2 as each node affects only half of the pipe length connected
                A(i) = sum(L(pidxs)/2*2*pi.*d(pidxs)/2); % Area of pipe exposed to groundwater that needs applying at node i
            end
            q = groundwater_st(q, x, V_gw, A);
            
            Outputs.('q') = q;


        otherwise
            error("Source terms could not be set properly.");
%             % default ST
%             [q_in, q_out, Tr, Tf_ini, no, nn, Ho] = varargin{1:end};
%
%             % Inflow temperature for all the nodes
%             T_st = zeros(no+nn,1) + Tr;
    end

    Outputs = Outputs;
end
