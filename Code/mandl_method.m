%% Reset
clc
close all
format compact
clearvars

%% Main loop
dataset_names = [ "1979-Mandl", "1991-Baaj", "1994-Shih", "2002-Chakroborty", "2009-Fan(a)", "2009-Fan(b)", "2010-Fan", "2010-Zhang(a)", "2010-Zhang(b)", "2010-Zhang(c)", "2011-Bagloee", "2012-Chew", "2013-Kilic", "2013-Nikolic", "2013-Yan", "2014-Kechagiopoulos", "2014-Kilic(a)", "2014-Kilic(b)", "2014-Nayeem", "2014-Nikolic(a)", "2014-Nikolic(b)", "2015-Amiripour", "2015-Arbex", "2016-Buba", "2018-Buba", "2019-Jha", "2019-Moghaddam", "2020-Capali", "2020-Katsaragakis", "2021-Ahern(a)", "2021-Ahern(b)", "2021-Ahern(d)" ];
dataset_sizes = [ "-4", "-6", "-7", "-8", "-10", "-12" ];

% dataset_names = [ "_Rosentreter" ];

for size = dataset_sizes
    for name = dataset_names

        dataset = "_input/mandls_network/others_solutions/" + name + size;

        if (isfile(dataset + ".csv"))
            evaluate_network(dataset, 0);
        else
%             disp(dataset + ",-,-,-,-,-");
            disp(dataset + ",-");
        end
    end
end

%% Function 
function evaluate_network(dataset, option)
%% Variable Set
num_buses = 20;

%% Import Network - Mandl
nodes = readmatrix("_input/mandls_network/nodes.csv");
demand = readmatrix("_input/mandls_network/demand.csv");
demand_m = zeros(height(nodes));
for i = 1:height(demand)
    demand_m(demand(i, 1), demand(i, 2)) = demand(i, 3);
end
demand = demand_m;

edges = readmatrix("_input/mandls_network/links.csv");
edges_m = zeros(height(nodes)) + 1000;
for i = 1:height(edges)
    edges_m(edges(i, 1), edges(i, 2)) = edges(i, 3);
end
edges = edges_m;

clear demand_m edges_m i

%% Import Routes
routes = [struct("route", [])];

% be sure about this, better to overestimate than underestimate
maxNumberOfElementPerLine = 15;

% build a reading format which can accomodate the longest line
readFormat = repmat('%f',1,maxNumberOfElementPerLine) ;

fidcsv = fopen(dataset + ".csv",'r') ;

M = textscan( fidcsv , readFormat , Inf ,...
    'delimiter',',',...
    'EndOfLine','\r\n',...
    'CollectOutput',true) ;

fclose(fidcsv) ;
M = cell2mat(M) ; % convert to numerical matrix

for i = 1:height(M)
    routes(i) = struct("route", M(i, ~isnan(M(i, :))) + 1);
end


clear i ans M fidcsv maxNumberOfElementPerLine readFormat

%% Route Lengths
lengths = [];

for i = 1:length(routes)
    route = routes(i).route;

    len = 0;
    for j = 1:(length(route) - 1)
        len = len + edges(route(j), route(j + 1));
    end

    lengths = [lengths, len];
end

if (option == 0)
    waiting_time = (sum(lengths) / num_buses) + 0.0001;
else
    waiting_time = 5.0001;
end

clear i j len route

%% Transfer Network
%% SP1: Changing Points
num_routes_connected = zeros(height(nodes), 1);
for i = 1:length(routes)
    route = routes(i).route;

    for j = 1:length(route)
        num_routes_connected(route(j)) = num_routes_connected(route(j)) + 1;
    end
end
mask = num_routes_connected > 1;

I = 1:15;
I = I(mask);

Q = 1:15;
Q = Q(~mask);

clear i j mask route num_routes_connected

%% SP2: Nearest Changing Point
Y_nodes = [];
for i = 1:length(I)
    for j = 1:length(routes)
        route = routes(j).route;

        if (sum(route == I(i)) == 1)
            Y_nodes = [Y_nodes; [I(i), j]];
        end
    end
end

E_edges = [];
for i = 1:length(I)
    for j = 1:length(routes)
        for k = 1:length(routes)
            if (j < k)
                route1 = routes(j).route;
                route2 = routes(k).route;
        
                if (sum(route1 == I(i)) == 1)
                    if (sum(route2 == I(i)) == 1)
                        E_edges = [E_edges; [I(i), j, I(i), k, waiting_time]];
                        E_edges = [E_edges; [I(i), k, I(i), j, waiting_time]];
                    end
                end
            end
        end
    end
end

for i = 1:length(I)
    for l = 1:length(I)
        if (i ~= l)
            for j = 1:length(routes)
                route = routes(j).route;
        
                if (sum(route == I(i)) == 1)
                    if (sum(route == I(l)) == 1)
                        travel_time = 0;
                        
                        i_index = find(route == I(i));
                        l_index = find(route == I(l));

                        o_index = min(i_index, l_index);
                        d_index = max(i_index, l_index);

                        for k = o_index:(d_index - 1)
                            travel_time = travel_time + edges(route(k), route(k + 1));
                        end

                        E_edges = [E_edges; [I(i), j, I(l), j, travel_time]];
                    end
                end
            end
        end
    end
end

nearest_forward = [];
for a = 1:length(Q)
    q = Q(a);

    nearest_forward = [nearest_forward; struct("q", q, "closest", [], "time", 100);];

    for j = 1:length(routes)
        route = routes(j).route;

        if (sum(route == q) == 1)
            
            for b = 1:length(I)
                i = I(b);

                if (sum(route == i) == 1)
                    q_index = find(route == q);
                    i_index = find(route == i);
                    
                    if (q_index < i_index)
                        o_index = min(q_index, i_index);
                        d_index = max(q_index, i_index);
    
                        travel_time = 0;
                        for k = o_index:(d_index - 1)
                            travel_time = travel_time + edges(route(k), route(k + 1));
                        end
    
                        if (travel_time < nearest_forward(end).time)
                            nearest_forward(end) = struct("q", q, "closest", [i, j], "time", travel_time);
                        elseif (travel_time == nearest_forward(end).time)
                            nearest_forward(end).closest = [nearest_forward(end).closest; [i, j]];
                        end
                    end
                end
            end

        end
    end
end

nearest_reverse = [];
for a = 1:length(Q)
    q = Q(a);

    nearest_reverse = [nearest_reverse; struct("q", q, "closest", [], "time", 100);];

    for j = 1:length(routes)
        route = routes(j).route;

        if (sum(route == q) == 1)
            
            for b = 1:length(I)
                i = I(b);

                if (sum(route == i) == 1)
                    q_index = find(route == q);
                    i_index = find(route == i);
                    
                    if (q_index > i_index)
                        o_index = max(q_index, i_index);
                        d_index = min(q_index, i_index);
    
                        travel_time = 0;
                        for k = d_index:(o_index - 1)
                            travel_time = travel_time + edges(route(k), route(k + 1));
                        end
    
                        if (travel_time < nearest_reverse(end).time)
                            nearest_reverse(end) = struct("q", q, "closest", [i, j], "time", travel_time);
                        elseif (travel_time == nearest_reverse(end).time)
                            nearest_reverse(end).closest = [nearest_reverse(end).closest; [i, j]];
                        end
                    end
                end
            end

        end
    end
end


clear i j k route i_index l_index o_index d_index l travel_time route1 route2
clear a b q q_index

%% Changing Index
t1 = E_edges(:, 1) == Y_nodes(:, 1)';
t2 = E_edges(:, 2) == Y_nodes(:, 2)';
t3 = t1 + t2 == 2;
E_1 = find(t3');
for i = 1:length(E_1)
    E_1(i) = E_1(i) - length(Y_nodes) * (i - 1);
end

t1 = E_edges(:, 3) == Y_nodes(:, 1)';
t2 = E_edges(:, 4) == Y_nodes(:, 2)';
t3 = t1 + t2 == 2;
E_2 = find(t3');
for i = 1:length(E_2)
    E_2(i) = E_2(i) - length(Y_nodes) * (i - 1);
end

E_edges_dash = [E_1, E_2, E_edges(:, 5)];

clear E_1 E_2 i t1 t2 t3

%% SP3: Shortest Paths in H
shortest_paths_h = zeros(length(Y_nodes)) + 1000;

for i = 1:length(E_edges_dash)
    shortest_paths_h(E_edges_dash(i, 1), E_edges_dash(i, 2)) = E_edges_dash(i, 3);
end

for i = 1:length(shortest_paths_h)
    shortest_paths_h(i, i) = 0;
end

for k = 1:length(shortest_paths_h)
    for i = 1:length(shortest_paths_h)
        for j = 1:length(shortest_paths_h)
            if (shortest_paths_h(i, j) > shortest_paths_h(i, k) + shortest_paths_h(k, j))
                shortest_paths_h(i, j) = shortest_paths_h(i, k) + shortest_paths_h(k, j);
            end
        end
    end
end

clear i j k

%% SP4: All Shorest Paths
shortest_paths = zeros(height(nodes)) + 1000;

for o = 1:length(nodes)
    for d = 1:length(nodes)
        if (o ~= d)

            % Travel Time in Transfer Network (i -> n(i) -> n(j) -> j)
            origin_on_network = false;
            destination_on_network = false;

            % i -> n(i)
            if (sum(Q == o) == 1)
                origin_index = find([nearest_forward.q] == o);

                origin_edge_list = [];
                if (~isempty(nearest_forward(origin_index).closest))
                    origin_edge_list = [origin_edge_list; nearest_forward(origin_index).closest, nearest_forward(origin_index).time];
                end
                if (~isempty(nearest_reverse(origin_index).closest))
                    origin_edge_list = [origin_edge_list; nearest_reverse(origin_index).closest, nearest_reverse(origin_index).time];
                end
            else
                origin_on_network = true;
            end

            % n(j) -> j
            if (sum(Q == d) == 1)
                destination_index = find([nearest_forward.q] == d);

                destination_edge_list = [];
                if (~isempty(nearest_forward(destination_index).closest))
                    destination_edge_list = [destination_edge_list; nearest_forward(destination_index).closest, nearest_forward(destination_index).time];
                end
                if (~isempty(nearest_reverse(destination_index).closest))
                    destination_edge_list = [destination_edge_list; nearest_reverse(destination_index).closest, nearest_reverse(destination_index).time];
                end
            else
                destination_on_network = true;
            end


            if (origin_on_network + destination_on_network == 2)
                % Search through network to find shorest path
                origin_edge_index = find((Y_nodes(:, 1) == o) == 1);
                destination_edge_index = find((Y_nodes(:, 1) == d) == 1);

                for i = 1:length(origin_edge_index)
                    for j = 1:length(destination_edge_index)
                        travel_time = shortest_paths_h(origin_edge_index(i), destination_edge_index(j));
                        if (travel_time < shortest_paths(o, d))
                            shortest_paths(o, d) = travel_time;
                        end
                    end
                end
            elseif (origin_on_network == 1)
                for i = 1:height(destination_edge_list)
                    destination_edge_index = find(((Y_nodes(:, 1) == destination_edge_list(i, 1)) + (Y_nodes(:, 2) == destination_edge_list(i, 2))) == 2);
%                     origin_edge_index = find(((Y_nodes(:, 1) == o) + (Y_nodes(:, 2) == destination_edge_list(i, 2))) == 2);
%                     origin_edge_index = find((Y_nodes(:, 2) == destination_edge_list(i, 2)) == 1);
                    origin_edge_index = find((Y_nodes(:, 1) == o) == 1);

                    for j = 1:length(origin_edge_index)
                        travel_time = destination_edge_list(i, 3) + shortest_paths_h(origin_edge_index(j), destination_edge_index);
                        if (travel_time < shortest_paths(o, d))
                            shortest_paths(o, d) = travel_time;
                        end
                    end
                end
            elseif (destination_on_network == 1)
                for i = 1:height(origin_edge_list)
                    origin_edge_index = find(((Y_nodes(:, 1) == origin_edge_list(i, 1)) + (Y_nodes(:, 2) == origin_edge_list(i, 2))) == 2);
%                     destination_edge_index = find(((Y_nodes(:, 1) == d) + (Y_nodes(:, 2) == origin_edge_list(i, 2))) == 2);
%                     destination_edge_index = find((Y_nodes(:, 2) == origin_edge_list(i, 2)) == 1);
                    destination_edge_index = find((Y_nodes(:, 1) == d) == 1);

                    for j = 1:length(destination_edge_index)
                        travel_time = origin_edge_list(i, 3) + shortest_paths_h(origin_edge_index, destination_edge_index(j));
                        if (travel_time < shortest_paths(o, d))
                            shortest_paths(o, d) = travel_time;
                        end
                    end
                end
            else
                for i = 1:height(origin_edge_list)
                    for j = 1:height(destination_edge_list)
                        origin_edge_index = find(((Y_nodes(:, 1) == origin_edge_list(i, 1)) + (Y_nodes(:, 2) == origin_edge_list(i, 2))) == 2);
                        destination_edge_index = find(((Y_nodes(:, 1) == destination_edge_list(j, 1)) + (Y_nodes(:, 2) == destination_edge_list(j, 2))) == 2);
                    
                        travel_time = origin_edge_list(i, 3) + shortest_paths_h(origin_edge_index, destination_edge_index) + destination_edge_list(j, 3);
                        if (travel_time < shortest_paths(o, d))
                            shortest_paths(o, d) = travel_time;
                        end
                    end
                end
            end

            % Shorest Time On a Single Route
            for r = 1:length(routes)
                route = routes(r).route;

                if (sum(route == o) == 1)
                    if (sum(route == d) == 1)
                        n1 = find(route == o);
                        n2 = find(route == d);

                        o_index = min(n1, n2);
                        d_index = max(n2, n1);

                        travel_time = 0;
                        for i = o_index:(d_index - 1)
                            travel_time = travel_time + edges(route(i), route(i + 1));
                        end

                        if (travel_time <= shortest_paths(o, d))
                            shortest_paths(o, d) = travel_time;
                        end
                    end
                end
            end
        end
    end
end

shortest_paths = shortest_paths + waiting_time;

clear d d_index destination_edge_index destination_edge_list destination_index destination_on_network
clear o o_index origin_edge_index origin_edge_list origin_index origin_on_network
clear i j n1 n2 r route travel_time

%% Total System Time
% total_time = sum(sum(demand .* (floor(shortest_paths .* 100) ./ 100)));
% 
% if (option == 0)
%     total_time = sum(sum(demand .* shortest_paths));
% else
%     total_time = sum(sum(demand .* (floor(shortest_paths * 1000) / 1000)));
% end

total_time = sum(sum(demand .* (floor(shortest_paths * 1000) / 1000)));

% disp("CTTT: " + num2str(total_time));
% disp("PTTT: " + '260485.85');
% temp = total_time - 260485.85;
% disp("diff: " + "     " + num2str(temp));
% 
% disp("CATT: " + "    " + num2str(total_time / sum(sum(demand))));
% disp("PATT: " + "    "  + num2str(260485.85 / sum(sum(demand))));

%% Percentages
transfers = round(mod(shortest_paths * 100, 1) * 100);

% if (option == 0)
%     transfers = round(mod(shortest_paths * 10, 1) * 10);
% end

d_0 = sum(sum((transfers == 1) .* demand)) / sum(sum(demand));
d_1 = sum(sum((transfers == 2) .* demand)) / sum(sum(demand));
d_2 = sum(sum((transfers == 3) .* demand)) / sum(sum(demand));
d_u = sum(sum((transfers > 3) .* demand)) / sum(sum(demand));

% disp("Cd_0: " + "    " + num2str(round(d_0 * 100 * 100) / 100));
% disp("Cd_1: " + "    " + num2str(round(d_1 * 100 * 100) / 100));
% disp("Cd_2: " + "     " + num2str(round(d_2 * 100 * 100) / 100));
% disp("Cd_u: " + "     " + num2str(round(d_u * 100 * 100) / 100));
% disp("CATT: " + "    " + num2str((total_time / sum(sum(demand)) - waiting_time)));

% disp(num2str(round(d_0 * 100 * 100) / 100 / 100));
% disp(num2str(round(d_1 * 100 * 100) / 100 / 100));
% disp(num2str(round(d_2 * 100 * 100) / 100 / 100));
% disp(num2str(round(d_u * 100 * 100) / 100 / 100));
% disp(num2str((total_time / sum(sum(demand)) - waiting_time)));
% disp (dataset);

% disp(dataset + "," + num2str(waiting_time - 0.0001));

disp(dataset + "," + ...
     num2str(round(d_0 * 100 * 100) / 100 / 100) + "," + ...
     num2str(round(d_1 * 100 * 100) / 100 / 100) + "," + ...
     num2str(round(d_2 * 100 * 100) / 100 / 100) + "," + ...
     num2str(round(d_u * 100 * 100) / 100 / 100) + "," + ...
     num2str((total_time / sum(sum(demand))) - (waiting_time - 0.0001)));

end