function [ net2_in_net1, net1_to_net2, net1_in_net2, net2_to_net1 ] = map_nuclei( net1_aa, net1_zz, net2_aa, net2_zz )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% net1_to_net2
[net2_in_net1,net1_to_net2] = ismember( [net2_aa,net2_zz], [net1_aa,net1_zz], 'rows' );
% net1_to_net2 = zeros(size(net2_aa));
% k = 1;
% for i = 1:length(net2_aa)
%     j = find( net2_aa(i) == net1_aa & net2_zz(i) == net1_zz );
%     if ~isempty(j)
%         net1_to_net2(k) = j;
%         k = k + 1;
%     end
% end
net1_to_net2(net1_to_net2 == 0) = [];

%% net2_to_net1
[net1_in_net2,net2_to_net1] = ismember( [net1_aa,net1_zz], [net2_aa,net2_zz], 'rows' );
% net2_to_net1 = zeros(size(net1_aa));
% k = 1;
% for i = 1:length(net1_aa)
%     j = find( net1_aa(i) == net2_aa & net1_zz(i) == net2_zz );
%     if ~isempty(j)
%         net2_to_net1(k) = j;
%         k = k + 1;
%     end
% end
net2_to_net1(net2_to_net1 == 0) = [];

end

