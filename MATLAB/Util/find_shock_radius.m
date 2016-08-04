function [ shock_radius, shock_theta ] = find_shock_radius( radius, theta, bounce_index, p_shocked, max_bins )
% Estimate the shock radius

shock_radius = zeros(size(radius,1),max_bins);
shock_theta  = zeros(size(radius,1),max_bins);
% tic
for i = bounce_index:size(radius,1)
    if length(unique(theta(i,p_shocked{i}))) < max_bins
        nbins = length(unique(theta(i,p_shocked{i})));
    else
        nbins = max_bins;
    end
%     [theta_sort, ~, ~] = unique(theta(i,p_shocked{i}),'sorted');
    [theta_sort, ~] = sort(theta(i,p_shocked{i}));
    edges = interp1(1:length(theta_sort)+2,[0,theta_sort,pi],linspace(1,length(theta_sort)+2,nbins+1));
    
    for j = 1:nbins
        shock_theta(i,j) = (edges(j) + edges(j+1)) / 2;
    end
 
    [~,bin]=histc(theta(i,p_shocked{i}),edges);
    for j = 1:nbins
        tmp = max( radius(i,p_shocked{i}(bin==j)));
        if ~isempty(tmp)
            if tmp > shock_radius(i-1,j)
                shock_radius(i,j) = tmp;
            else
                shock_radius(i,j) = shock_radius(i-1,j);
            end
        else
            shock_radius(i,j) = shock_radius(i-1,j);
        end
    end
    shock_radius(i,:) = smooth( shock_radius(i,:), 3 );
end
for j = 1:nbins
    shock_radius(:,j) = smooth( shock_radius(:,j), 20 );
end
% toc

end