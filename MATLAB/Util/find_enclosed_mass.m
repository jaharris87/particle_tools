function [ M_r ] = find_enclosed_mass( radius, M_inner, M_tracer )
% Calcuate enclosed mass for all particle assuming spherical mass distribution

    tic
    
    % Initialize array
    M_r = zeros(size(radius));
    
    % Sort the radii for each time step
    sort_rad = sort( radius, 2, 'ascend' );

    for i = 1:size(radius,1)

        % Count the number of particles interior to each particle
        [~,p_inside] = histc( radius(i,:), sort_rad(i,:) );
        
        % Calculate the total mass interior to each particle
        M_r(i,:) = p_inside*M_tracer + M_inner;
        
        % Alternate methods listed below (much slower)
        
%         M_r(i,:) = M_inner;
%         for j = 1:nj
%             for k = 1:nk
%                 if radius(i,k) <= radius(i,j)
%                     M_r(i,j) = M_r(i,j) + M_tracer;
%                 end
%             end
%         end
        
%         for j = 1:nj
%             mask = radius(i,:) <= radius(i,j);
%             M_r(i,j) = nnz(mask)*M_tracer;
%         end
%         M_r(i,:) = M_r(i,:) + M_inner;
        
    end
    
    toc

end
