function [ unui, dunui ] = nu_egrid( unumn, unumx, nnugpmx )
%% Calculate CHIMERA neutrino energy grid at infinity

rngpm1 = nnugpmx - 1;
runul = log(unumx/unumn)/rngpm1; 
runu = exp(runul);
runuh = sqrt(runu);
unui(1:nnugpmx) = 0.0;
unubi(1:nnugpmx) = 0.0;
unui(1) = unumn;
unubi(1) = unumn/runuh;
for k = 2:nnugpmx
    unui(k) = runu * unui(k-1);
    unubi(k) = runu * unubi(k-1);
end
unubi(1) = 0;
unubi(nnugpmx+1) = runu*unubi(nnugpmx);
unui(1:nnugpmx) = sqrt(unubi(1:nnugpmx).^2 + unubi(1:nnugpmx).*unubi(2:nnugpmx+1) + unubi(2:nnugpmx+1).^2)/sqrt(3);
dunui = unubi(2:nnugpmx+1)-unui(1:nnugpmx);

end

