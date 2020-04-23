function dd = delta_ab(attack1,defense1,attack2,defense2)
%   Compute the Rue-Salvesen adjustment.
%   from Rue & Salvesen (1997), section 2.2.

dd = (attack1 + defense1 - attack2 - defense2) ./ 2 ;
end

