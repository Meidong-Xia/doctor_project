function [a_half,result] = i_hp_sqrt(a)
    if ~isreal(a)
        [a_half,result] = c_hp_sqrt(a);
    else
        [a_half,result] = r_hp_sqrt(a);
    end
end