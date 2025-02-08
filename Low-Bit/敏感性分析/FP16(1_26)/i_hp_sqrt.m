function [a_half,result] = i_hp_sqrt(a)

     % 输出两个参数或一个参数
    if nargout == 2
        if ~isreal(a)
            [a_half,result] = c_hp_sqrt(a);
        else
            [a_half,result] = r_hp_sqrt(a);
        end
    else
        if ~isreal(a)
            [a_half] = c_hp_sqrt(a);
        else
            [a_half] = r_hp_sqrt(a);
        end
    end
end