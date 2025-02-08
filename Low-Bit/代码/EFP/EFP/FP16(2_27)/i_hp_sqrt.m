function [a_half,result] = i_hp_sqrt(a)

     % 输出两个参数或一个参数
    if nargout == 2
            a_half = hp_turn(double(a));
            result = hp_turn(double(sqrt(a_half)));
    else
            a_half = hp_turn(double(a));
            a_half = hp_turn(double(sqrt(a_half)));
    end
end