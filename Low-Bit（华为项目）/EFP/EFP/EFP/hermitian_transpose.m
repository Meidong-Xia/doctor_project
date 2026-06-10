function [result_fp8Binary, result_config] = hermitian_transpose(fp8Binary1,config1)
    [result] = arrayfun(@Hermitian, fp8Binary1, config1,'UniformOutput',false);
    result_fp8Binary = result';
    result_config = config1';
end

function result = Hermitian(fp8Binary1,config1) 
    config1 = config1{1};
    fp8Binary1 = fp8Binary1{1};
    if bin2dec(fp8Binary1(1+config1(1)+config1(2)+1:end)) ~= 0
        fp8Binary1(1+config1(1)+config1(2)+1) = dec2bin(1-bin2dec(fp8Binary1(1+config1(1)+config1(2)+1)));
    end
    result = fp8Binary1;
end