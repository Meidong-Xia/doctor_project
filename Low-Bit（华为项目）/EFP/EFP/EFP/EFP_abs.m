function [result_fp8Binary, result_config] = EFP_abs(fp8Binary1,config1)
    [result] = arrayfun(@abs, fp8Binary1, config1,'UniformOutput',false);
    result_fp8Binary = result;
    result_config = config1;
end

function result = abs(fp8Binary1,config1) 
    config1 = config1{1};
    fp8Binary1 = fp8Binary1{1};
    if ~strcmp(fp8Binary1(1:1+config1(1)+config1(2)),['1',dec2bin(0,config1(1)+config1(2))])
        fp8Binary1(1) = '0';
    end
    if ~strcmp(fp8Binary1(1+config1(1)+config1(2)+1:1+config1(1)+config1(2)+1+config1(1)+config1(2)),['1',dec2bin(0,config1(1)+config1(2))])
        fp8Binary1(1+config1(1)+config1(2)+1) = '0';
    end
    result = fp8Binary1;
end