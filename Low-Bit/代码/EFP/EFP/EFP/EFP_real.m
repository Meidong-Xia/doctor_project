function [result_fp8Binary, result_config] = EFP_real(fp8Binary1,config1)
    [result_fp8Binary] = arrayfun(@real, fp8Binary1, config1,'UniformOutput',false);
    
    result_config = config1;
end

function result = real(fp8Binary1,config1) 
    config1 = config1{1};
    fp8Binary1 = fp8Binary1{1};
    fp8Binary1(1+config1(1)+config1(2)+1:end) = dec2bin(0,1+config1(1)+config1(2));
    result = fp8Binary1;
end