function [result_fp8Binary, config1] = EFP_neg(fp8Binary1,config1)
    [result_fp8Binary] = arrayfun(@neg, fp8Binary1, config1,'UniformOutput',false);
end

function result = neg(fp8Binary1,config1) 
    config1 = config1{1};
    fp8Binary1 = fp8Binary1{1};
    if bin2dec(fp8Binary1(1+config1(1)+config1(2)+1:end)) ~= 0
        fp8Binary1(1+config1(1)+config1(2)+1) = dec2bin(1-bin2dec(fp8Binary1(1+config1(1)+config1(2)+1)));
    end
    if bin2dec(fp8Binary1(1:1+config1(1)+config1(2))) ~= 0
        fp8Binary1(1) = dec2bin(1-bin2dec(fp8Binary1(1)));
    end
    result = fp8Binary1;
end