function [fp8Binary,result_config] = decToEFP_auto(decimalValue,m_bit, base, fraction_tables)
     config_1 = [2 m_bit 2 0 0 2 m_bit 2 0 0];
     [sizem,sizen] = size(decimalValue);	
     config = repmat({config_1},sizem,sizen);
     [fp8Binary,result_config] = decToEFP(decimalValue,config,base,fraction_tables);
end