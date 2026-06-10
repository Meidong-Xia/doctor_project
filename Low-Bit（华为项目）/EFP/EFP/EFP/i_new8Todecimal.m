function decimalValue = i_new8Todecimal(fp8Binary,config, base, fraction_tables)
    % fp8Binary = '0000001101101100';
    % config = [2,5,4,2,5,-2];
    % base = 10;
    % config = config{1};
    % fp8Binary = fp8Binary{1};
    fp8Binary_real = fp8Binary(1:1+config(1)+config(2));
    fp8Binary_imag = fp8Binary(1+config(1)+config(2)+1:end);
    dec_real = new8Todecimal(fp8Binary_real,config(1:3), base, fraction_tables);
    dec_imag = new8Todecimal(fp8Binary_imag,config(6:8), base, fraction_tables);
    decimalValue = dec_real + 1i*dec_imag;
end