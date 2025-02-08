function result=nx(v)
result=i_hp_matrix_mul(complex(v'),complex(v));
result=real(result);
result=i_hp_sqrt(result);
result=real(result);
end