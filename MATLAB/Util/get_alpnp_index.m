function [ inuc_alpha, inuc_alpnp ] = get_alpnp_index( aa, zz )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

inuc_nn   = find(aa==1  & zz==0);
inuc_pp   = find(aa==1  & zz==1);
inuc_he4  = find(aa==4  & zz==2);
inuc_c12  = find(aa==12 & zz==6);
inuc_o16  = find(aa==16 & zz==8);
inuc_ne20 = find(aa==20 & zz==10);
inuc_mg24 = find(aa==24 & zz==12);
inuc_si28 = find(aa==28 & zz==14);
inuc_s32  = find(aa==32 & zz==16);
inuc_ar36 = find(aa==36 & zz==18);
inuc_ca40 = find(aa==40 & zz==20);
inuc_fe52 = find(aa==52 & zz==26);
inuc_cr48 = find(aa==48 & zz==24);
inuc_ti44 = find(aa==44 & zz==22);
inuc_ni56 = find(aa==56 & zz==28);
inuc_zn60 = find(aa==60 & zz==30);

inuc_fe56 = find(aa==56 & zz==26);

inuc_alpha = [inuc_he4,inuc_c12,inuc_o16,inuc_ne20,inuc_mg24,inuc_si28,inuc_s32,inuc_ar36,inuc_ca40,inuc_ti44,inuc_cr48,inuc_fe52,inuc_ni56,inuc_zn60];
inuc_alpnp = [inuc_nn,inuc_pp,inuc_he4,inuc_c12,inuc_o16,inuc_ne20,inuc_mg24,inuc_si28,inuc_s32,inuc_ar36,inuc_ca40,inuc_ti44,inuc_cr48,inuc_fe52,inuc_fe56,inuc_ni56,inuc_zn60];

end

