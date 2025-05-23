function floor_friction_g1_partial_dfr_dfn_custom_gd = floor_friction_g1_partial_dfr_dfn_func_custom_gd(in1)
%FLOOR_FRICTION_G1_PARTIAL_DFR_DFN_FUNC_CUSTOM_GD
%    FLOOR_FRICTION_G1_PARTIAL_DFR_DFN_CUSTOM_GD = FLOOR_FRICTION_G1_PARTIAL_DFR_DFN_FUNC_CUSTOM_GD(IN1)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    19-Dec-2024 20:57:51

dt = in1(:,11);
fn1 = in1(:,7);
fn2 = in1(:,8);
fn3 = in1(:,9);
mu = in1(:,10);
n_gd1 = in1(:,13);
n_gd2 = in1(:,14);
n_gd3 = in1(:,15);
x1s_x = in1(:,1);
x1s_y = in1(:,2);
x1s_z = in1(:,3);
x1s_x0 = in1(:,4);
x1s_y0 = in1(:,5);
x1s_z0 = in1(:,6);
t2 = abs(fn1);
t3 = abs(fn2);
t4 = abs(fn3);
t5 = sign(fn1);
t6 = sign(fn2);
t7 = sign(fn3);
t11 = 1.0./dt;
t12 = -x1s_x0;
t13 = -x1s_y0;
t14 = -x1s_z0;
t8 = t2.^2;
t9 = t3.^2;
t10 = t4.^2;
t15 = t12+x1s_x;
t16 = t13+x1s_y;
t17 = t14+x1s_z;
t18 = t11.*t15;
t19 = t11.*t16;
t20 = t11.*t17;
t24 = t8+t9+t10;
t21 = n_gd1.*t18;
t22 = n_gd2.*t19;
t23 = n_gd3.*t20;
t25 = 1.0./sqrt(t24);
t26 = t21+t22+t23;
t27 = n_gd1.*t26;
t28 = n_gd2.*t26;
t29 = n_gd3.*t26;
t30 = -t27;
t31 = -t28;
t32 = -t29;
t33 = t18+t30;
t34 = t19+t31;
t35 = t20+t32;
t36 = abs(t33);
t37 = abs(t34);
t38 = abs(t35);
t39 = t36.^2;
t40 = t37.^2;
t41 = t38.^2;
t42 = t39+t40+t41;
t43 = 1.0./sqrt(t42);
floor_friction_g1_partial_dfr_dfn_custom_gd = reshape([-mu.*t2.*t5.*t25.*t33.*t43,-mu.*t2.*t5.*t25.*t34.*t43,-mu.*t2.*t5.*t25.*t35.*t43,-mu.*t3.*t6.*t25.*t33.*t43,-mu.*t3.*t6.*t25.*t34.*t43,-mu.*t3.*t6.*t25.*t35.*t43,-mu.*t4.*t7.*t25.*t33.*t43,-mu.*t4.*t7.*t25.*t34.*t43,-mu.*t4.*t7.*t25.*t35.*t43],[3,3]);
end
