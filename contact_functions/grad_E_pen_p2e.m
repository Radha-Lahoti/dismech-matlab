function grad_E1_p2e = grad_E_pen_p2e(in1)
%grad_E_pen_p2e
%    grad_E1_p2e = grad_E_pen_p2e(IN1)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    04-Sep-2024 17:55:23

Delta = in1(:,13);
h = in1(:,14);
x2e1 = in1(:,10);
x2e2 = in1(:,11);
x2e3 = in1(:,12);
x1s1 = in1(:,1);
x1s2 = in1(:,2);
x1s3 = in1(:,3);
x2s1 = in1(:,7);
x2s2 = in1(:,8);
x2s3 = in1(:,9);
t2 = Delta.*2.0;
t3 = h.*4.0;
t5 = -x1s1;
t6 = -x1s2;
t7 = -x1s3;
t8 = -x2s1;
t9 = -x2s2;
t10 = -x2s3;
t4 = -t3;
t11 = t5+x2e1;
t12 = t6+x2e2;
t13 = t7+x2e3;
t14 = t8+x2e1;
t15 = t9+x2e2;
t16 = t10+x2e3;
t17 = t8+x1s1;
t18 = t9+x1s2;
t19 = t10+x1s3;
t20 = abs(t14);
t21 = abs(t15);
t22 = abs(t16);
t23 = sign(t14);
t24 = sign(t15);
t25 = sign(t16);
t26 = t2+t4;
t30 = t14.*t18;
t31 = t15.*t17;
t32 = t14.*t19;
t33 = t16.*t17;
t34 = t15.*t19;
t35 = t16.*t18;
t27 = t20.^2;
t28 = t21.^2;
t29 = t22.^2;
t36 = -t31;
t37 = -t33;
t38 = -t35;
t39 = t27+t28+t29;
t40 = t30+t36;
t41 = t32+t37;
t42 = t34+t38;
t43 = abs(t40);
t44 = abs(t41);
t45 = abs(t42);
t46 = sign(t40);
t47 = sign(t41);
t48 = sign(t42);
t52 = 1.0./sqrt(t39);
t49 = t43.^2;
t50 = t44.^2;
t51 = t45.^2;
t53 = t52.^3;
t54 = t49+t50+t51;
t55 = sqrt(t54);
t56 = 1.0./t55;
t57 = t20.*t23.*t53.*t55;
t58 = t21.*t24.*t53.*t55;
t59 = t22.*t25.*t53.*t55;
grad_E1_p2e = [t26.*t52.*t56.*(t15.*t43.*t46.*2.0+t16.*t44.*t47.*2.0).*(-1.0./2.0),(t26.*t52.*t56.*(t14.*t43.*t46.*2.0-t16.*t45.*t48.*2.0))./2.0,(t26.*t52.*t56.*(t14.*t44.*t47.*2.0+t15.*t45.*t48.*2.0))./2.0,0.0,0.0,0.0,t26.*(t57+(t52.*t56.*(t12.*t43.*t46.*2.0+t13.*t44.*t47.*2.0))./2.0),t26.*(t58-(t52.*t56.*(t11.*t43.*t46.*2.0-t13.*t45.*t48.*2.0))./2.0),t26.*(t59-(t52.*t56.*(t11.*t44.*t47.*2.0+t12.*t45.*t48.*2.0))./2.0),-t26.*(t57-(t52.*t56.*(t18.*t43.*t46.*2.0+t19.*t44.*t47.*2.0))./2.0),-t26.*(t58+(t52.*t56.*(t17.*t43.*t46.*2.0-t19.*t45.*t48.*2.0))./2.0),-t26.*(t59+(t52.*t56.*(t17.*t44.*t47.*2.0+t18.*t45.*t48.*2.0))./2.0)];
