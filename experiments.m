% experiments.m   This script runs the experiments in [Sec. 5, 1].
%
% Requirements:
%   CPFloat (https://github.com/north-numerical-computing/cpfloat/).
%
% References:
%   [1] T. Mary and M. Mikaitis.
%       Error Analysis of Matrix Multiplication with Narrow Range
%       Floating-Point Arithmetic. hal-04671474. Aug. 2024.

clear all; rng(1);

input_format = 'fp8-e4m3';
accum_format = 'binary16';
p = 1;
matmul_test;
p = 2;

matmul_test;
p = 3;
matmul_test;

input_format = 'fp8-e5m2';
accum_format = 'binary16';
p = 1;
matmul_test;
p = 2;
matmul_test;
p = 3;
matmul_test;

input_format = 'fp8-e4m3';
accum_format = 'binary32';
p = 1;
matmul_test;
p = 2;
matmul_test;
p = 3;
matmul_test;

input_format = 'fp8-e5m2';
accum_format = 'binary32';
p = 1;
matmul_test;
p = 2;
matmul_test;
p = 3;
matmul_test;

input_format = 'binary16';
accum_format = 'binary32';
p = 1;
matmul_test;
p = 2;
matmul_test;
p = 3;
matmul_test;
