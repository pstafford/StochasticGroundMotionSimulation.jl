
# Definition of Boore & Thompson (2012) constant values to be used within subsequent functions
const m_ii_bt15 = [ 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0 ]
const r_jj_bt15 = [ 2.00, 3.17, 5.02, 7.96, 12.62, 20.00, 31.70, 50.24, 79.62, 126.20, 200.01, 317.00, 502.41, 796.26, 1262.00 ]
const num_m_ii_bt15 = 13
const num_r_jj_bt15 = 15


const coefs_wna_bt15 = [ 2.00	2.00	9.38E-01	-1.09E-02	2.00E+00	1.00E+00	1.86E+00	2.00E+00	1.07E+00	1.03E+00	1.00E+00;
2.50	2.00	9.06E-01	-1.38E-02	2.00E+00	1.00E+00	1.67E+00	2.00E+00	1.07E+00	1.05E+00	1.04E+00;
3.00	2.00	9.55E-01	-1.74E-02	2.00E+00	1.00E+00	1.17E+00	2.00E+00	1.28E+00	1.04E+00	1.03E+00;
3.50	2.00	9.55E-01	-2.27E-02	2.00E+00	1.00E+00	9.35E-01	2.00E+00	1.29E+00	1.04E+00	1.03E+00;
4.00	2.00	9.85E-01	-3.71E-02	2.00E+00	1.00E+00	6.72E-01	2.00E+00	1.37E+00	1.03E+00	1.02E+00;
4.50	2.00	1.01E+00	-7.99E-02	2.00E+00	1.00E+00	5.18E-01	2.00E+00	1.39E+00	1.04E+00	1.02E+00;
5.00	2.00	1.04E+00	-1.10E-01	2.00E+00	1.00E+00	4.31E-01	2.00E+00	1.34E+00	1.03E+00	1.01E+00;
5.50	2.00	1.03E+00	-1.28E-01	2.00E+00	1.00E+00	3.20E-01	2.00E+00	1.34E+00	1.05E+00	1.01E+00;
6.00	2.00	9.54E-01	-7.01E-02	2.00E+00	1.00E+00	2.23E-01	2.00E+00	1.16E+00	1.05E+00	1.00E+00;
6.50	2.00	9.52E-01	-6.39E-02	2.00E+00	1.00E+00	1.86E-01	2.00E+00	1.16E+00	1.05E+00	9.96E-01;
7.00	2.00	8.48E-01	2.92E-02	2.00E+00	1.00E+00	1.15E-01	2.00E+00	1.15E+00	1.07E+00	9.88E-01;
7.50	2.00	8.49E-01	3.80E-02	2.00E+00	1.00E+00	1.02E-01	2.00E+00	1.15E+00	1.07E+00	9.65E-01;
8.00	2.00	8.30E-01	3.78E-02	2.00E+00	1.00E+00	8.68E-02	2.00E+00	1.05E+00	1.08E+00	9.50E-01;
2.00	3.17	9.17E-01	-9.48E-03	2.00E+00	1.00E+00	3.35E+00	2.00E+00	1.01E+00	1.05E+00	1.02E+00;
2.50	3.17	9.17E-01	-2.05E-02	2.00E+00	1.00E+00	3.08E+00	2.00E+00	1.01E+00	1.05E+00	1.04E+00;
3.00	3.17	9.37E-01	-1.12E-02	2.00E+00	1.00E+00	1.89E+00	2.00E+00	1.19E+00	1.04E+00	1.03E+00;
3.50	3.17	9.36E-01	-1.44E-02	2.00E+00	1.00E+00	1.48E+00	2.00E+00	1.19E+00	1.04E+00	1.03E+00;
4.00	3.17	9.62E-01	-3.49E-02	2.00E+00	1.00E+00	9.96E-01	2.00E+00	1.27E+00	1.04E+00	1.02E+00;
4.50	3.17	9.66E-01	-4.16E-02	2.00E+00	1.00E+00	6.78E-01	2.00E+00	1.30E+00	1.04E+00	1.03E+00;
5.00	3.17	1.01E+00	-9.19E-02	2.00E+00	1.00E+00	5.24E-01	2.00E+00	1.34E+00	1.04E+00	1.02E+00;
5.50	3.17	1.02E+00	-1.15E-01	2.00E+00	1.00E+00	3.74E-01	2.00E+00	1.29E+00	1.05E+00	1.01E+00;
6.00	3.17	9.80E-01	-8.70E-02	2.00E+00	1.00E+00	2.67E-01	2.00E+00	1.15E+00	1.05E+00	9.85E-01;
6.50	3.17	9.45E-01	-5.79E-02	2.00E+00	1.00E+00	1.90E-01	2.00E+00	1.17E+00	1.06E+00	9.84E-01;
7.00	3.17	8.36E-01	4.35E-02	2.00E+00	1.00E+00	1.15E-01	2.00E+00	1.16E+00	1.07E+00	9.83E-01;
7.50	3.17	8.45E-01	3.80E-02	2.00E+00	1.00E+00	1.02E-01	2.00E+00	1.13E+00	1.08E+00	9.61E-01;
8.00	3.17	8.28E-01	3.80E-02	2.00E+00	1.00E+00	7.62E-02	2.00E+00	1.06E+00	1.08E+00	9.41E-01;
2.00	5.02	9.13E-01	-4.79E-03	2.00E+00	1.00E+00	5.72E+00	2.00E+00	9.90E-01	1.05E+00	1.04E+00;
2.50	5.02	9.02E-01	-3.31E-03	2.00E+00	1.00E+00	5.19E+00	2.00E+00	1.02E+00	1.05E+00	1.04E+00;
3.00	5.02	9.03E-01	-2.62E-03	2.00E+00	1.00E+00	3.71E+00	2.00E+00	1.04E+00	1.05E+00	1.05E+00;
3.50	5.02	9.09E-01	-2.44E-03	2.00E+00	1.00E+00	2.44E+00	2.00E+00	1.12E+00	1.05E+00	1.05E+00;
4.00	5.02	9.39E-01	-2.64E-02	2.00E+00	1.00E+00	1.49E+00	2.00E+00	1.21E+00	1.05E+00	1.05E+00;
4.50	5.02	9.36E-01	-2.14E-02	2.00E+00	1.00E+00	1.05E+00	2.00E+00	1.19E+00	1.04E+00	1.03E+00;
5.00	5.02	9.76E-01	-6.53E-02	2.00E+00	1.00E+00	6.74E-01	2.00E+00	1.25E+00	1.05E+00	1.03E+00;
5.50	5.02	9.95E-01	-9.77E-02	2.00E+00	1.00E+00	4.58E-01	2.00E+00	1.23E+00	1.06E+00	1.02E+00;
6.00	5.02	9.96E-01	-1.02E-01	2.00E+00	1.00E+00	3.34E-01	2.00E+00	1.21E+00	1.06E+00	1.01E+00;
6.50	5.02	8.88E-01	-1.34E-02	2.00E+00	1.00E+00	1.90E-01	2.00E+00	1.15E+00	1.07E+00	1.00E+00;
7.00	5.02	8.24E-01	5.87E-02	2.00E+00	1.00E+00	1.20E-01	2.00E+00	1.16E+00	1.07E+00	9.91E-01;
7.50	5.02	8.41E-01	3.74E-02	2.00E+00	1.00E+00	9.77E-02	2.00E+00	1.08E+00	1.08E+00	9.59E-01;
8.00	5.02	8.04E-01	6.20E-02	2.00E+00	1.00E+00	7.63E-02	2.00E+00	1.06E+00	1.08E+00	9.48E-01;
2.00	7.96	9.10E-01	-1.24E-02	2.00E+00	1.00E+00	1.07E+01	2.00E+00	9.03E-01	1.05E+00	1.03E+00;
2.50	7.96	9.03E-01	6.10E-04	2.00E+00	1.00E+00	7.60E+00	2.00E+00	1.01E+00	1.05E+00	1.05E+00;
3.00	7.96	8.93E-01	-4.72E-03	2.00E+00	1.00E+00	6.14E+00	2.00E+00	1.00E+00	1.06E+00	1.06E+00;
3.50	7.96	8.93E-01	-4.55E-03	2.00E+00	1.00E+00	3.78E+00	2.00E+00	1.00E+00	1.05E+00	1.05E+00;
4.00	7.96	9.15E-01	-2.36E-02	2.00E+00	1.00E+00	2.52E+00	2.00E+00	1.07E+00	1.06E+00	1.05E+00;
4.50	7.96	9.29E-01	-2.47E-02	2.00E+00	1.00E+00	1.51E+00	2.00E+00	1.16E+00	1.06E+00	1.05E+00;
5.00	7.96	9.71E-01	-6.76E-02	2.00E+00	1.00E+00	9.17E-01	2.00E+00	1.21E+00	1.06E+00	1.04E+00;
5.50	7.96	9.69E-01	-6.62E-02	2.00E+00	1.00E+00	6.06E-01	2.00E+00	1.20E+00	1.06E+00	1.03E+00;
6.00	7.96	9.49E-01	-6.62E-02	2.00E+00	1.00E+00	3.53E-01	2.00E+00	1.15E+00	1.06E+00	1.01E+00;
6.50	7.96	8.31E-01	4.55E-02	2.00E+00	1.00E+00	1.75E-01	2.00E+00	1.15E+00	1.07E+00	1.02E+00;
7.00	7.96	8.31E-01	4.58E-02	2.00E+00	1.00E+00	1.38E-01	2.00E+00	1.13E+00	1.07E+00	9.90E-01;
7.50	7.96	8.31E-01	4.14E-02	2.00E+00	1.00E+00	1.12E-01	2.00E+00	1.12E+00	1.08E+00	9.75E-01;
8.00	7.96	8.04E-01	6.82E-02	2.00E+00	1.00E+00	7.01E-02	2.00E+00	1.05E+00	1.08E+00	9.48E-01;
2.00	12.62	9.01E-01	-2.38E-04	2.00E+00	1.00E+00	1.39E+01	2.00E+00	9.03E-01	1.05E+00	1.05E+00;
2.50	12.62	8.78E-01	-2.16E-03	2.00E+00	1.00E+00	1.20E+01	2.00E+00	8.87E-01	1.06E+00	1.06E+00;
3.00	12.62	8.67E-01	-2.53E-03	2.00E+00	1.00E+00	9.21E+00	2.00E+00	8.98E-01	1.07E+00	1.06E+00;
3.50	12.62	8.88E-01	-4.74E-03	2.00E+00	1.00E+00	4.96E+00	2.00E+00	1.00E+00	1.07E+00	1.06E+00;
4.00	12.62	8.89E-01	-1.24E-02	2.00E+00	1.00E+00	3.28E+00	2.00E+00	1.02E+00	1.06E+00	1.06E+00;
4.50	12.62	8.90E-01	-1.27E-02	2.00E+00	1.00E+00	1.71E+00	2.00E+00	1.02E+00	1.06E+00	1.05E+00;
5.00	12.62	9.15E-01	-2.13E-02	2.00E+00	1.00E+00	9.17E-01	2.00E+00	1.14E+00	1.06E+00	1.05E+00;
5.50	12.62	9.69E-01	-6.63E-02	2.00E+00	1.00E+00	6.29E-01	2.00E+00	1.20E+00	1.06E+00	1.03E+00;
6.00	12.62	9.99E-01	-1.15E-01	2.00E+00	1.00E+00	5.22E-01	2.00E+00	1.16E+00	1.06E+00	1.02E+00;
6.50	12.62	8.21E-01	5.22E-02	2.00E+00	1.00E+00	1.75E-01	2.00E+00	1.08E+00	1.07E+00	1.00E+00;
7.00	12.62	8.42E-01	3.63E-02	2.00E+00	1.00E+00	1.49E-01	2.00E+00	1.10E+00	1.07E+00	9.84E-01;
7.50	12.62	8.18E-01	5.02E-02	2.00E+00	1.00E+00	1.12E-01	2.00E+00	1.09E+00	1.08E+00	9.69E-01;
8.00	12.62	8.02E-01	6.36E-02	2.00E+00	1.00E+00	7.30E-02	2.00E+00	1.07E+00	1.08E+00	9.47E-01;
2.00	20.00	8.67E-01	1.11E-03	2.00E+00	1.00E+00	2.31E+01	2.00E+00	8.44E-01	1.07E+00	1.05E+00;
2.50	20.00	8.79E-01	-2.62E-03	2.00E+00	1.00E+00	1.54E+01	2.00E+00	9.00E-01	1.07E+00	1.06E+00;
3.00	20.00	8.78E-01	-2.33E-03	2.00E+00	1.00E+00	1.24E+01	2.00E+00	9.00E-01	1.06E+00	1.05E+00;
3.50	20.00	8.78E-01	-9.54E-03	2.00E+00	1.00E+00	7.24E+00	2.00E+00	9.47E-01	1.07E+00	1.06E+00;
4.00	20.00	8.75E-01	-8.07E-03	2.00E+00	1.00E+00	4.40E+00	2.00E+00	9.51E-01	1.06E+00	1.06E+00;
4.50	20.00	8.82E-01	-1.08E-02	2.00E+00	1.00E+00	2.28E+00	2.00E+00	1.00E+00	1.07E+00	1.06E+00;
5.00	20.00	8.90E-01	-1.53E-02	2.00E+00	1.00E+00	1.35E+00	2.00E+00	1.01E+00	1.06E+00	1.04E+00;
5.50	20.00	8.86E-01	-1.56E-02	2.00E+00	1.00E+00	7.89E-01	2.00E+00	1.01E+00	1.06E+00	1.03E+00;
6.00	20.00	9.69E-01	-8.21E-02	2.00E+00	1.00E+00	5.22E-01	2.00E+00	1.15E+00	1.07E+00	1.03E+00;
6.50	20.00	9.68E-01	-8.34E-02	2.00E+00	1.00E+00	3.54E-01	2.00E+00	1.17E+00	1.06E+00	1.01E+00;
7.00	20.00	8.30E-01	3.98E-02	2.00E+00	1.00E+00	1.57E-01	2.00E+00	1.11E+00	1.08E+00	1.00E+00;
7.50	20.00	8.36E-01	4.00E-02	2.00E+00	1.00E+00	1.19E-01	2.00E+00	1.10E+00	1.07E+00	9.68E-01;
8.00	20.00	8.01E-01	6.36E-02	2.00E+00	1.00E+00	7.10E-02	2.00E+00	1.05E+00	1.08E+00	9.49E-01;
2.00	31.70	8.77E-01	-2.85E-03	2.00E+00	1.00E+00	2.84E+01	2.00E+00	8.36E-01	1.07E+00	1.06E+00;
2.50	31.70	8.72E-01	-2.76E-03	2.00E+00	1.00E+00	2.74E+01	2.00E+00	8.37E-01	1.07E+00	1.06E+00;
3.00	31.70	8.76E-01	-2.14E-03	2.00E+00	1.00E+00	1.88E+01	2.00E+00	9.14E-01	1.07E+00	1.06E+00;
3.50	31.70	8.55E-01	-9.99E-03	2.00E+00	1.00E+00	1.05E+01	2.00E+00	8.84E-01	1.09E+00	1.07E+00;
4.00	31.70	8.76E-01	-1.01E-02	2.00E+00	1.00E+00	6.42E+00	2.00E+00	9.42E-01	1.07E+00	1.06E+00;
4.50	31.70	8.68E-01	-7.45E-03	2.00E+00	1.00E+00	3.65E+00	2.00E+00	9.52E-01	1.07E+00	1.06E+00;
5.00	31.70	8.84E-01	-2.10E-02	2.00E+00	1.00E+00	1.97E+00	2.00E+00	9.67E-01	1.07E+00	1.05E+00;
5.50	31.70	8.91E-01	-2.11E-02	2.00E+00	1.00E+00	1.06E+00	2.00E+00	9.86E-01	1.06E+00	1.03E+00;
6.00	31.70	9.70E-01	-8.86E-02	2.00E+00	1.00E+00	6.97E-01	2.00E+00	1.13E+00	1.07E+00	1.03E+00;
6.50	31.70	9.66E-01	-8.88E-02	2.00E+00	1.00E+00	4.23E-01	2.00E+00	1.12E+00	1.07E+00	1.00E+00;
7.00	31.70	7.93E-01	7.72E-02	2.00E+00	1.00E+00	1.57E-01	2.00E+00	1.07E+00	1.07E+00	1.00E+00;
7.50	31.70	8.07E-01	6.50E-02	2.00E+00	1.00E+00	1.26E-01	2.00E+00	1.08E+00	1.08E+00	9.72E-01;
8.00	31.70	7.98E-01	6.55E-02	2.00E+00	1.00E+00	8.50E-02	2.00E+00	1.07E+00	1.08E+00	9.51E-01;
2.00	50.24	8.70E-01	-2.75E-03	2.00E+00	1.00E+00	3.34E+01	2.00E+00	8.31E-01	1.08E+00	1.06E+00;
2.50	50.24	8.63E-01	-1.52E-02	2.00E+00	1.00E+00	2.83E+01	2.00E+00	8.21E-01	1.08E+00	1.07E+00;
3.00	50.24	8.76E-01	-1.96E-03	2.00E+00	1.00E+00	1.80E+01	2.00E+00	9.11E-01	1.07E+00	1.07E+00;
3.50	50.24	8.54E-01	2.04E-03	2.00E+00	1.00E+00	1.37E+01	2.00E+00	8.84E-01	1.08E+00	1.07E+00;
4.00	50.24	8.79E-01	-1.08E-02	2.00E+00	1.00E+00	6.94E+00	2.00E+00	9.35E-01	1.07E+00	1.07E+00;
4.50	50.24	8.75E-01	-1.25E-02	2.00E+00	1.00E+00	4.72E+00	2.00E+00	9.33E-01	1.07E+00	1.06E+00;
5.00	50.24	8.79E-01	-1.76E-02	2.00E+00	1.00E+00	2.58E+00	2.00E+00	9.66E-01	1.07E+00	1.06E+00;
5.50	50.24	8.76E-01	-7.74E-03	2.00E+00	1.00E+00	1.30E+00	2.00E+00	9.89E-01	1.07E+00	1.05E+00;
6.00	50.24	8.97E-01	-3.42E-02	2.00E+00	1.00E+00	7.04E-01	2.00E+00	9.98E-01	1.07E+00	1.02E+00;
6.50	50.24	9.07E-01	-4.26E-02	2.00E+00	1.00E+00	4.49E-01	2.00E+00	1.04E+00	1.07E+00	1.01E+00;
7.00	50.24	9.07E-01	-4.26E-02	2.00E+00	1.00E+00	3.05E-01	2.00E+00	1.05E+00	1.07E+00	1.01E+00;
7.50	50.24	8.06E-01	6.49E-02	2.00E+00	1.00E+00	1.32E-01	2.00E+00	1.08E+00	1.07E+00	9.82E-01;
8.00	50.24	8.00E-01	6.46E-02	2.00E+00	1.00E+00	9.83E-02	2.00E+00	1.08E+00	1.08E+00	9.62E-01;
2.00	79.62	8.77E-01	-1.54E-02	2.00E+00	1.00E+00	2.13E+01	2.00E+00	8.24E-01	1.08E+00	1.06E+00;
2.50	79.62	8.74E-01	-1.34E-02	2.00E+00	1.00E+00	1.89E+01	2.00E+00	8.41E-01	1.07E+00	1.07E+00;
3.00	79.62	8.58E-01	-8.54E-03	2.00E+00	1.00E+00	1.51E+01	2.00E+00	8.63E-01	1.08E+00	1.08E+00;
3.50	79.62	8.70E-01	-1.69E-02	2.00E+00	1.00E+00	1.17E+01	2.00E+00	8.59E-01	1.07E+00	1.07E+00;
4.00	79.62	8.69E-01	-1.02E-02	2.00E+00	1.00E+00	6.79E+00	2.00E+00	9.17E-01	1.07E+00	1.06E+00;
4.50	79.62	8.64E-01	-9.08E-03	2.00E+00	1.00E+00	4.01E+00	2.00E+00	9.22E-01	1.08E+00	1.07E+00;
5.00	79.62	8.80E-01	-2.34E-02	2.00E+00	1.00E+00	2.25E+00	2.00E+00	9.66E-01	1.07E+00	1.06E+00;
5.50	79.62	8.77E-01	-1.97E-02	2.00E+00	1.00E+00	1.26E+00	2.00E+00	9.69E-01	1.07E+00	1.05E+00;
6.00	79.62	8.77E-01	-1.15E-02	2.00E+00	1.00E+00	7.06E-01	2.00E+00	1.05E+00	1.07E+00	1.04E+00;
6.50	79.62	8.93E-01	-2.12E-02	2.00E+00	1.00E+00	4.34E-01	2.00E+00	1.08E+00	1.07E+00	1.01E+00;
7.00	79.62	8.91E-01	-2.06E-02	2.00E+00	1.00E+00	3.06E-01	2.00E+00	1.11E+00	1.07E+00	1.00E+00;
7.50	79.62	8.87E-01	-1.92E-02	2.00E+00	1.00E+00	1.99E-01	2.00E+00	1.12E+00	1.07E+00	9.88E-01;
8.00	79.62	7.90E-01	6.96E-02	2.00E+00	1.00E+00	1.06E-01	2.00E+00	1.09E+00	1.08E+00	9.64E-01;
2.00	126.20	8.67E-01	-1.54E-02	2.00E+00	1.00E+00	1.06E+01	2.00E+00	8.24E-01	1.08E+00	1.07E+00;
2.50	126.20	8.70E-01	-1.99E-02	2.00E+00	1.00E+00	1.11E+01	2.00E+00	8.25E-01	1.08E+00	1.07E+00;
3.00	126.20	8.68E-01	-2.25E-02	2.00E+00	1.00E+00	9.66E+00	2.00E+00	8.31E-01	1.08E+00	1.07E+00;
3.50	126.20	8.69E-01	-2.90E-02	2.00E+00	1.00E+00	9.08E+00	2.00E+00	8.38E-01	1.08E+00	1.08E+00;
4.00	126.20	8.66E-01	-1.41E-02	2.00E+00	1.00E+00	5.26E+00	2.00E+00	9.05E-01	1.08E+00	1.08E+00;
4.50	126.20	8.67E-01	-1.24E-02	2.00E+00	1.00E+00	3.47E+00	2.00E+00	9.22E-01	1.07E+00	1.06E+00;
5.00	126.20	8.83E-01	-3.09E-02	2.00E+00	1.00E+00	2.31E+00	2.00E+00	9.33E-01	1.07E+00	1.05E+00;
5.50	126.20	8.87E-01	-3.49E-02	2.00E+00	1.00E+00	1.33E+00	2.00E+00	9.58E-01	1.07E+00	1.04E+00;
6.00	126.20	8.85E-01	-4.67E-03	2.00E+00	1.00E+00	6.42E-01	2.00E+00	1.11E+00	1.06E+00	1.04E+00;
6.50	126.20	8.83E-01	-7.45E-03	2.00E+00	1.00E+00	4.14E-01	2.00E+00	1.11E+00	1.06E+00	1.02E+00;
7.00	126.20	9.14E-01	-3.72E-02	2.00E+00	1.00E+00	3.05E-01	2.00E+00	1.13E+00	1.06E+00	9.97E-01;
7.50	126.20	8.54E-01	1.79E-02	2.00E+00	1.00E+00	1.91E-01	2.00E+00	1.12E+00	1.07E+00	9.80E-01;
8.00	126.20	8.06E-01	6.96E-02	2.00E+00	1.00E+00	1.13E-01	2.00E+00	1.17E+00	1.07E+00	9.62E-01;
2.00	200.01	8.56E-01	-1.56E-02	2.00E+00	1.00E+00	1.03E+01	2.00E+00	8.08E-01	1.09E+00	1.07E+00;
2.50	200.01	8.54E-01	-2.40E-02	2.00E+00	1.00E+00	1.13E+01	2.00E+00	7.70E-01	1.09E+00	1.09E+00;
3.00	200.01	8.57E-01	-2.40E-02	2.00E+00	1.00E+00	1.16E+01	2.00E+00	7.70E-01	1.09E+00	1.07E+00;
3.50	200.01	8.67E-01	-2.90E-02	2.00E+00	1.00E+00	9.08E+00	2.00E+00	8.30E-01	1.09E+00	1.08E+00;
4.00	200.01	8.69E-01	-2.93E-02	2.00E+00	1.00E+00	7.53E+00	2.00E+00	8.32E-01	1.08E+00	1.07E+00;
4.50	200.01	8.73E-01	-3.34E-02	2.00E+00	1.00E+00	5.56E+00	2.00E+00	8.32E-01	1.08E+00	1.07E+00;
5.00	200.01	8.81E-01	-3.48E-02	2.00E+00	1.00E+00	3.46E+00	2.00E+00	8.72E-01	1.08E+00	1.06E+00;
5.50	200.01	9.15E-01	-6.55E-02	2.00E+00	1.00E+00	2.27E+00	2.00E+00	9.25E-01	1.07E+00	1.05E+00;
6.00	200.01	9.16E-01	-5.45E-02	2.00E+00	1.00E+00	1.37E+00	2.00E+00	9.66E-01	1.07E+00	1.05E+00;
6.50	200.01	8.12E-01	6.67E-02	2.00E+00	1.00E+00	5.05E-01	2.00E+00	1.06E+00	1.06E+00	1.03E+00;
7.00	200.01	8.10E-01	7.00E-02	2.00E+00	1.00E+00	3.05E-01	2.00E+00	1.11E+00	1.06E+00	1.02E+00;
7.50	200.01	8.40E-01	3.43E-02	2.00E+00	1.00E+00	2.54E-01	2.00E+00	1.17E+00	1.07E+00	1.01E+00;
8.00	200.01	8.39E-01	3.23E-02	2.00E+00	1.00E+00	1.87E-01	2.00E+00	1.17E+00	1.07E+00	9.74E-01;
2.00	317.00	9.43E-01	-1.02E-01	2.00E+00	1.00E+00	1.21E+01	2.00E+00	7.87E-01	1.08E+00	1.08E+00;
2.50	317.00	8.98E-01	-6.21E-02	2.00E+00	1.00E+00	1.02E+01	2.00E+00	7.73E-01	1.09E+00	1.08E+00;
3.00	317.00	8.87E-01	-5.91E-02	2.00E+00	1.00E+00	1.05E+01	2.00E+00	7.62E-01	1.09E+00	1.09E+00;
3.50	317.00	8.88E-01	-5.91E-02	2.00E+00	1.00E+00	9.81E+00	2.00E+00	7.62E-01	1.09E+00	1.08E+00;
4.00	317.00	8.79E-01	-2.97E-02	2.00E+00	1.00E+00	7.52E+00	2.00E+00	8.23E-01	1.08E+00	1.08E+00;
4.50	317.00	8.67E-01	-3.08E-02	2.00E+00	1.00E+00	5.96E+00	2.00E+00	8.07E-01	1.08E+00	1.07E+00;
5.00	317.00	8.69E-01	-3.36E-02	2.00E+00	1.00E+00	5.00E+00	2.00E+00	8.24E-01	1.08E+00	1.08E+00;
5.50	317.00	9.22E-01	-6.81E-02	2.00E+00	1.00E+00	3.09E+00	2.00E+00	9.07E-01	1.08E+00	1.07E+00;
6.00	317.00	9.16E-01	-5.09E-02	2.00E+00	1.00E+00	2.01E+00	2.00E+00	9.57E-01	1.07E+00	1.05E+00;
6.50	317.00	9.05E-01	-5.09E-02	2.00E+00	1.00E+00	1.30E+00	2.00E+00	9.57E-01	1.07E+00	1.05E+00;
7.00	317.00	7.84E-01	9.62E-02	2.00E+00	1.00E+00	4.18E-01	2.00E+00	1.10E+00	1.06E+00	1.04E+00;
7.50	317.00	8.33E-01	3.42E-02	2.00E+00	1.00E+00	2.99E-01	2.00E+00	1.12E+00	1.07E+00	1.02E+00;
8.00	317.00	8.48E-01	3.23E-02	2.00E+00	1.00E+00	1.87E-01	2.00E+00	1.17E+00	1.07E+00	9.92E-01;
2.00	502.41	9.28E-01	-9.26E-02	2.00E+00	1.00E+00	6.71E+00	2.00E+00	7.45E-01	1.08E+00	1.08E+00;
2.50	502.41	9.23E-01	-9.30E-02	2.00E+00	1.00E+00	6.52E+00	2.00E+00	7.44E-01	1.09E+00	1.08E+00;
3.00	502.41	9.27E-01	-9.51E-02	2.00E+00	1.00E+00	6.52E+00	2.00E+00	7.57E-01	1.09E+00	1.07E+00;
3.50	502.41	9.27E-01	-9.25E-02	2.00E+00	1.00E+00	6.01E+00	2.00E+00	7.57E-01	1.09E+00	1.08E+00;
4.00	502.41	9.18E-01	-9.25E-02	2.00E+00	1.00E+00	6.01E+00	2.00E+00	7.57E-01	1.09E+00	1.08E+00;
4.50	502.41	8.63E-01	-1.46E-02	2.00E+00	1.00E+00	4.55E+00	2.00E+00	8.12E-01	1.07E+00	1.07E+00;
5.00	502.41	8.64E-01	-3.31E-02	2.00E+00	1.00E+00	3.53E+00	2.00E+00	8.25E-01	1.09E+00	1.08E+00;
5.50	502.41	9.22E-01	-6.59E-02	2.00E+00	1.00E+00	3.32E+00	2.00E+00	9.05E-01	1.07E+00	1.07E+00;
6.00	502.41	9.16E-01	-5.09E-02	2.00E+00	1.00E+00	2.01E+00	2.00E+00	9.57E-01	1.07E+00	1.06E+00;
6.50	502.41	8.75E-01	9.81E-03	2.00E+00	1.00E+00	1.25E+00	2.00E+00	9.97E-01	1.06E+00	1.05E+00;
7.00	502.41	7.79E-01	9.10E-02	2.00E+00	1.00E+00	5.57E-01	2.00E+00	1.07E+00	1.07E+00	1.05E+00;
7.50	502.41	7.96E-01	9.16E-02	2.00E+00	1.00E+00	2.64E-01	2.00E+00	1.16E+00	1.06E+00	1.04E+00;
8.00	502.41	8.53E-01	3.23E-02	2.00E+00	1.00E+00	1.92E-01	2.00E+00	1.17E+00	1.06E+00	1.00E+00;
2.00	796.26	9.17E-01	-8.91E-02	2.00E+00	1.00E+00	3.60E+00	2.00E+00	7.48E-01	1.09E+00	1.07E+00;
2.50	796.26	9.17E-01	-7.64E-02	2.00E+00	1.00E+00	3.60E+00	2.00E+00	7.79E-01	1.08E+00	1.08E+00;
3.00	796.26	9.14E-01	-7.40E-02	2.00E+00	1.00E+00	3.33E+00	2.00E+00	7.79E-01	1.08E+00	1.08E+00;
3.50	796.26	9.16E-01	-7.70E-02	2.00E+00	1.00E+00	3.36E+00	2.00E+00	7.86E-01	1.08E+00	1.08E+00;
4.00	796.26	9.18E-01	-7.88E-02	2.00E+00	1.00E+00	3.51E+00	2.00E+00	7.89E-01	1.08E+00	1.08E+00;
4.50	796.26	8.72E-01	-3.31E-02	2.00E+00	1.00E+00	2.60E+00	2.00E+00	8.12E-01	1.08E+00	1.08E+00;
5.00	796.26	8.72E-01	-3.09E-02	2.00E+00	1.00E+00	2.55E+00	2.00E+00	8.06E-01	1.08E+00	1.08E+00;
5.50	796.26	9.24E-01	-8.56E-02	2.00E+00	1.00E+00	2.12E+00	2.00E+00	8.66E-01	1.09E+00	1.08E+00;
6.00	796.26	9.02E-01	-5.29E-02	2.00E+00	1.00E+00	1.89E+00	2.00E+00	8.98E-01	1.08E+00	1.07E+00;
6.50	796.26	8.50E-01	9.81E-03	2.00E+00	1.00E+00	1.17E+00	2.00E+00	9.79E-01	1.07E+00	1.06E+00;
7.00	796.26	8.49E-01	8.98E-03	2.00E+00	1.00E+00	6.98E-01	2.00E+00	9.83E-01	1.07E+00	1.05E+00;
7.50	796.26	8.61E-01	2.80E-02	2.00E+00	1.00E+00	3.00E-01	2.00E+00	1.14E+00	1.06E+00	1.04E+00;
8.00	796.26	8.65E-01	2.74E-02	2.00E+00	1.00E+00	8.62E-02	2.00E+00	1.23E+00	1.06E+00	1.02E+00;
2.00	1262.00	7.65E-01	8.02E-02	2.00E+00	1.00E+00	8.90E-01	2.00E+00	8.18E-01	1.08E+00	1.07E+00;
2.50	1262.00	7.63E-01	8.02E-02	2.00E+00	1.00E+00	7.95E-01	2.00E+00	8.14E-01	1.08E+00	1.07E+00;
3.00	1262.00	7.64E-01	8.05E-02	2.00E+00	1.00E+00	7.48E-01	2.00E+00	8.18E-01	1.08E+00	1.07E+00;
3.50	1262.00	7.50E-01	8.04E-02	2.00E+00	1.00E+00	7.40E-01	2.00E+00	8.03E-01	1.09E+00	1.08E+00;
4.00	1262.00	9.87E-01	-1.40E-01	2.00E+00	1.00E+00	1.68E+00	2.00E+00	8.39E-01	1.08E+00	1.07E+00;
4.50	1262.00	9.82E-01	-1.45E-01	2.00E+00	1.00E+00	1.35E+00	2.00E+00	8.20E-01	1.09E+00	1.07E+00;
5.00	1262.00	9.93E-01	-1.40E-01	2.00E+00	1.00E+00	1.68E+00	2.00E+00	8.50E-01	1.07E+00	1.08E+00;
5.50	1262.00	9.72E-01	-1.21E-01	2.00E+00	1.00E+00	1.20E+00	2.00E+00	8.71E-01	1.08E+00	1.07E+00;
6.00	1262.00	9.00E-01	-5.20E-02	2.00E+00	1.00E+00	1.10E+00	2.00E+00	9.00E-01	1.08E+00	1.07E+00;
6.50	1262.00	8.51E-01	-8.69E-03	2.00E+00	1.00E+00	2.15E-01	2.00E+00	9.31E-01	1.08E+00	1.07E+00;
7.00	1262.00	8.74E-01	-4.90E-03	2.00E+00	1.00E+00	2.06E-03	2.00E+00	1.03E+00	1.07E+00	1.06E+00;
7.50	1262.00	9.11E-01	-2.45E-02	2.00E+00	1.00E+00	5.54E-06	2.00E+00	1.17E+00	1.06E+00	1.05E+00;
8.00	1262.00	9.16E-01	-2.36E-02	2.00E+00	1.00E+00	1.51E-09	2.00E+00	1.21E+00	1.05E+00	1.04E+00 ]




const coefs_ena_bt15 = [ 2.00	2.00	9.29E-01	-1.27E-02	2.00E+00	1.00E+00	1.41E+00	1.90E+00	1.15E+00	1.05E+00	1.05E+00;
2.50	2.00	9.20E-01	-2.34E-02	2.00E+00	1.00E+00	1.19E+00	1.89E+00	1.11E+00	1.05E+00	1.05E+00;
3.00	2.00	9.38E-01	-2.74E-02	2.00E+00	1.00E+00	8.36E-01	1.91E+00	1.19E+00	1.05E+00	1.05E+00;
3.50	2.00	9.80E-01	-5.66E-02	2.00E+00	1.00E+00	5.75E-01	1.89E+00	1.35E+00	1.05E+00	1.03E+00;
4.00	2.00	1.03E+00	-1.01E-01	2.00E+00	1.00E+00	5.04E-01	1.85E+00	1.34E+00	1.06E+00	1.02E+00;
4.50	2.00	1.04E+00	-1.14E-01	2.00E+00	1.00E+00	4.30E-01	1.79E+00	1.37E+00	1.05E+00	1.01E+00;
5.00	2.00	1.06E+00	-1.41E-01	2.00E+00	1.00E+00	3.77E-01	1.76E+00	1.43E+00	1.05E+00	9.81E-01;
5.50	2.00	1.06E+00	-1.41E-01	2.00E+00	1.00E+00	3.48E-01	1.72E+00	1.43E+00	1.06E+00	9.55E-01;
6.00	2.00	9.48E-01	-5.23E-02	2.00E+00	1.00E+00	2.23E-01	1.77E+00	1.29E+00	1.07E+00	9.54E-01;
6.50	2.00	9.48E-01	-5.22E-02	2.00E+00	1.00E+00	2.06E-01	1.77E+00	1.28E+00	1.07E+00	9.36E-01;
7.00	2.00	8.56E-01	2.30E-02	2.00E+00	1.00E+00	1.34E-01	1.88E+00	1.22E+00	1.08E+00	9.23E-01;
7.50	2.00	8.52E-01	2.46E-02	2.00E+00	1.00E+00	1.24E-01	1.88E+00	1.22E+00	1.08E+00	9.08E-01;
8.00	2.00	8.45E-01	2.16E-02	2.00E+00	1.00E+00	1.13E-01	1.82E+00	1.13E+00	1.08E+00	9.05E-01;
2.00	3.17	9.11E-01	-3.46E-02	2.00E+00	1.00E+00	2.70E+00	2.06E+00	9.58E-01	1.05E+00	1.06E+00;
2.50	3.17	9.10E-01	-3.94E-02	2.00E+00	1.00E+00	1.92E+00	2.03E+00	9.58E-01	1.07E+00	1.06E+00;
3.00	3.17	9.50E-01	-6.92E-02	2.00E+00	1.00E+00	1.21E+00	2.10E+00	1.07E+00	1.05E+00	1.06E+00;
3.50	3.17	9.73E-01	-4.57E-02	2.00E+00	1.00E+00	7.38E-01	2.05E+00	1.28E+00	1.06E+00	1.04E+00;
4.00	3.17	9.96E-01	-1.01E-01	2.00E+00	1.00E+00	5.96E-01	1.99E+00	1.19E+00	1.06E+00	1.03E+00;
4.50	3.17	9.99E-01	-9.17E-02	2.00E+00	1.00E+00	4.77E-01	1.90E+00	1.28E+00	1.05E+00	1.02E+00;
5.00	3.17	1.05E+00	-1.29E-01	2.00E+00	1.00E+00	4.23E-01	1.84E+00	1.41E+00	1.06E+00	9.81E-01;
5.50	3.17	9.97E-01	-8.61E-02	2.00E+00	1.00E+00	3.50E-01	1.73E+00	1.35E+00	1.07E+00	9.76E-01;
6.00	3.17	9.51E-01	-6.70E-02	2.00E+00	1.00E+00	2.86E-01	1.76E+00	1.30E+00	1.06E+00	9.45E-01;
6.50	3.17	9.04E-01	-1.56E-02	2.00E+00	1.00E+00	2.05E-01	1.78E+00	1.27E+00	1.07E+00	9.32E-01;
7.00	3.17	8.44E-01	3.24E-02	2.00E+00	1.00E+00	1.34E-01	1.88E+00	1.22E+00	1.08E+00	9.28E-01;
7.50	3.17	8.44E-01	3.26E-02	2.00E+00	1.00E+00	1.16E-01	1.89E+00	1.22E+00	1.07E+00	8.98E-01;
8.00	3.17	8.28E-01	3.91E-02	2.00E+00	1.00E+00	9.04E-02	1.87E+00	1.12E+00	1.09E+00	8.98E-01;
2.00	5.02	9.03E-01	-2.31E-02	2.00E+00	1.00E+00	4.70E+00	2.03E+00	9.16E-01	1.06E+00	1.06E+00;
2.50	5.02	9.09E-01	-1.77E-02	2.00E+00	1.00E+00	2.89E+00	2.06E+00	1.01E+00	1.07E+00	1.06E+00;
3.00	5.02	9.08E-01	-3.10E-02	2.00E+00	1.00E+00	2.07E+00	2.04E+00	1.00E+00	1.06E+00	1.05E+00;
3.50	5.02	9.40E-01	-4.10E-02	2.00E+00	1.00E+00	1.21E+00	2.08E+00	1.11E+00	1.07E+00	1.06E+00;
4.00	5.02	9.56E-01	-6.43E-02	2.00E+00	1.00E+00	7.86E-01	2.04E+00	1.15E+00	1.06E+00	1.04E+00;
4.50	5.02	9.53E-01	-6.41E-02	2.00E+00	1.00E+00	6.75E-01	1.92E+00	1.15E+00	1.07E+00	1.04E+00;
5.00	5.02	9.77E-01	-7.32E-02	2.00E+00	1.00E+00	5.21E-01	1.81E+00	1.28E+00	1.06E+00	1.01E+00;
5.50	5.02	1.02E+00	-1.24E-01	2.00E+00	1.00E+00	4.02E-01	1.81E+00	1.31E+00	1.06E+00	9.69E-01;
6.00	5.02	9.36E-01	-4.62E-02	2.00E+00	1.00E+00	2.86E-01	1.78E+00	1.28E+00	1.07E+00	9.61E-01;
6.50	5.02	8.92E-01	-1.15E-02	2.00E+00	1.00E+00	2.10E-01	1.80E+00	1.18E+00	1.07E+00	9.36E-01;
7.00	5.02	8.21E-01	5.64E-02	2.00E+00	1.00E+00	1.34E-01	1.89E+00	1.19E+00	1.07E+00	9.40E-01;
7.50	5.02	8.21E-01	4.97E-02	2.00E+00	1.00E+00	1.11E-01	1.92E+00	1.17E+00	1.08E+00	9.23E-01;
8.00	5.02	8.09E-01	5.00E-02	2.00E+00	1.00E+00	9.04E-02	1.90E+00	1.11E+00	1.08E+00	8.94E-01;
2.00	7.96	8.92E-01	-7.58E-03	2.00E+00	1.00E+00	8.62E+00	2.04E+00	9.17E-01	1.07E+00	1.07E+00;
2.50	7.96	8.96E-01	-1.84E-02	2.00E+00	1.00E+00	5.22E+00	2.06E+00	9.44E-01	1.07E+00	1.05E+00;
3.00	7.96	8.97E-01	-3.22E-02	2.00E+00	1.00E+00	3.75E+00	2.07E+00	9.16E-01	1.08E+00	1.06E+00;
3.50	7.96	8.92E-01	-3.79E-02	2.00E+00	1.00E+00	2.45E+00	2.05E+00	9.15E-01	1.08E+00	1.06E+00;
4.00	7.96	9.00E-01	-2.92E-02	2.00E+00	1.00E+00	1.44E+00	1.95E+00	9.95E-01	1.06E+00	1.04E+00;
4.50	7.96	9.52E-01	-6.98E-02	2.00E+00	1.00E+00	7.96E-01	2.06E+00	1.11E+00	1.06E+00	1.04E+00;
5.00	7.96	9.51E-01	-7.59E-02	2.00E+00	1.00E+00	5.52E-01	1.97E+00	1.11E+00	1.07E+00	1.02E+00;
5.50	7.96	9.58E-01	-7.26E-02	2.00E+00	1.00E+00	4.54E-01	1.88E+00	1.20E+00	1.07E+00	9.80E-01;
6.00	7.96	9.61E-01	-6.65E-02	2.00E+00	1.00E+00	3.30E-01	1.84E+00	1.18E+00	1.07E+00	9.59E-01;
6.50	7.96	8.89E-01	-1.09E-02	2.00E+00	1.00E+00	2.11E-01	1.85E+00	1.17E+00	1.07E+00	9.38E-01;
7.00	7.96	8.21E-01	5.64E-02	2.00E+00	1.00E+00	1.34E-01	1.89E+00	1.13E+00	1.08E+00	9.44E-01;
7.50	7.96	8.21E-01	5.26E-02	2.00E+00	1.00E+00	1.10E-01	1.95E+00	1.16E+00	1.07E+00	9.13E-01;
8.00	7.96	8.09E-01	5.02E-02	2.00E+00	1.00E+00	9.48E-02	1.89E+00	1.06E+00	1.08E+00	9.04E-01;
2.00	12.62	8.80E-01	-1.09E-02	2.00E+00	1.00E+00	1.44E+01	1.97E+00	8.61E-01	1.08E+00	1.07E+00;
2.50	12.62	8.77E-01	-1.29E-02	2.00E+00	1.00E+00	1.09E+01	1.98E+00	8.50E-01	1.08E+00	1.07E+00;
3.00	12.62	8.94E-01	-2.13E-02	2.00E+00	1.00E+00	6.61E+00	2.11E+00	9.07E-01	1.07E+00	1.07E+00;
3.50	12.62	8.92E-01	-3.71E-02	2.00E+00	1.00E+00	4.01E+00	2.16E+00	9.08E-01	1.08E+00	1.07E+00;
4.00	12.62	9.01E-01	-3.09E-02	2.00E+00	1.00E+00	2.09E+00	2.06E+00	9.64E-01	1.08E+00	1.06E+00;
4.50	12.62	9.45E-01	-7.00E-02	2.00E+00	1.00E+00	1.35E+00	2.15E+00	1.03E+00	1.07E+00	1.04E+00;
5.00	12.62	9.77E-01	-1.02E-01	2.00E+00	1.00E+00	9.01E-01	2.07E+00	1.06E+00	1.07E+00	1.01E+00;
5.50	12.62	1.01E+00	-1.29E-01	2.00E+00	1.00E+00	5.91E-01	2.00E+00	1.11E+00	1.07E+00	9.96E-01;
6.00	12.62	9.42E-01	-6.54E-02	2.00E+00	1.00E+00	3.90E-01	1.88E+00	1.14E+00	1.07E+00	9.74E-01;
6.50	12.62	9.40E-01	-7.12E-02	2.00E+00	1.00E+00	3.00E-01	1.87E+00	1.12E+00	1.07E+00	9.46E-01;
7.00	12.62	7.88E-01	8.40E-02	2.00E+00	1.00E+00	1.30E-01	1.94E+00	1.13E+00	1.08E+00	9.46E-01;
7.50	12.62	7.88E-01	8.09E-02	2.00E+00	1.00E+00	8.08E-02	2.11E+00	1.15E+00	1.07E+00	9.23E-01;
8.00	12.62	7.92E-01	7.05E-02	2.00E+00	1.00E+00	4.83E-02	2.28E+00	1.07E+00	1.08E+00	9.09E-01;
2.00	20.00	8.54E-01	-3.91E-03	2.00E+00	1.00E+00	6.77E+01	1.83E+00	7.47E-01	1.09E+00	1.08E+00;
2.50	20.00	8.48E-01	-3.92E-03	2.00E+00	1.00E+00	4.46E+01	1.85E+00	7.70E-01	1.09E+00	1.08E+00;
3.00	20.00	8.43E-01	-1.93E-03	2.00E+00	1.00E+00	2.99E+01	1.88E+00	7.54E-01	1.09E+00	1.07E+00;
3.50	20.00	8.38E-01	-5.17E-03	2.00E+00	1.00E+00	1.71E+01	1.87E+00	7.58E-01	1.08E+00	1.08E+00;
4.00	20.00	8.28E-01	2.97E-03	2.00E+00	1.00E+00	8.98E+00	1.85E+00	7.70E-01	1.08E+00	1.06E+00;
4.50	20.00	9.01E-01	-6.28E-02	2.00E+00	1.00E+00	5.92E+00	2.30E+00	8.47E-01	1.08E+00	1.06E+00;
5.00	20.00	8.99E-01	-5.38E-02	2.00E+00	1.00E+00	2.99E+00	2.14E+00	8.59E-01	1.08E+00	1.04E+00;
5.50	20.00	8.99E-01	-4.81E-02	2.00E+00	1.00E+00	1.45E+00	2.03E+00	9.10E-01	1.08E+00	1.02E+00;
6.00	20.00	8.89E-01	-4.06E-02	2.00E+00	1.00E+00	7.91E-01	1.94E+00	9.21E-01	1.07E+00	1.01E+00;
6.50	20.00	8.99E-01	-3.99E-02	2.00E+00	1.00E+00	5.11E-01	1.92E+00	1.02E+00	1.08E+00	9.93E-01;
7.00	20.00	9.00E-01	-4.37E-02	2.00E+00	1.00E+00	3.20E-01	1.92E+00	1.05E+00	1.07E+00	9.73E-01;
7.50	20.00	9.00E-01	-4.30E-02	2.00E+00	1.00E+00	2.34E-01	1.91E+00	1.08E+00	1.08E+00	9.42E-01;
8.00	20.00	7.78E-01	7.59E-02	2.00E+00	1.00E+00	4.83E-02	2.56E+00	1.06E+00	1.08E+00	9.28E-01;
2.00	31.70	8.45E-01	2.02E-03	2.00E+00	1.00E+00	2.00E+02	1.76E+00	7.32E-01	1.09E+00	1.09E+00;
2.50	31.70	8.45E-01	-5.58E-03	2.00E+00	1.00E+00	1.08E+02	1.70E+00	7.30E-01	1.10E+00	1.09E+00;
3.00	31.70	8.41E-01	-1.08E-02	2.00E+00	1.00E+00	8.64E+01	1.77E+00	6.97E-01	1.09E+00	1.09E+00;
3.50	31.70	8.33E-01	-5.48E-03	2.00E+00	1.00E+00	4.94E+01	1.81E+00	7.33E-01	1.09E+00	1.08E+00;
4.00	31.70	8.34E-01	-5.48E-03	2.00E+00	1.00E+00	2.98E+01	1.80E+00	7.21E-01	1.09E+00	1.08E+00;
4.50	31.70	8.34E-01	-5.48E-03	2.00E+00	1.00E+00	1.37E+01	1.86E+00	7.39E-01	1.09E+00	1.07E+00;
5.00	31.70	9.04E-01	-6.29E-02	2.00E+00	1.00E+00	8.88E+00	2.32E+00	8.18E-01	1.08E+00	1.06E+00;
5.50	31.70	8.92E-01	-5.15E-02	2.00E+00	1.00E+00	3.99E+00	2.16E+00	8.12E-01	1.08E+00	1.05E+00;
6.00	31.70	8.92E-01	-4.99E-02	2.00E+00	1.00E+00	2.05E+00	2.11E+00	8.49E-01	1.08E+00	1.02E+00;
6.50	31.70	8.97E-01	-5.50E-02	2.00E+00	1.00E+00	1.02E+00	2.04E+00	8.86E-01	1.07E+00	1.00E+00;
7.00	31.70	9.12E-01	-5.74E-02	2.00E+00	1.00E+00	5.88E-01	2.00E+00	9.85E-01	1.08E+00	9.96E-01;
7.50	31.70	8.98E-01	-4.83E-02	2.00E+00	1.00E+00	4.22E-01	1.88E+00	9.88E-01	1.08E+00	9.61E-01;
8.00	31.70	7.14E-01	1.34E-01	2.00E+00	1.00E+00	4.68E-02	2.85E+00	9.96E-01	1.08E+00	9.41E-01;
2.00	50.24	8.23E-01	9.71E-03	2.00E+00	1.00E+00	1.99E+02	1.45E+00	6.53E-01	1.10E+00	1.09E+00;
2.50	50.24	8.17E-01	4.57E-03	2.00E+00	1.00E+00	1.92E+02	1.52E+00	6.37E-01	1.10E+00	1.09E+00;
3.00	50.24	8.16E-01	4.75E-03	2.00E+00	1.00E+00	1.21E+02	1.53E+00	6.53E-01	1.09E+00	1.09E+00;
3.50	50.24	8.31E-01	-1.80E-03	2.00E+00	1.00E+00	6.94E+01	1.64E+00	7.11E-01	1.09E+00	1.08E+00;
4.00	50.24	8.28E-01	-1.80E-03	2.00E+00	1.00E+00	4.47E+01	1.70E+00	7.11E-01	1.09E+00	1.08E+00;
4.50	50.24	8.40E-01	-8.55E-03	2.00E+00	1.00E+00	2.47E+01	1.86E+00	7.43E-01	1.08E+00	1.08E+00;
5.00	50.24	8.41E-01	-1.05E-02	2.00E+00	1.00E+00	1.42E+01	1.87E+00	7.48E-01	1.09E+00	1.07E+00;
5.50	50.24	8.81E-01	-4.36E-02	2.00E+00	1.00E+00	7.52E+00	2.11E+00	8.01E-01	1.09E+00	1.06E+00;
6.00	50.24	9.08E-01	-6.64E-02	2.00E+00	1.00E+00	3.80E+00	2.20E+00	8.44E-01	1.08E+00	1.06E+00;
6.50	50.24	9.03E-01	-6.34E-02	2.00E+00	1.00E+00	1.81E+00	2.17E+00	8.72E-01	1.09E+00	1.04E+00;
7.00	50.24	8.99E-01	-6.06E-02	2.00E+00	1.00E+00	9.85E-01	2.07E+00	9.32E-01	1.08E+00	1.01E+00;
7.50	50.24	9.06E-01	-5.68E-02	2.00E+00	1.00E+00	5.25E-01	2.06E+00	9.78E-01	1.09E+00	9.83E-01;
8.00	50.24	9.02E-01	-5.65E-02	2.00E+00	1.00E+00	4.38E-01	1.85E+00	9.78E-01	1.09E+00	9.61E-01;
2.00	79.62	8.13E-01	9.64E-03	2.00E+00	1.00E+00	1.62E+02	1.47E+00	6.28E-01	1.10E+00	1.09E+00;
2.50	79.62	8.15E-01	4.60E-03	2.00E+00	1.00E+00	1.40E+02	1.52E+00	6.36E-01	1.10E+00	1.09E+00;
3.00	79.62	8.15E-01	3.78E-03	2.00E+00	1.00E+00	9.38E+01	1.57E+00	6.57E-01	1.09E+00	1.09E+00;
3.50	79.62	8.31E-01	-1.75E-03	2.00E+00	1.00E+00	5.71E+01	1.64E+00	7.08E-01	1.09E+00	1.08E+00;
4.00	79.62	8.31E-01	-1.58E-03	2.00E+00	1.00E+00	3.60E+01	1.73E+00	7.29E-01	1.09E+00	1.08E+00;
4.50	79.62	8.47E-01	-1.34E-02	2.00E+00	1.00E+00	2.20E+01	1.91E+00	7.56E-01	1.08E+00	1.08E+00;
5.00	79.62	8.40E-01	-1.03E-02	2.00E+00	1.00E+00	1.29E+01	1.85E+00	7.45E-01	1.08E+00	1.07E+00;
5.50	79.62	8.80E-01	-4.34E-02	2.00E+00	1.00E+00	6.95E+00	2.11E+00	8.01E-01	1.09E+00	1.06E+00;
6.00	79.62	9.08E-01	-6.79E-02	2.00E+00	1.00E+00	3.68E+00	2.18E+00	8.43E-01	1.08E+00	1.06E+00;
6.50	79.62	9.02E-01	-6.39E-02	2.00E+00	1.00E+00	1.76E+00	2.15E+00	8.73E-01	1.09E+00	1.04E+00;
7.00	79.62	8.99E-01	-6.06E-02	2.00E+00	1.00E+00	9.65E-01	2.07E+00	9.32E-01	1.08E+00	1.02E+00;
7.50	79.62	9.06E-01	-5.68E-02	2.00E+00	1.00E+00	5.25E-01	2.05E+00	9.78E-01	1.08E+00	9.87E-01;
8.00	79.62	9.06E-01	-6.15E-02	2.00E+00	1.00E+00	4.35E-01	1.87E+00	9.92E-01	1.09E+00	9.67E-01;
2.00	126.20	8.07E-01	1.13E-02	2.00E+00	1.00E+00	9.27E+01	1.48E+00	6.28E-01	1.10E+00	1.10E+00;
2.50	126.20	7.99E-01	1.14E-02	2.00E+00	1.00E+00	7.95E+01	1.51E+00	6.28E-01	1.10E+00	1.10E+00;
3.00	126.20	8.15E-01	1.60E-03	2.00E+00	1.00E+00	6.07E+01	1.56E+00	6.60E-01	1.10E+00	1.09E+00;
3.50	126.20	8.21E-01	4.21E-03	2.00E+00	1.00E+00	4.81E+01	1.68E+00	6.96E-01	1.09E+00	1.09E+00;
4.00	126.20	8.21E-01	4.21E-03	2.00E+00	1.00E+00	2.97E+01	1.68E+00	7.06E-01	1.09E+00	1.08E+00;
4.50	126.20	8.39E-01	-8.12E-03	2.00E+00	1.00E+00	1.85E+01	1.84E+00	7.39E-01	1.09E+00	1.08E+00;
5.00	126.20	8.36E-01	-9.62E-03	2.00E+00	1.00E+00	1.02E+01	1.85E+00	7.46E-01	1.08E+00	1.07E+00;
5.50	126.20	8.37E-01	2.48E-03	2.00E+00	1.00E+00	4.92E+00	1.91E+00	8.11E-01	1.08E+00	1.06E+00;
6.00	126.20	9.19E-01	-7.39E-02	2.00E+00	1.00E+00	3.31E+00	2.14E+00	8.67E-01	1.08E+00	1.05E+00;
6.50	126.20	9.04E-01	-6.36E-02	2.00E+00	1.00E+00	1.55E+00	2.12E+00	8.89E-01	1.08E+00	1.04E+00;
7.00	126.20	9.10E-01	-6.30E-02	2.00E+00	1.00E+00	9.99E-01	2.01E+00	9.37E-01	1.08E+00	1.01E+00;
7.50	126.20	9.10E-01	-5.98E-02	2.00E+00	1.00E+00	4.80E-01	2.11E+00	9.96E-01	1.09E+00	9.99E-01;
8.00	126.20	9.11E-01	-6.10E-02	2.00E+00	1.00E+00	3.88E-01	1.98E+00	1.03E+00	1.08E+00	9.65E-01;
2.00	200.01	8.03E-01	1.14E-02	2.00E+00	1.00E+00	5.52E+01	1.54E+00	6.35E-01	1.10E+00	1.10E+00;
2.50	200.01	8.05E-01	1.04E-02	2.00E+00	1.00E+00	5.48E+01	1.56E+00	6.42E-01	1.10E+00	1.10E+00;
3.00	200.01	7.96E-01	9.85E-03	2.00E+00	1.00E+00	4.57E+01	1.61E+00	6.42E-01	1.10E+00	1.09E+00;
3.50	200.01	8.17E-01	4.13E-03	2.00E+00	1.00E+00	3.42E+01	1.67E+00	6.94E-01	1.09E+00	1.09E+00;
4.00	200.01	8.21E-01	3.96E-03	2.00E+00	1.00E+00	2.26E+01	1.69E+00	7.14E-01	1.09E+00	1.09E+00;
4.50	200.01	8.36E-01	-5.16E-03	2.00E+00	1.00E+00	1.55E+01	1.81E+00	7.40E-01	1.08E+00	1.09E+00;
5.00	200.01	8.58E-01	-1.69E-02	2.00E+00	1.00E+00	8.78E+00	1.94E+00	7.98E-01	1.08E+00	1.08E+00;
5.50	200.01	8.57E-01	-2.11E-02	2.00E+00	1.00E+00	5.28E+00	1.92E+00	8.12E-01	1.08E+00	1.08E+00;
6.00	200.01	9.16E-01	-7.26E-02	2.00E+00	1.00E+00	3.31E+00	2.12E+00	8.67E-01	1.07E+00	1.05E+00;
6.50	200.01	9.12E-01	-6.39E-02	2.00E+00	1.00E+00	1.72E+00	2.09E+00	9.05E-01	1.07E+00	1.04E+00;
7.00	200.01	9.13E-01	-5.82E-02	2.00E+00	1.00E+00	9.01E-01	2.07E+00	9.64E-01	1.07E+00	1.01E+00;
7.50	200.01	9.09E-01	-4.69E-02	2.00E+00	1.00E+00	5.33E-01	2.03E+00	1.01E+00	1.07E+00	1.01E+00;
8.00	200.01	9.11E-01	-6.10E-02	2.00E+00	1.00E+00	3.73E-01	1.98E+00	1.06E+00	1.08E+00	9.65E-01;
2.00	317.00	7.96E-01	1.29E-02	2.00E+00	1.00E+00	3.56E+01	1.54E+00	6.16E-01	1.10E+00	1.10E+00;
2.50	317.00	8.04E-01	1.04E-02	2.00E+00	1.00E+00	3.17E+01	1.56E+00	6.42E-01	1.10E+00	1.09E+00;
3.00	317.00	7.98E-01	9.68E-03	2.00E+00	1.00E+00	2.78E+01	1.60E+00	6.45E-01	1.10E+00	1.09E+00;
3.50	317.00	8.17E-01	6.15E-03	2.00E+00	1.00E+00	2.27E+01	1.66E+00	6.91E-01	1.09E+00	1.08E+00;
4.00	317.00	8.07E-01	1.58E-02	2.00E+00	1.00E+00	1.68E+01	1.60E+00	7.06E-01	1.09E+00	1.09E+00;
4.50	317.00	8.36E-01	-2.58E-03	2.00E+00	1.00E+00	1.51E+01	1.77E+00	7.40E-01	1.08E+00	1.08E+00;
5.00	317.00	8.62E-01	-1.75E-02	2.00E+00	1.00E+00	9.30E+00	1.84E+00	7.97E-01	1.08E+00	1.08E+00;
5.50	317.00	8.50E-01	-1.08E-02	2.00E+00	1.00E+00	5.28E+00	1.90E+00	8.12E-01	1.08E+00	1.07E+00;
6.00	317.00	9.11E-01	-6.30E-02	2.00E+00	1.00E+00	3.75E+00	2.16E+00	8.89E-01	1.07E+00	1.06E+00;
6.50	317.00	9.01E-01	-4.47E-02	2.00E+00	1.00E+00	1.85E+00	2.09E+00	9.26E-01	1.07E+00	1.05E+00;
7.00	317.00	8.90E-01	-4.04E-02	2.00E+00	1.00E+00	1.17E+00	2.12E+00	9.31E-01	1.07E+00	1.03E+00;
7.50	317.00	9.09E-01	-4.71E-02	2.00E+00	1.00E+00	6.17E-01	2.27E+00	1.04E+00	1.08E+00	1.02E+00;
8.00	317.00	9.47E-01	-8.14E-02	2.00E+00	1.00E+00	5.15E-01	2.08E+00	1.13E+00	1.08E+00	9.94E-01;
2.00	502.41	7.98E-01	1.55E-02	2.00E+00	1.00E+00	1.93E+01	1.55E+00	6.24E-01	1.10E+00	1.09E+00;
2.50	502.41	8.09E-01	7.96E-03	2.00E+00	1.00E+00	1.93E+01	1.64E+00	6.51E-01	1.09E+00	1.09E+00;
3.00	502.41	7.97E-01	1.72E-02	2.00E+00	1.00E+00	1.70E+01	1.59E+00	6.58E-01	1.10E+00	1.09E+00;
3.50	502.41	7.99E-01	2.10E-02	2.00E+00	1.00E+00	1.63E+01	1.56E+00	6.81E-01	1.09E+00	1.09E+00;
4.00	502.41	7.99E-01	2.48E-02	2.00E+00	1.00E+00	1.17E+01	1.60E+00	7.12E-01	1.09E+00	1.09E+00;
4.50	502.41	8.36E-01	-3.24E-03	2.00E+00	1.00E+00	1.12E+01	1.77E+00	7.52E-01	1.09E+00	1.09E+00;
5.00	502.41	8.35E-01	-4.26E-03	2.00E+00	1.00E+00	8.83E+00	1.78E+00	7.54E-01	1.09E+00	1.09E+00;
5.50	502.41	8.53E-01	-1.36E-02	2.00E+00	1.00E+00	5.88E+00	1.89E+00	8.12E-01	1.08E+00	1.07E+00;
6.00	502.41	9.13E-01	-6.23E-02	2.00E+00	1.00E+00	4.09E+00	2.11E+00	8.79E-01	1.08E+00	1.07E+00;
6.50	502.41	9.03E-01	-4.36E-02	2.00E+00	1.00E+00	2.47E+00	2.11E+00	9.26E-01	1.08E+00	1.06E+00;
7.00	502.41	9.03E-01	-4.38E-02	2.00E+00	1.00E+00	1.52E+00	2.16E+00	9.33E-01	1.07E+00	1.05E+00;
7.50	502.41	1.14E+00	-2.86E-01	2.00E+00	1.00E+00	1.83E+00	2.33E+00	1.02E+00	1.08E+00	1.03E+00;
8.00	502.41	1.16E+00	-2.90E-01	2.00E+00	1.00E+00	1.29E+00	2.17E+00	1.10E+00	1.07E+00	1.02E+00;
2.00	796.26	9.55E-01	-1.31E-01	2.00E+00	1.00E+00	1.76E+01	2.00E+00	6.96E-01	1.10E+00	1.10E+00;
2.50	796.26	7.94E-01	2.00E-02	2.00E+00	1.00E+00	9.44E+00	1.60E+00	6.70E-01	1.10E+00	1.10E+00;
3.00	796.26	8.02E-01	1.36E-02	2.00E+00	1.00E+00	9.54E+00	1.62E+00	6.63E-01	1.10E+00	1.09E+00;
3.50	796.26	8.02E-01	1.34E-02	2.00E+00	1.00E+00	9.54E+00	1.54E+00	6.63E-01	1.10E+00	1.09E+00;
4.00	796.26	7.68E-01	5.02E-02	2.00E+00	1.00E+00	7.31E+00	1.51E+00	6.90E-01	1.10E+00	1.09E+00;
4.50	796.26	8.23E-01	6.32E-03	2.00E+00	1.00E+00	8.01E+00	1.74E+00	7.33E-01	1.10E+00	1.09E+00;
5.00	796.26	8.28E-01	1.78E-04	2.00E+00	1.00E+00	6.91E+00	1.75E+00	7.51E-01	1.09E+00	1.08E+00;
5.50	796.26	8.55E-01	-1.67E-02	2.00E+00	1.00E+00	5.63E+00	1.85E+00	7.99E-01	1.08E+00	1.07E+00;
6.00	796.26	8.52E-01	-1.20E-02	2.00E+00	1.00E+00	4.44E+00	2.00E+00	8.20E-01	1.08E+00	1.08E+00;
6.50	796.26	9.03E-01	-4.83E-02	2.00E+00	1.00E+00	3.06E+00	2.26E+00	9.14E-01	1.08E+00	1.07E+00;
7.00	796.26	9.35E-01	-7.41E-02	2.00E+00	1.00E+00	2.38E+00	2.05E+00	9.55E-01	1.08E+00	1.06E+00;
7.50	796.26	1.15E+00	-2.85E-01	2.00E+00	1.00E+00	2.49E+00	2.45E+00	1.03E+00	1.08E+00	1.05E+00;
8.00	796.26	1.18E+00	-3.10E-01	2.00E+00	1.00E+00	1.39E+00	2.42E+00	1.15E+00	1.07E+00	1.04E+00;
2.00	1262.00	1.06E+00	-2.31E-01	2.00E+00	1.00E+00	1.08E+01	1.97E+00	7.07E-01	1.10E+00	1.09E+00;
2.50	1262.00	1.06E+00	-2.32E-01	2.00E+00	1.00E+00	1.04E+01	1.94E+00	7.08E-01	1.10E+00	1.09E+00;
3.00	1262.00	1.06E+00	-2.32E-01	2.00E+00	1.00E+00	1.02E+01	1.94E+00	7.08E-01	1.09E+00	1.09E+00;
3.50	1262.00	1.06E+00	-2.36E-01	2.00E+00	1.00E+00	9.57E+00	1.96E+00	7.05E-01	1.10E+00	1.09E+00;
4.00	1262.00	1.05E+00	-2.36E-01	2.00E+00	1.00E+00	1.03E+01	1.96E+00	7.01E-01	1.10E+00	1.09E+00;
4.50	1262.00	1.07E+00	-2.37E-01	2.00E+00	1.00E+00	1.12E+01	2.07E+00	7.37E-01	1.10E+00	1.09E+00;
5.00	1262.00	8.30E-01	4.11E-03	2.00E+00	1.00E+00	4.26E+00	1.62E+00	7.53E-01	1.10E+00	1.08E+00;
5.50	1262.00	8.48E-01	-2.00E-02	2.00E+00	1.00E+00	4.02E+00	1.83E+00	7.85E-01	1.08E+00	1.08E+00;
6.00	1262.00	8.53E-01	-1.28E-02	2.00E+00	1.00E+00	3.56E+00	1.94E+00	8.23E-01	1.09E+00	1.09E+00;
6.50	1262.00	9.03E-01	-4.83E-02	2.00E+00	1.00E+00	3.17E+00	2.25E+00	9.14E-01	1.08E+00	1.07E+00;
7.00	1262.00	9.52E-01	-8.82E-02	2.00E+00	1.00E+00	3.14E+00	2.59E+00	9.79E-01	1.07E+00	1.06E+00;
7.50	1262.00	1.16E+00	-2.85E-01	2.00E+00	1.00E+00	2.56E+00	2.52E+00	1.04E+00	1.08E+00	1.06E+00;
8.00	1262.00	1.18E+00	-3.07E-01	2.00E+00	1.00E+00	3.10E+00	3.15E+00	1.15E+00	1.07E+00	1.05E+00 ]