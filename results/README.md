# Here we present the result obtain by running the code for different inputs. 

Note: In all examples for the outer range we set *B = 400* and *M = 200 000*. The results we give below are the solution to the first level of the proposed hierarchy. 


## Univariate example (Bertsimas Popescu)
First consider the data from Table 1 in the paper:

strikes = [95 100 110 115 120]

prices = [12.875 8.375 1.875 0.625 0.25]

weights = [1]

For *K = 105* we get an *upper bound* of 5.125 and a *lower bound* of 3.875.

## Explicit example with two assets

strikes = [100 110;
          102 107]

prices = [12 3;
          10 6]

weights = [1 / 2 1 / 2]

For *K = 105* we get an *upper bound* of 7.4 and a *lower bound* of 2.387.

## Varying strikes 

strikes = [90 95 100 110 120; 90 96 102 107 115]

prices = [20 15.5 12 5.5 1; 20.5 15 10 6 0.75]

weights = [1 / 2 1 / 2]

We present the results in the following form *(K, lower bound, upper bound)*

(90, 16.875, 20.25)

(95, 12.792, 15.7)

(100, 8.708, 11.55)

(105, 4.625, 8.016)

(110, 1.675, 4.75)

(115, 0.0, 2)

## Currency basket option

strikes = [135.5 138.5; 116 119]

prices = [2.77 1.17; 2.21 0.67]

weights = [2 / 3 1 / 3]


We present the results in the following form *(K, lower bound, upper bound)*

(100, 1.4933, 31.5834)

(105, 1.2599, 26.5833)

(110, 1.0266, 21.5833)

(115, 0.7933, 16.5833)

(120, 0.56, 11.5833)

## Basket option on tech stocks 

strikes = [120 130 145 160 170; 155 170 180 190 200; 175 180 190 195 227.5; 130 145 157.5 167.5 175]

prices = [45.2 35.7 21.75 9.1 3.35; 52.7 38.5 29.85 22 14.75; 57.9 53.2 43.85 39.35 10.75; 35.35 20.5 8.8 2.32 0.47]

weights = [1 / 4 1 / 4 1 / 4 1 / 4]

We present the results in the following form *(K, upper bound, lower bound)*

(140, 52.79, 46.26)

(150, 42.89, 36.26)

(160, 33.48, 26.27)

(170, 24.53, 16.28)

(180, 15.68, 6.28)

(190, 8.51, 0.0)

(200, 6.99, 0.0)

## Boyle and Lin example

For the means given by [44.21, 44.21, 44.21] and covariance matrix given by [184.04 164.88 164.88; 164.88 184.04 164.88; 164.88 164.88 184.04]
we find for the smallest upper bound on the price of the call on max option for different values of *K* the following values:

We present the results in the following form *(K, upper bound, lower bound)*

(30, 21.51, 14.21)

(35, 17.17, 9.21)

(40, 13.2, 4.21)

(45, 9.84, 0.0)

(50, 7.3, 0.0)

## Example inner range

Here we give the result for the example using the inner range. We give results for the levels *r = 2,...,7* of the hierarchy.

strikes = [100, 110]

prices = [8.375, 1.875]

k = 105

We present the results in the following form *(level, epsilon, upper bound, lower bound)*

(2, 0.0273, 5.1279, 5.122)

(3, 0.02525, 5.1366, 5.1136)

(4, 0.022125, 5.1288, 5.1221)

(5, 0.01755, 5.1264, 5.1251)

(6, 0.0161, - ,  4.224)

(7, 0.0161, - , 3.3522)

