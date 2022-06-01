# Data for replication

Here is the data needed for replicating the experiments. 


## Artificial example 1 ("Explicit example with two assets" in paper)

strikes = [100 110; 102 107]


prices = [12 3; 10 6]

weights = [1 / 2 1 / 2]


## Artificial example 2

strikes = [95 100 110 115;
    96 102 107 112]

prices = [15.5 11 3 1.5;
    14.5 9 6 4]

weights = [1 / 2 1 / 2]


## Artificial example ("Varying strikes" in paper)

strikes = [90 95 100 110 120;
    90 96 102 107 115]
    
prices = [20 15.5 12 5.5 1;
    20.5 15 10 6 0.75]
  
weights = [1 / 2 1 / 2]


## Artificial exampel ("Currency basket" in paper)

strikes = [135.5 138.5;
    116 119]
    
prices = [2.77 1.17;
    2.21 0.67]
    
weights = [2 / 3 1 / 3]


## Example Bertsimas Popescu (data from The Wall Street Journal, July 7, 1998, Microsoft Call options)

strikes = [95 100 110 115 120]

prices = [12.875 8.375 1.875 0.625 0.25]

weights = [1]

## Strike and prices of tech stock call options (data from https://www.nasdaq.com/market-activity)

#AAPL

strikes = [120 130 145 160 170]

prices = [45.2 35.7 21.75 9.1 3.35]

weights = [1]

#FB

strikes = [155 170 180 190 200 210]

prices = [52.7 38.5 29.85 22 14.75 9.15]

weights = [1]

#NVDA

strikes = [175 180 190 195 227.5]

prices = [57.9 53.2 43.85 39.35 10.75]

weights = [1]

#QCOM

strikes = [130 145 157.2 167.5 175]

prices = [35.35 20.5 8.8 2.32 0.47]

weights = [1]

## Data from Table 6 in paper to produce Table 7

#AAPL, FB, NVDA, QCOM

strikes = [120 130 145 160 170;
    155 170 180 190 200;
    175 180 190 195 227.5;
    130 145 157.5 167.5 175]
    
prices = [45.2 35.7 21.75 9.1 3.35;
    52.7 38.5 29.85 22 14.75;
    57.9 53.2 43.85 39.35 10.75;
    35.35 20.5 8.8 2.32 0.47]

weights = [1 / 4 1 / 4 1 / 4 1 / 4]
