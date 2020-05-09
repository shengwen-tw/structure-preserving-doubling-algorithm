# Structure-preserving Doubling Algorithm (MATLAB)

A fast CARE (Continuous Algebraic Riccati Equation) solving algorithm implemented in MATLAB.

## Reference paper

* [A structure-preserving doubling algorithm for continuous-time algebraic Riccati equations](https://www.sciencedirect.com/science/article/pii/S0024379504004434)

## Demonstrations

### 1. sda.m

A simple implementation of Structure-preserving Doubling Algorithm (SDA) compare to MATLAB's CARE function.

### 2. pendulum_matlab_care.m

A pendulum LQR control simulation implemented with MATLAB's CARE function.

### 3. pendulum_sda.m

A demonstration of speed enhancement using SDA for pendulum LQR control problem.

In second demonstration, the simulation took 41.49 seconds on author's computer, and only took 6.72 seconds for SDA version. The result came out to be 6.1 times faster.

For detailed deriviation of pendulum LQR control, please refer to this [pdf](https://github.com/shengwen-tw/structure-preserving-doubling-algorithm/raw/master/lqr.pdf).

 <img src="https://github.com/shengwen-tw/structure-preserving-doubling-algorithm/blob/master/math.png?raw=true" width="90%" height="90%">


