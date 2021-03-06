# genosolver

gensolver is the official solver of the GENO (GENeric Optimization) framework. Using GENO requires a working installation of genosolver. Solvers for your optimization problem can be obtained from the official geno4neurips website:

  http://geno4neurips.pythonanywhere.com/

### Requirements
- python3
- pip
- c++ compiler that supports the  c++11 standard.

Tested:
- **Linux**: Ubuntu 14/16/18 with anaconda python3.6/3.7 and g++7.4.0, clang++ 3.8.0
- **Windows**: Windows 10 with with anaconda python3.6/3.7 and MSVS 17/18/19
- **Mac**: macOS Mojave with anaconda python3.6/3.7 and clang++ 3.8


### Installation via pip
Clone and move to the genosolver repository:

```sh
pip install -r requirements.txt && python setup.py install
```

You can test your installation as follows:
```sh
python scripts/sample-02e.py
```

### Usage
**You can download a solver from the official geno4neurips website:**
  http://geno4neurips.pythonanywhere.com/

#### Example - l2-regularized Logistic Regression

Modeling a l2-regularized Logistic Regression can be achieved in the following way:

```
parameters
   Matrix X
   Scalar c
   Vector y
variables
   Vector w
min
   0.5 * w' * w 
   + c * sum(log(exp((-y) .* (X * w))
   + vector(1)))
```
Enter this into the input field on http://geno4neurips.pythonanywhere.com/ and click the "Download Solver" button. Once the download has completed, you can use the solver via 

```sh
python solver.py
````

You can feed your own data into the solver by replacing the line:

```python
X, y = generateRandomData()
```

#### Example - Support Vector Machine

```
parameters
   Matrix K symmetric
   Scalar c
   Vector y
variables
   Vector a
min
   0.5 * (a.*y)' * K * (a.*y) - sum(a)
st
   a >= 0
   y' * a == 0
```
#### Example - Non-negative Least Squares

```
parameters
   Matrix A
   Vector b
variables
   Vector x
min
   norm2(A*x - b).^2
st
   x > 0
```
