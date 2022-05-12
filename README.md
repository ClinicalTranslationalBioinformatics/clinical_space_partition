# CSP
CSP (Clinical Space Partition) is a method for comparing simultaneously an arbitrary number of predictors across the clinical space.
The program provides the optimal division of the clinical space and for each predictor, gives the clinical space it occupies.

There are two versions of the program, one that takes into account the coverage of the predictor (CSP-rej) and another that does not consider it (CSP-norej).

In the case of CSP-norej, each predictor is defined by 2 parameters: sensitivity and specificity;
whereas in the case of CSP-rej, each predictor is defined by 3 parameters: sensitivity, specificity and coverage.

Additionally, in both programs, there is another parameter called *rho* (see accompanying paper), 
related to the frequency of the minority class (which corresponds to the pathogenic variants in the article). 
We use 0.5 as a default value, but the user can modify it between 1·10<sup>-5</sup> and 1.

# Requirements
We recommend using Python 3.9 from http://www.python.org. CSP is supported and tested on it.

CSP requires [Numpy](https://numpy.org/),  [Shapely](https://pypi.org/project/Shapely/) and [SymPy](https://www.sympy.org) packages which can be installed with pip:

```
pip3 install numpy shapely sympy
```

# Usage

## CSP-rej

To get the fraction of the clinical space corresponding to a set of predictors, run the following command:

```
python3 csp_rej.py file.config
```

For example, in the case of 3 predictors, the config file should include:

```
[rho]
rho=rho

[predictors]
predictor1=sensitivity1,specificity1,coverage1
predictor2=sensitivity2,specificity2,coverage2
predictor3=sensitivity3,specificity3,coverage3
```

There is an example of config file in the demo directory. 
To run it, type in your terminal:

```
python3 csp_rej.py ../demo/csp-rej.config
```
The results will be prompted in your terminal. You can find a copy of them in the file demo/csp-rej.txt

## CSP-norej

To get the fraction of the clinical space corresponding to a set of predictors, run the following command:

```
python3 csp_norej.py file.config
```

For example, in the case of 3 predictors, the config file should include:

```
[rho]
rho=rho

[predictors]
predictor1=sensitivity1,specificity1
predictor2=sensitivity2,specificity2
predictor3=sensitivity3,specificity3
```

There is an example of config file in the demo directory. 
To run it, type in your terminal:

```
python3 csp_norej.py ../demo/csp-norej.config
```
The results will be prompted in your terminal. You can find a copy of them in the file demo/csp-norej.txt