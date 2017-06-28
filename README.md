# R packages 'nonlinearICP' and 'CondIndTests'
R Code for 'nonlinearICP' and 'CondIndTests'. 

## CondIndTests

Code for a variety of nonlinear conditional independence tests: 

  * Kernel conditional independence test (Zhang et al., UAI 2011),
  * Residual Prediction test (based on Shah and Buehlmann, arXiv 2015),
  * Invariant environment prediction,
  * Invariant target prediction,
  * Invariant residual distribution test,
  * Invariant conditional quantile prediction (all from Heinze-Deml et al., <arXiv:1706.08576>).
 
### Installation

#### From Github with `devtools`
```r
devtools::install_github("christinaheinze/nonlinearICP-and-CondIndTests/CondIndTests")
```

## nonlinearICP

Code for 'nonlinear Invariant Causal Prediction' to estimate the 
    causal parents of a given target variable from data collected in
    different experimental or environmental conditions, extending
    'Invariant Causal Prediction' from Peters, Buehlmann and Meinshausen (2016)
    to nonlinear settings. For more details, see C. Heinze-Deml, J. Peters and 
    N. Meinshausen: 'Invariant Causal Prediction for Nonlinear Models', 
    <arXiv:1706.08576>.
 
### Installation

#### From Github with `devtools`
```r
devtools::install_github("christinaheinze/nonlinearICP-and-CondIndTests/nonlinearICP")
```

## References
 If you are using these packages, please cite C. Heinze-Deml, J. Peters and 
    N. Meinshausen: 'Invariant Causal Prediction for Nonlinear Models', [arXiv:1706.08576](http://arxiv.org/abs/1706.08576).
