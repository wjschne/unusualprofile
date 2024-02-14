# unusualprofile 0.1.3 _2022-12-07_

* Added `zSEE` to output list from `cond_maha`. It contains the standardized standard error of the estimate.
* Fixed bug where `SEE` slot in output lists from `cond_maha` was incorrectly reporting standardized standard error of the estimate---what what is now `zSEE`. The `SEE` slot now reports the (unstandardized) standard error of the estimate.


# unusualprofile 0.1.2 _2022-12-07_

* Fixed plot.cond_maha bugs related to having scores with different means and standard deviations


# unusualprofile 0.1.1 _2022-12-01_

* Fixed a bug that prevented independent variables and independent composites to be used together as predictors
* removed reduction statistics in plot.cond_maha


# unusualprofile 0.1.0 _2021-05-12_

* First submission to CRAN
