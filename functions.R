#################################################################
## Functions
#################################################################
# density estimate of Student t distribution
tdistro <- function(X, Mu, sigma, a, b){
  -1* sum(( ((2*a+1.)/2.) * log(b + ((X - Mu)^2.)/(2.*sigma)) + .5 * log(sigma) ))
}   
#################################################################
# Function for linear interpolation in C++ to improve performance especially for large datasets
# Function assumes x is sorted in ascending order and there are no missing values in x
#################################################################
cppFunction('
NumericVector approxima(NumericVector x, NumericVector y, NumericVector xout) {
  int n = x.size();
  int nout = xout.size();
  NumericVector yout(nout, NA_REAL);
  
  int j = 0; // Index for the current interval
  
  for (int i = 0; i < nout; ++i) {
    // Find the interval containing xout[i]
    while (j < n - 1 && xout[i] > x[j + 1]) {
      ++j;
    }
    
    if (j < n - 1) {
      // Perform linear interpolation
      double slope = (y[j + 1] - y[j]) / (x[j + 1] - x[j]);
      double intercept = y[j] - slope * x[j];
      yout[i] = slope * xout[i] + intercept;
    } else if (xout[i] == x[j]) {
      yout[i] = y[j]; // Exact match found
    }
  }
  
  return yout;
}

NumericVector interpolate(NumericVector x, NumericVector y, NumericVector xout) {
  // Remove missing values
  NumericVector x_cleaned = x[!is_na(x) & !is_na(y)];
  NumericVector y_cleaned = y[!is_na(x) & !is_na(y)];
  
  // Sort x_cleaned and rearrange y_cleaned accordingly
  IntegerVector order = match(x_cleaned, clone(x).sort());
  NumericVector sorted_y = y_cleaned[order - 1];
  
  // Perform interpolation
  return approxima(x_cleaned, sorted_y, xout);
}'
)
#################################################################

#################################################################
# Function for linear extrapolation in C++
#################################################################
cppFunction('
NumericVector approxima_extrap(NumericVector x, NumericVector y, NumericVector xout) {
  int n = x.size();
  int nout = xout.size();
  NumericVector yout(nout, NA_REAL);
  
  int j = 0; // Index for the current interval
  
  for (int i = 0; i < nout; ++i) {
    // Find the interval containing xout[i]
    while (j < n - 1 && xout[i] > x[j + 1]) {
      ++j;
    }
    
    if (j < n - 1 && xout[i] >= x[0] && xout[i] <= x[n - 1]) {
      // Perform linear interpolation
      double slope = (y[j + 1] - y[j]) / (x[j + 1] - x[j]);
      double intercept = y[j] - slope * x[j];
      yout[i] = slope * xout[i] + intercept;
    } else if (xout[i] < x[0]) {
      // Perform linear extrapolation for values below the range of x
      double slope = (y[1] - y[0]) / (x[1] - x[0]);
      yout[i] = slope * (xout[i] - x[0]) + y[0];
    } else if (xout[i] > x[n - 1]) {
      // Perform linear extrapolation for values above the range of x
      double slope = (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]);
      yout[i] = slope * (xout[i] - x[n - 1]) + y[n - 1];
    } else if (xout[i] == x[j]) {
      yout[i] = y[j]; // Exact match found
    }
  }
  
  return yout;
}

NumericVector interpolate(NumericVector x, NumericVector y, NumericVector xout) {
  // Remove missing values
  NumericVector x_cleaned = x[!is_na(x) & !is_na(y)];
  NumericVector y_cleaned = y[!is_na(x) & !is_na(y)];
  
  // Sort x_cleaned and rearrange y_cleaned accordingly
  IntegerVector order = match(x_cleaned, clone(x).sort());
  NumericVector sorted_y = y_cleaned[order - 1];
  
  // Perform interpolation
  return approxima_extrap(x_cleaned, sorted_y, xout);
}'
)
#################################################################

