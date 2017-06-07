# ScagnosticsSharp

This project is C# migrated version of original scagnostics R package :

https://cran.r-project.org/web/packages/scagnostics/index.html  
https://github.com/cran/scagnostics

## Theoretical Background (Paper)

Wilkinson L., Anand, A., and Grossman, R. (2006). High-Dimensional visual analytics:  
Interactive exploration guided by pairwise views of point distributions.  
IEEE Transactions on Visualization and Computer Graphics, November/December 2006 (Vol. 12, No. 6) pp. 1363-1372.

## Requirements

#### For using pre-compiled binary :

* ≥ .Net Framework 4.0

#### For building project :

* ≥ .Net Framework 4.0

* ≥ Visual Studio 2013

## Usage

1. Create "Scagnostics" class object.
2. Run "Compute()" method of Scagnostcs object.
3. Get a double type array as a result.

* To get exactly the same results as the R package, you have to run "LoadJavaRandomNumber()" method before running Compute(). However, theoretically, there is no problem without doing this.

* You can see more detailed usage by referring to "ScagnosticsSharp.Test" project.
