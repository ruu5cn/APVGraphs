# APVGraphs
APV Graphs

The code APV5 has two important methods - Sigma, and main. Sigma calculates the cross section for a given $x, Q^2$, and $s$. Main calculates all the APV data points given the settings it has and dumps it in a text file.

Things to be aware of in the code
  -Asymmetry values can differ greatly by slight changes in $sin^2(\theta_w)$ (SIN2THETAW). 
  -To change the polarization, change the variable $pol$ that's listed at the beginning of the main menu to "0","L", "H", or "LH". 
  -Beam energies are defined when calculating $s$, pay attention to these.
  -The data file is created and read with this pattern
      -First line is the number of x-points for every $Q^2$ value
      -Every other line after follows the paper \{$Q^2$ value\},\{x value\},\{Asymmetry value\}
