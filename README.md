# NMFDAD
NMFDAD can be used for the decomposition of LC-DAD data, in oder to obtain the individual UV-vis spectrum for each substance in complex dissolved mixtures detected by UHPLC-DAD-UHRMS (Ultra High Performance Liquid Chromatography-Diode Array Detector-Ultra High Resolution Mass Spectrometry). Overlapping spectra are resolved by a cnmf (Constrained Non-negative Matrix Factorization) model. 

## Overview
LC-DAD presents absorption of the eluates in both wavelength and retention time domain. The acquired data are matrixes with n rows and m columns, where m refers to retention time points, n is wavelength points. The data matrix (X) can be decomposed into a spectrographic (S) matrix factor and an eluting peak (P) matrix factor:  
X(n×m) =S(n×r) × P(r×m),  
where r is corresponding to the number of absorbing sub-stances. Each column of S represents the absorption spectrum of a certain substances, while each row of P represents its eluting peak. Since all the elements in X, P, and S matrices are all nonnegative, the factorization is indeed a NMF problem.  

For more infomation please refer to our paper:
Xu, N.; Hu, M.; Li, X.; Song, K.; Qiu, Y.; Sun, H. X.; Wang, Y.; Zeng, L.; Li, M.; Wang, H.; Hu, S.; Zong, T.; Bai, Y.; Zhang, Z.; Li, S.; Shuai, S.; Chen, Y.; Guo, S., Resolving Ultraviolet–Visible Spectra for Complex Dissolved Mixtures of Multitudinous Organic Matters in Aerosols. Anal. Chem. 2024, 96, (5), 1834-1842.
DOI: 10.1021/acs.analchem.3c02700

## Operation
1.Prepare data as the main.m required.  
2.Add the function in cnmf.m to MATLAB.  
3.Run the script in main.m.  
