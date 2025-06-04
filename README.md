# Quick Loss Estimation Tool (QLET)

**QLET** is a MATLAB-based tool designed for rapid seismic risk assessment of Canadian buildings. It enables users to compute seismic loss exceedance curves and Value-at-Risk (VaR) metrics efficiently for both site-specific and regional scenarios. QLET integrates:

- National seismic hazard models from the Geological Survey of Canada (GSC)
- Building exposure and vulnerability models
- Global and regional Vs30 data
- Monte Carlo simulation-based loss estimation
- Tail approximation methods for hazard curve estimation

This tool was developed to facilitate regional-scale risk evaluations within practical timeframes, making it particularly applicable to emergency preparedness, urban planning, and insurance analysis.

## Reference

For detailed methodology and case studies using QLET, please refer to the following publication:

**Momeni, P.; Goda, K.; Sirous, N.; Molnar, S. (2025).**  
*Rapid Computation of Seismic Loss Curves for Canadian Buildings Using Tail Approximation Method*.  
**GeoHazards, 6(2), 26.**  
[https://doi.org/10.3390/geohazards6020026]

A copy of the published paper is also included in the root directory of this repository.

## Repository Structure

- `Matlab_codes/`: Core MATLAB scripts for hazard and risk computations
- `App/`: GUI-based implementation of QLET using MATLAB
- `memo_for_data.txt`: Data access note (please read before using)

## Important Notes

- Some datasets (e.g., NRCan exposure data and Nastev et al. 2016 Vs30 data) have been removed from the repository. If you are interested in accessing these datasets, please contact:
- Dr. Katsuichiro Goda (kgoda2@uwo.ca)
- Dr. Payam Momeni (pmomeni@uwo.ca)
