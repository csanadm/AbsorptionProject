# AbsorptionProject
Calculating sound absorption in air, water and various other materials. It requires gcc (tested with 5.4), GNU Make (tested with 4.1), ROOT (tested with 6.15), the latter also including the Minuit2 component (for the fitting).

## Description
This project is based on sound absorption formulas. For the case of air, it is based on:
- ANSI Standard [S1-26:1995 (ISO 9613-1:1996)](https://puc.sd.gov/commission/dockets/electric/2019/el19-003/KMExhibit9.pdf)
- H.E. Bass et al [J.Acoust.Soc.Am., 97(1), 680 (1995)](https://calhoun.nps.edu/handle/10945/62134)
For water, it is based on:
- Francois & Garrison, [J.Acoust.Soc.Am., 72(3), 896](https://asa.scitation.org/doi/10.1121/1.388170) and [J.Acoust.Soc.Am., 72(6), 1879](https://asa.scitation.org/doi/10.1121/1.388673), both 1982.
- Fisher & Simmons, [J.Acoust.Soc.Am., 62, 558, 1977](https://asa.scitation.org/doi/10.1121/1.381574).
- Ainslie & McColm, [J.Acoust.Soc.Am., 103(3), 1671, 1998](https://asa.scitation.org/doi/10.1121/1.421258).
For other, construction related materials, it is based on [measurements](https://www.acoustic.ua/st/web_absorption_data_eng.pdf) of a Kiev-based company Acoustic Trafic.

It consists of codes of two kind:
- fitting the original formulas with a simplified version and determining the parameters
- using the simple formulas

## File content
- `README.md`: This README file
- `Makefile`: Using make `<basename>.exe`, it will create an executable from any `<basename>.cc`
- `base_formulas_original.h`: This contains all basic formulas from the original papers
- `air_absorption_4component_fitter.cc`: Fitting the original air absorption formulas with a simplified, piecewise formula, with 4 components
- `air_absorption_5component_fitter.cc`: Fitting the original air absorption formulas with a simplified, piecewise formula, with 5 components
- `water_absorption_general.cc`: Fitting the original water absorption formulas with a simplified, piecewise formula, with 3 components
- `water_air_simpifield_formula.h`: Calculation of the piecewise simple formula for water and air
- `other_materials_ai.h`: Calculation of the piecewise simple formula for other (construction-related) materials
- `plot_simple_example.cc`: Creating example plots using the simplified formulas
- `air_approximation_T20_H30_P1.png`: Example plot for air
- `water_approximation_T8_S35_D500_pH8.png`: Example plot for water
- `othermat_approximation_nmat5.png`: Example plot for a construction material

## Example plots as results
![air_approximation_T20_H30_P1](https://user-images.githubusercontent.com/38218165/176395940-dab9d8be-c048-4e18-9d11-c1f9f550754b.png|width=300)
![water_approximation_T8_S35_D500_pH8](https://user-images.githubusercontent.com/38218165/176395949-d7313fea-6fef-4019-bde3-64d96b9587ab.png|width=300)
![othermat_approximation_nmat5](https://user-images.githubusercontent.com/38218165/176506240-391cf60d-4407-462d-8d44-5e1257e6a594.png|width=300)
