# uSRD

This repository contains the code the for the "Utility function-based generalization of sum of ranking differences--country-wise analysis of greenhouse gas emissions" scientific article published in Ecological Indicators.


## Abstract:
The utility function-based sum of ranking differences (uSRD) method is proposed as a utility function-based multi-criteria decision analysis tool. Our idea is that the transformation functions can be represented by a utility function that can be aggregated with multi-attribute utility functions. We present a framework incorporating utility values as the basis for three different but interconnected analyses. The exemplary application focuses on greenhouse gas emissions and economic indicators of 147 countries. First, the uSRD is applied to the utility values to uncover the hidden relationships of the 40 indicators. A ranking of countries is established to see which sample performs the best and the worst in both emissions and economy. Lastly, mitigation actions are delegated to countries through a three-stage assignment that connects emissions to utilities, sectors, and mitigation actions. The results show that the uSRD excels as a support tool for decision-making.

## What it does:
 - Preprocesses CAIT dataset (results in [climate_data.xlsx](./climate_data.xlsx)), originating from [here](https://www.climatewatchdata.org/data-explorer/historical-emissions?historical-emissions-data-sources=All%20Selected&historical-emissions-gases=&historical-emissions-regions=&historical-emissions-sectors=&page=1) 
 - Performs SRD on the dataset for exploratory data analysis
 - Uses [Derringer transformation](./derringer.m) and ABC analysis to convert GHG emissions and economical data to perference values (in 0-1 range).
 - Performs SRD on the transformed data
 - Ranks country performance in GHG and ecological performance (better economical indicators with less GHG emissions), referred to as Clean development utility score in the publication. Generated data ([mapData.xlsx](./mapData.xlsx)) for maps with the scores ()
 - Links utility values with mitigation actions to form a recommender system to country-specific emissions portfolio.




## Contact
If you have any questions or proposals concerning this work, reach out to one of the following contracts:
 - Ádám Ipkovich - ipkovichadam@gmail.com
 - János Abonyi - janos@abonyilab.com
 - Viktor Sebestyén - sebestyenv@almos.uni-pannon.hu 
 - Károly Héberger - heberger.karoly@ttk.hu

## Cite as
If the work has been used in an application, we kindly ask you to cite as follows (note, this was uploaded during proof reading where whe have not yet been given doi.):
### General
Á. Ipkovich, K. Héberger, V. Sebestyén and J. Abonyi, Utility function-based generalization of sum of ranking differences--country-wise analysis of greenhouse gas emissions. Ecological Indicators. 2024.

[Citation]

### Bibtext
 @article{uSRD,
 author = { Ádám Ipkovich and Károly Héberger and Viktor Sebestyén and János Abonyi},
 title = {Utility function-based generalization of sum of ranking differences--country-wise analysis of greenhouse gas emissions},
 journal = {Ecological Indicators},
 volume = {},
 number = {},
 pages = {},
 year  = {2024},
 doi = {}
 }
 



## Previous version
Previous version can be found here, which was inaccessable during the proofing procedure. We recommend to check out other works of HUN-REN-PE Complex Systems Monitoring Research Group, here as well:
https://github.com/abonyilab/Preproc_SRD.git
If you used the previous code, please cite it as shown above.

