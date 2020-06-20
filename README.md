# AKI Phenogroup External Validation Open Source Implementation
The external validation experiment of aki phenogroup research

## Prerequisite
### Data
We assume the you are already able to access the [MIMIC-III (V1.4)](https://mimic.physionet.org/about/mimic/), which is a publicily available electronic health record (EHR) databse which contains over 40 thousand patients who stayed in critical care units of the Beth Israel Deaconess Medical Center between 2001 and 2012.
If you are not able to access MIMIC-III, please complete a required training course and then submit an application in advance at [request link](https://mimic.physionet.org/gettingstarted/access/). Once the application is approved, you can download the MIMIC-III freely.
Please ensure your computer has enough storage space because MIMIC-III is a large dataset (45 GiB space is required in this study)

### Environment
Please install Python 3.7 environment, as well as lifelines (Version 0.24.8), matplotlib (3.2.1), numpy (1.18.4), pandas (1.0.4), scikit-learn (0.23.1), and scipy (1.4.1) packages in advance.
