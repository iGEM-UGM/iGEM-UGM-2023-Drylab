iGEM UGM 2023: ColDBlu Project
Dry Lab Section
==============================

Kinetic Modelling also can be accessed through this [Tableau Visualization](https://public.tableau.com/app/profile/kayla2917/vizzes)

Project Organization
------------

    ├── README.md          <- The description of this project.
    ├── Kinetic-Modelling
    │   ├── Kinetic Modeling Equation for LIRA.pdf  <- Equation of Kinetic Modelling and constanta
    │   ├── Kinetic_Modeling_for_LIRA.ipynb         <- Kinetic modelling for LIRA with visualization
    │   └── kinetic_modelling_for_tableau.ipynb     <- Kinetic modelling for combination setup
    │  
    └── LiRA-Optimization
        ├── autogluon-nupack.ipynb                              <- Processing using AutoGluon
        ├── final.csv                                           <- Gathering parameter result
        ├── NUPACK Design Parameter Settings for LIRA.pdf       <- Guide to setup Design Job
        ├── pipeline.ipynb                                      <- Gathering Sequence, Parameters, and Linear Regression Multivariate Regression
        └── table_data_a10996f7-3a7c-49eb-ba3c-bf682effc883.csv <- Raw data

--------
# Installations
Some of codes used are run on Google Colaboratory and Kaggle, so some of parts (might include pip install) have been included the installations process or using pre-installed package in cloud-based server. However, NUPACK can't run on cloud-based server such as Colab/Kaggle due to its requirement to be installed on top of Linux-based Environment. Follow these guide to install main library used in this repo:
* NUPACK: https://docs.nupack.org/start/
* AutoGluon: https://auto.gluon.ai/stable/install.html
