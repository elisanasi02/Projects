Here is a brief description of the projects I worked on.
# Speech Emotion Recognition
This was a group project for ST456 - Deep Learning at LSE. I worked in collaboration with Michele Bergami, Simone Moawad and Lucy Malabar.

This projected explores the realm of speech emotion recognition (SER), investigating the influence of data augmentation, model architectures, 
and dataset characteristics on classification accuracy. Our study introduces CNN and parallel CNN-RNN models for SER, examining their performance 
across various datasets and augmentation techniques. Our results demonstrate that data augmentation techniques, particularly pitch manipulation, 
significantly improve classification accuracy, particularly within the CNN framework. Also, we used Grad-CAM to provide insights into the models’ 
decision-making processes, revealing patterns in audio features associated with different emotions. Finally, we investigated the impact of gender 
bias and emotion intensity on model performance, revealing that models are more accurate under conditions of stronger emotino intensity and when 
the speaker is female.

My technical contribution was mainly on the literary review and the general outline of the project. 
In addition, in collaboration of my groupmates, I researched the appropriate datasets, I worked on the interpretation of the Grad-CAM results, on the 
implementation of the models and I curated the final report.

# Statistics Practitioners' Challenge
This was a group project for the LSE Statistics Practitioners' Challenge in collaboration of an important insurance company. I worked in collaboration with Michele Bergami, Simone Moawad, Barath Raaj and Anushka Agrawal. We analysed data about an insurance company, the main issue was related to imbalanced data. We explored the python library resreg, containing the possibility to implement several resampling techniques, from random undersampling to SMOTER and WERCS. At the same time we employed different predictive models: Poisson Regression, Random Forests, XGBoost and REBAGG.

# Telecom Churn + Coordinate Descent
This was a group project for ST443 - Machine Learning and Data Mining at LSE. I worked in collaboration with Michele Bergami, Simone Moawad, Omar Almutoteh and Giorgi Kvinikadze.

We assessed the performance of kNN, LDA, QDA, and Logistic Regression on a dataset from an Iranian telecommunication company, predicting customer churn. Additionally, we implemented and tested coordinate-descent algorithms for LASSO and Elastic Net across various simulated scenarios.

My technical contributions on this project are mainly on the first part of the project. More specifically, I worked on the implementation of the proposed classification techniques and I curated the final report.

# Attraction-Repulsion clustering
The present work aims to analyze a new fair clustering technique called Attraction Repulsion clustering.
In the first chapter, the concepts, tools, and methodologies typical of cluster analysis are outlined. Cluster analysis is applied in numerous fields and generally aims to identify homogeneous groups within a reference population.
In the second chapter, the focus is on the scenario where one or more of the considered variables are defined as sensitive; that is, variables for which homogeneity is not a desirable aspect. Examples include gender, ethnicity, or nationality. In this context, the Attraction Repulsion clustering technique, inspired by the principle of attraction and repulsion (hence the name) in electromagnetism, seeks to promote diversity, translated in terms of demographic parity, with respect to the sensitive variable. The chapter focuses on defining new aspects and tools, particularly the proposal of new dissimilarities, which characterize this technique.
After outlining this theoretical framework, chapter three addresses a case study: the technique is applied to a dataset containing information related to schools in the state of Massachusetts. The goal is to determine groups of schools that are geographically compact while also promoting diversity in terms of the ethnic composition of the students. The chapter concludes with reflections on the results obtained and a comparison with the application of a traditional cluster analysis technique.
