# Cross hierarchy propagation

## Paper 
- Gu, Y., Sainburg, L. E., Kuang, S., Han, F., Williams, J. W., Liu, Y., Zhang, N., Zhang, X., Leopold, D. A., & Liu, X. (2021). Brain Activity Fluctuations Propagate as Waves Traversing the Cortical Hierarchy. Cerebral Cortex , 31(9), 3986–4005. [[Paper](https://academic.oup.com/cercor/article-abstract/31/9/3986/6210040) , [github](https://github.com/YamengGu/the-cross-hierarchy-propagation)]

 
## Dataset
- Download
    - [Input data](https://drive.google.com/drive/folders/1KP_41R_qzuClfyd7r5CBQ6ULgfl-vfXd?usp=sharing)
    - [Simulation data](https://github.com/YamengGu/the-cross-hierarchy-propagation/tree/master/Data)


## Code 
- [projecting_onto_the_PG.py](https://github.com/wnsgus-SKKU/fMRI_project_hierarchy-propagation/blob/master/src/projecting_onto_the_PG.py) / [projecting_onto_the_PG.ipynb](https://github.com/wnsgus-SKKU/fMRI_project_hierarchy-propagation/blob/master/src/projecting_onto_the_PG.ipynb) : the code to project human input data onto a direction and calculate time-position correlation across time segments.

- [principal_delay.py](https://github.com/wnsgus-SKKU/fMRI_project_hierarchy-propagation/blob/master/src/principal_delay.py) / [principal_delay.ipynb](https://github.com/wnsgus-SKKU/fMRI_project_hierarchy-propagation/blob/master/src/principal_delay.ipynb): the code to calculate the principal delay profile.

- [finction.py](https://github.com/wnsgus-SKKU/fMRI_project_hierarchy-propagation/blob/master/src/function.py) : prerequisites

- Using as a [Colaboratory](https://colab.research.google.com/notebooks/welcome.ipynb?hl=ko-kr) / [how to?](https://ndb796.tistory.com/562)
## Simulation
- The [Connectome Workbench](http://www.humanconnectome.org/software/connectome-workbench) software is used to display the human brain surface.
- [Download data](https://github.com/YamengGu/the-cross-hierarchy-propagation/tree/master/Data)

## Presentation
### Notion
- [view](https://rectangular-flag-d72.notion.site/Cross-Hierarchy-Propagation-cb868980d647468b82195f2b58527103) 

### Contents

1. Introduction
    1. functional connectivity
    2. rsfMRI
    3. connectivity in rsfMRI
    4. about PG (background)  
      
2. Methods and results
    1. projecting the rsfMRI signals onto the PG direction [projecting_onto_the_PG.py](https://github.com/wnsgus-SKKU/fMRI_project_hierarchy-propagation/blob/master/src/projecting_onto_the_PG.py) / [projecting_onto_the_PG.ipynb](https://github.com/wnsgus-SKKU/fMRI_project_hierarchy-propagation/blob/master/src/projecting_onto_the_PG.ipynb)
    2. principal propagating direction in the human rsfMRI signals [principal_delay.py](https://github.com/wnsgus-SKKU/fMRI_project_hierarchy-propagation/blob/master/src/principal_delay.py) / [principal_delay.ipynb](https://github.com/wnsgus-SKKU/fMRI_project_hierarchy-propagation/blob/master/src/principal_delay.ipynb)
    3. simulation of rsfMRI signals with artificial propagations
    4. principal propagating direction in the monkey ECoG data
    5. fine-scale propagations within sensory modalities in the rsfMRI signals
    6. subcortical coactivations/deactivations associated with rsfMRI propagations
    7. modulation of cross-hierarchy propagations across different brain states
