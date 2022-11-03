# Joint Sparse PCA

### Paper Presentation and Implementation

***Original Paper: *** https://dl.acm.org/doi/10.1016/j.patcog.2016.08.025 

JSPCA is designed by relaxing the orthogonal constraint of transformation matrix $Q$, introducing another transformation matrix $P$. It is imposing joint $ℓ_{2,1}$-norms on both loss term and regularization term. The proposed method has more freedom to jointly select the useful features for a low-dimensional representation and is robust to outliers. 

This work includes:

- An expositionary report of the original paper and the presentation of the talk given in Chennai Mathematical Insitute for the course Linear Algebra and Its Applications. 
- An implementation of the algorithm described in the paper on Breast Cancer data and comaparing it with older PCA variants.