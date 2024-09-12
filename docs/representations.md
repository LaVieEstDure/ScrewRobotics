---

marp: true

---


# Screw Motions
```
```

---
## Representation

In 3D
$$
\hat{\mathcal{V}}=
\begin{bmatrix}
\hat{\omega} \\ \hat{v}
\end{bmatrix}=

\begin{bmatrix}
\hat{\omega}_{x} & \hat{\omega}_{y} & \hat{\omega}_{z} & \hat{v}_{x} & \hat{v}_{y} & \hat{v}_{z}
\end{bmatrix}^{\top}\in\R^{6}
$$
Given $\hat{s}$, $q$ and $h$:
- For a revolute joint:
    - $\hat{\omega} = \hat{s},\quad\hat{v}=-[\hat{s}]_{\times}q$

- For a prismatic joint:
    - $\hat{\omega}=0,\quad\quad\hat{v}=\hat{s}$

- For a helical joint:
    - $\hat{\omega} = \hat{s},\quad\hat{v}=-\hat{s}\times q+h\hat{s}$

---

In 2D
$$
\hat{\mathcal{V}}=
\begin{bmatrix}
\hat{\omega} \\ \hat{v}
\end{bmatrix}=

\begin{bmatrix}
\hat{\omega} \\ \hat{v}_{x} \\ \hat{v}_{y}
\end{bmatrix}\in\R^{3}
$$

- For a revolute joint at $q$:
    - $\hat{\omega}\in\{-1, 1\},\quad\hat{v}=-\hat{\omega}[1]_{\times}q=\hat{\omega}\begin{bmatrix}0 & -1 \\ 1 & 0\end{bmatrix}\begin{bmatrix}q_{x} \\ q_{y}\end{bmatrix}=\hat{\omega}\begin{bmatrix}-q_{y} \\ q_{x}\end{bmatrix}$

- For a prismatic joint with offset $p$:
    - $\hat{\omega}=0,\quad\quad\hat{v}=\frac{p}{\Vert p\Vert}$


---

