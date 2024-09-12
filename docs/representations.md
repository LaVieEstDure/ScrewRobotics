## Representations

### In 3D
```math
\hat{\mathcal{V}}=
\begin{bmatrix}
\hat{\omega} \\\ \hat{v}
\end{bmatrix}=
\begin{bmatrix}
\hat{\omega}_{x} \\\ \hat{\omega}_{y} \\\ \hat{\omega}_{z} \\\ \hat{v}_{x} \\\ \hat{v}_{y} \\\ \hat{v}_{z}
\end{bmatrix}\in\mathbb{R}^{6}
```

Given $\hat{s}$, $q$ and $h$:
- For a revolute joint:
```math
\hat{\omega} = \hat{s},\quad\hat{v}=-[\hat{s}]_{\times}q
```

- For a prismatic joint:
```math
\hat{\omega}=0,\quad\quad\hat{v}=\hat{s}
```

- For a helical joint:
```math
\hat{\omega} = \hat{s},\quad\hat{v}=-\hat{s}\times q+h\hat{s}
```

---

### In 2D
```math
\hat{\mathcal{V}}=
\begin{bmatrix}
\hat{\omega} \\ \hat{v}
\end{bmatrix}=
\begin{bmatrix}
\hat{\omega} \\ \hat{v}_{x} \\ \hat{v}_{y}
\end{bmatrix}\in\mathbb{R}^{3}
```

- For a revolute joint at $q$:
```math
\hat{\omega}\in\{-1, 1\},\quad\hat{v}=-\hat{\omega}[1]_{\times}q=\hat{\omega}\begin{bmatrix}0 & -1 \\ 1 & 0\end{bmatrix}\begin{bmatrix}q_{x} \\ q_{y}\end{bmatrix}=\hat{\omega}\begin{bmatrix}-q_{y} \\ q_{x}\end{bmatrix}
```

- For a prismatic joint with offset $p$:
```math
\hat{\omega}=0,\quad\quad\hat{v}=\frac{p}{\Vert p\Vert}
```


---

