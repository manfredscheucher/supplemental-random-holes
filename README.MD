# Holes in Random Point Sets

This website provides programs to verify the computer-assisted results from "Tight bounds on the expected number of holes in random point sets".

## Short description of the programs polyhedron.sage and ball.sage

The programs "polyhedron.sage" and "ball.sage" can be used to obtain estimates on the constants $$c_{3,4}^K=\lim_{n \to \infty} n^{-3}EH_{3,4}^K(n)$$ when $K$ is a Platonic solid or a ball. The idea is to repeatedly sample 3 points $p_1,p_2,p_3$ uniformly at random from $K$ and to compute the ratio of the so-called "Fischer triangle" spanned by $p_1,p_2,p_3$ which is contained inside $K$. This ratio coincides with $c_{3,4}^K$; see Section 3 in [our previous paper](https://doi.org/10.1002/rsa.21037). For the ball, run

```sage ball.sage```

For the tetrahedron, cube, octahedron, dodecahedron, or icosahedron, run

```sage polyhedron.sage [t/c/o/d/i]```

## Short description of the program test_planar_holes.py

The program samples $t$ sets of $n$ random points from a triangle, a square, or a disk and computes the average number of $k$-gons. More specifically, the points are sampled uniformly among the $[-g,+g] \times [-g,+g]$ grid. If the shape is "triangle", $x$-coordinate should be less than the $y$-coordinate, and if the shape is "disk", the distance to the origin is at most $g$. Sampled points which do not fulfill these conditions are omited. To run the program use

```python test_planar_holes.py g k n t [ball/triangle/square]```

## Downloads

*   The source code of the program "polyhedron.sage" [[download](polyhedron.sage)]
*   The source code of the program "ball.sage" [[download](ball.sage)]
*   The source code of the program "test_planar_holes.py" [[download](test_planar_holes.py)]