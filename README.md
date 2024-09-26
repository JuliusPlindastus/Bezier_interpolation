# Bezier_interpolation
Matlab code to make a Bezier interpolation, calculate curvature and improve
A lot of inspiration from this website: https://pomax.github.io/bezierinfo

## BezierFit 
It interpolates 2D curves with a bunch of known points using Bezier functions.
- QXY is the position of the known points.
  Preferably in complex notation, there are still some bugs with the two dimension matrices.
  Also check if your points are in columns or rows. This is also something to be changed in the future.
- options is a structure to give more info to the fit. the most important input:
  - options.degree: int. degree of Bezier polynomial (by default: 3) It makes sense to try up to 5-6.
  - options.anchor_start and options.anchor_last: boolean. If you're sure the first and last point are at the right position, the BezierFit will impose to stay on these points. (Usually in Bezier fitting, the first and last points are quite important)
  - options.optim: false, true, 'tangent'. By default: false
    - False, it will make an equidistant distance between the control points. This is the one that seems the clearest but does not necessarily converge to a nice fitting.
    - True, optimisation added to try to minimize the error between the known points and the interpolation. It fits better, but the control points are a bit less relevant to analyze. One may obtain very sharp curvature change, not ideal.
    - 'tangent', optimisation a bit more specific to fit the tangent at the extremities. If we are sure of the first two and last two points, we would like the Bezier curve to have the same tangent (slope of the first and last 2 control points) than what we have between the known points. May provide dodgy results sometimes.
  - options.weight: array the same length than QXY. If some of your points have less uncertainties, or more important to fit than others, you can increase their weight. (by default all of the weights are at one).
  
## BezierCurvature
It calculates the second derivative to compute the curvature of the curve. The Bezier interpolation has the advantage to normally provide a smooth curvature along the line.
It can be useful if you want to extrapolate the curve using the curvature evolution.

## BezierConstruction
Based on the control points Pz, you can recompute the Bezier interpolation. It can help to reduce large data with control points only.
Now the degree is in the options parameter, but it should be obvious to the code (not implemented yet) that the degree should be the length of the Pz-1.

