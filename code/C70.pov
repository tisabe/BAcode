#include "colors.inc"
#include "finish.inc"

global_settings {assumed_gamma 1 max_trace_level 6}
background {color White}
camera {orthographic
  right -18.51*x up 19.38*y
  direction 1.00*z
  location <0,0,50.00> look_at <0,0,0>}
light_source {<  2.00,   3.00,  40.00> color White
  area_light <0.70, 0, 0>, <0, 0.70, 0>, 3, 3
  adaptive 1 jitter}

#declare simple = finish {phong 0.7}
#declare pale = finish {ambient .5 diffuse .85 roughness .001 specular 0.200 }
#declare intermediate = finish {ambient 0.3 diffuse 0.6 specular 0.10 roughness 0.04 }
#declare vmd = finish {ambient .0 diffuse .65 phong 0.1 phong_size 40. specular 0.500 }
#declare jmol = finish {ambient .2 diffuse .6 specular 1 roughness .001 metallic}
#declare ase2 = finish {ambient 0.05 brilliance 3 diffuse 0.6 metallic specular 0.70 roughness 0.04 reflection 0.15}
#declare ase3 = finish {ambient .15 brilliance 2 diffuse .6 metallic specular 1. roughness .001 reflection .0}
#declare glass = finish {ambient .05 diffuse .3 specular 1. roughness .001}
#declare glass2 = finish {ambient .0 diffuse .3 specular 1. reflection .25 roughness .001}
#declare Rcell = 0.070;
#declare Rbond = 0.100;

#macro atom(LOC, R, COL, TRANS, FIN)
  sphere{LOC, R texture{pigment{color COL transmit TRANS} finish{FIN}}}
#end
#macro constrain(LOC, R, COL, TRANS FIN)
union{torus{R, Rcell rotate 45*z texture{pigment{color COL transmit TRANS} finish{FIN}}}
      torus{R, Rcell rotate -45*z texture{pigment{color COL transmit TRANS} finish{FIN}}}
      translate LOC}
#end

cylinder {< -6.48,  -6.66, -13.26>, <  6.79,  -6.66, -13.26>, Rcell pigment {Black}}
cylinder {< -6.48,   6.60, -13.26>, <  6.79,   6.60, -13.26>, Rcell pigment {Black}}
cylinder {< -6.48,   6.60,   0.00>, <  6.79,   6.60,   0.00>, Rcell pigment {Black}}
cylinder {< -6.48,  -6.66,   0.00>, <  6.79,  -6.66,   0.00>, Rcell pigment {Black}}
cylinder {< -6.48,  -6.66, -13.26>, < -6.48,   6.60, -13.26>, Rcell pigment {Black}}
cylinder {<  6.79,  -6.66, -13.26>, <  6.79,   6.60, -13.26>, Rcell pigment {Black}}
cylinder {<  6.79,  -6.66,   0.00>, <  6.79,   6.60,   0.00>, Rcell pigment {Black}}
cylinder {< -6.48,  -6.66,   0.00>, < -6.48,   6.60,   0.00>, Rcell pigment {Black}}
cylinder {< -6.48,  -6.66, -13.26>, < -6.48,  -6.66,   0.00>, Rcell pigment {Black}}
cylinder {<  6.79,  -6.66, -13.26>, <  6.79,  -6.66,   0.00>, Rcell pigment {Black}}
cylinder {<  6.79,   6.60, -13.26>, <  6.79,   6.60,   0.00>, Rcell pigment {Black}}
cylinder {< -6.48,   6.60, -13.26>, < -6.48,   6.60,   0.00>, Rcell pigment {Black}}
atom(<  0.08,  -3.64,  -9.88>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #0 
atom(< -1.22,  -1.66,  -0.05>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #1 
atom(< -1.06,   3.57,  -3.18>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #2 
atom(< -2.37,  -2.07,  -8.60>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #3 
atom(< -4.55,   6.24,  -1.20>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #4 
atom(<  0.95,  -2.50,  -1.48>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #5 
atom(<  2.47,  -1.48,  -6.54>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #6 
atom(<  0.47,   5.64,  -5.69>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #7 
atom(<  2.39,  -5.97,  -8.93>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #8 
atom(< -5.73,  -4.28,  -0.98>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #9 
atom(<  5.96,   2.81,  -3.55>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #10 
atom(< -0.35,   5.69,  -7.87>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #11 
atom(<  2.58,   1.96,  -5.32>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #12 
atom(< -3.53,   2.39,  -1.62>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #13 
atom(<  4.09,  -0.93,  -1.09>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #14 
atom(<  3.89,   2.96, -11.61>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #15 
atom(<  5.67,  -1.52,  -9.35>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #16 
atom(<  1.70,  -0.46,  -9.88>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #17 
atom(<  4.65,  -4.29, -10.25>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #18 
atom(<  2.27,  -0.31, -12.29>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #19 
atom(< -4.64,  -5.17,  -6.93>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #20 
atom(<  0.69,   6.04,  -5.52>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #21 
atom(<  0.51,   5.57,  -1.07>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #22 
atom(< -3.14,   2.39,  -0.47>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #23 
atom(< -3.47,   2.83,  -2.68>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #24 
atom(<  5.95,   4.66,  -7.67>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #25 
atom(< -6.41,  -6.20, -13.08>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #26 
atom(<  1.33,  -1.58, -11.74>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #27 
atom(<  4.99,   4.23,  -3.48>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #28 
atom(<  4.67,  -1.34,  -3.33>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #29 
atom(<  0.94,   1.47,  -9.31>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #30 
atom(<  1.51,  -1.02,  -6.12>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #31 
atom(< -0.02,  -0.03,  -1.09>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #32 
atom(< -2.99,  -6.00,  -7.13>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #33 
atom(< -1.06,  -0.20, -12.33>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #34 
atom(<  1.24,   6.43, -10.49>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #35 
atom(< -3.66,   6.08,  -3.89>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #36 
atom(<  1.39,  -1.61,  -4.79>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #37 
atom(<  0.19,  -5.47,  -9.32>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #38 
atom(<  4.51,  -4.34,  -1.27>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #39 
atom(<  1.32,  -6.44,  -8.36>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #40 
atom(< -5.92,  -4.98, -10.89>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #41 
atom(< -4.27,   5.06,  -1.25>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #42 
atom(< -3.36,  -1.26,  -7.01>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #43 
atom(<  4.53,  -2.79, -10.96>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #44 
atom(< -0.40,   5.31,  -5.37>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #45 
atom(< -1.01,  -1.28,  -9.17>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #46 
atom(< -2.26,  -1.16,  -4.80>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #47 
atom(<  4.07,  -6.17,  -3.52>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #48 
atom(<  3.77,  -3.42,  -0.20>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #49 
atom(<  2.47,   5.18, -10.04>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #50 
atom(< -5.79,   2.75,  -0.57>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #51 
atom(< -6.18,   5.45,  -3.60>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #52 
atom(<  6.19,  -3.04, -11.00>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #53 
atom(< -2.20,  -2.99,  -6.76>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #54 
atom(< -4.74,   3.13, -11.40>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #55 
atom(<  1.09,   0.41, -10.17>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #56 
atom(< -0.15,  -0.50,  -2.57>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #57 
atom(< -6.23,   2.79,  -9.93>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #58 
atom(<  2.28,   1.78, -12.51>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #59 
atom(<  3.06,  -2.16,  -5.97>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #60 
atom(<  3.39,  -4.30,  -9.94>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #61 
atom(< -4.70,  -4.66,  -0.81>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #62 
atom(<  6.41,  -4.77, -10.67>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #63 
atom(< -5.15,  -5.58,  -3.07>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #64 
atom(< -3.69,  -0.53, -12.15>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #65 
atom(< -2.08,   6.44,  -1.32>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #66 
atom(<  5.00,  -1.57,  -3.97>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #67 
atom(< -5.99,  -0.01, -11.75>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #68 
atom(<  3.25,  -0.87,  -8.93>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #69 
