# Nozzle_MOC
Design Algorithm of a supersonic nozzle, using the Method of characterictics

Nozzle_MOC V1.0 29-03-2022
===================================

This algorithm can be used to design a supersonic bell nozzle, and generate the Mach number contour plot, using the method of characteristics.

Follow these steps in order to use the tool,
1. Select Supersonic Mach Number (eg. Me = 2). As the Method of characteristics is applied to hyberbolic ODEs, the flow must be always supersonic (Me > 1). Very high Exit Mach number values may result in errors.
2. Select the value of γ (Ratio of Specific Heats) (eg. Gamma = 1.4). Gamma can take values < 1 however very low values may resulted in errors.
3. Select the number of characteristics (eg. n = 20) . The highest is the number of characteristics the more accurate are the nozzle results, which comes with the cost of higher computational time. Number of characteristics between 20 - 50 is fine. CHECK the the number of characteristics is always an integral number and it is not negative or zero.
4. Select the coordinates (x0,y0) at the upper point at the throat (eg. x0 = 0, y0 = 1)
5. Select the reference coordinates in order to generate dimensionless coordinates (eg. y_ref = 1).
6. Select any other plot properties' parameters as they are described in the functions of Nozzle_MOC.py file.
7. Check that Nozzle_MOC.py and PrandltMeyer.py are at the same path.
8. Run and generate the resulted plots.
