
# This file was automatically created by The UFO_usermod        

from object_library import all_lorentz, Lorentz
UUS1 = Lorentz(name = 'UUS1',
               spins = [-1, -1, 1],
               structure = '1')


UUV1 = Lorentz(name = 'UUV1',
               spins = [-1, -1, 3],
               structure = 'P(3,2) + P(3,3)')


SSS1 = Lorentz(name = 'SSS1',
               spins = [1, 1, 1],
               structure = '1')


FFS1 = Lorentz(name = 'FFS1',
               spins = [2, 2, 1],
               structure = 'ProjM(2,1)')


FFS2 = Lorentz(name = 'FFS2',
               spins = [2, 2, 1],
               structure = 'ProjM(2,1) - ProjP(2,1)')


FFS3 = Lorentz(name = 'FFS3',
               spins = [2, 2, 1],
               structure = 'ProjP(2,1)')


FFS4 = Lorentz(name = 'FFS4',
               spins = [2, 2, 1],
               structure = 'ProjM(2,1) + ProjP(2,1)')


FFV1 = Lorentz(name = 'FFV1',
               spins = [2, 2, 3],
               structure = 'Gamma(3,2,1)')


FFV2 = Lorentz(name = 'FFV2',
               spins = [2, 2, 3],
               structure = 'Gamma(3,2,-1)*ProjM(-1,1)')


FFV3 = Lorentz(name = 'FFV3',
               spins = [2, 2, 3],
               structure = 'Gamma(3,2,-1)*ProjP(-1,1)')


FFV4 = Lorentz(name = 'FFV4',
               spins = [2, 2, 3],
               structure = 'Gamma(3,2,-1)*ProjM(-1,1) + 4*Gamma(3,2,-1)*ProjP(-1,1)')


VSS1 = Lorentz(name = 'VSS1',
               spins = [3, 1, 1],
               structure = 'P(1,2) - P(1,3)')


VVS1 = Lorentz(name = 'VVS1',
               spins = [3, 3, 1],
               structure = 'Metric(1,2)')


VVV1 = Lorentz(name = 'VVV1',
               spins = [3, 3, 3],
               structure = 'P(3,1)*Metric(1,2) - P(3,2)*Metric(1,2) - P(2,1)*Metric(1,3) + P(2,3)*Metric(1,3) + P(1,2)*Metric(2,3) - P(1,3)*Metric(2,3)')


SSSS1 = Lorentz(name = 'SSSS1',
                spins = [1, 1, 1, 1],
                structure = '1')


VVSS1 = Lorentz(name = 'VVSS1',
                spins = [3, 3, 1, 1],
                structure = 'Metric(1,2)')


VVVV1 = Lorentz(name = 'VVVV1',
                spins = [3, 3, 3, 3],
                structure = 'Metric(1,4)*Metric(2,3) - Metric(1,3)*Metric(2,4)')


VVVV2 = Lorentz(name = 'VVVV2',
                spins = [3, 3, 3, 3],
                structure = 'Metric(1,4)*Metric(2,3) + Metric(1,3)*Metric(2,4) - 2*Metric(1,2)*Metric(3,4)')


VVVV3 = Lorentz(name = 'VVVV3',
                spins = [3, 3, 3, 3],
                structure = 'Metric(1,4)*Metric(2,3) - Metric(1,2)*Metric(3,4)')


VVVV4 = Lorentz(name = 'VVVV4',
                spins = [3, 3, 3, 3],
                structure = 'Metric(1,3)*Metric(2,4) - Metric(1,2)*Metric(3,4)')


VVVV5 = Lorentz(name = 'VVVV5',
                spins = [3, 3, 3, 3],
                structure = 'Metric(1,4)*Metric(2,3) - (Metric(1,3)*Metric(2,4))/2. - (Metric(1,2)*Metric(3,4))/2.')


FFS1__1 = Lorentz(name = 'FFS1__1',
                  spins = [2, 2, 1],
                  structure = 'P(-1,3)*Gamma(-1,2,-2)*ProjM(-2,1)')


FFSS1 = Lorentz(name = 'FFSS1',
                spins = [2, 2, 1, 1],
                structure = 'FFCT2((P(-3,3)+P(-3,4))*(P(-3,3)+P(-3,4))) *(P(-1,3)*Gamma(-1,2,-2)*ProjM(-2,1) - P(-1,4)*Gamma(-1,2,-2)*ProjM(-2,1))')


FFFF1 = Lorentz(name = 'FFFF1',
                spins = [2, 2, 2, 2],
                structure = 'Gamma(-1,2,-2)*Gamma(-1,4,-3)*ProjM(-3,3)*ProjM(-2,1)')


FFSSS1 = Lorentz(name = 'FFSSS1',
                 spins = [2, 2, 1, 1, 1],
                 structure = 'P(-1,3)*Gamma(-1,2,-2)*ProjM(-2,1) - P(-1,4)*Gamma(-1,2,-2)*ProjM(-2,1)')


FFSSS2 = Lorentz(name = 'FFSSS2',
                 spins = [2, 2, 1, 1, 1],
                 structure = 'FFCT3((P(-3,3)+P(-3,4)+P(-3,5))*(P(-3,3)+P(-3,4)+P(-3,5))) *FFCT3F1((P(-4,3)+P(-4,5))*(P(-4,3)+P(-4,5))) *(P(-1,3)*Gamma(-1,2,-2)*ProjM(-2,1) - P(-1,5)*Gamma(-1,2,-2)*ProjM(-2,1))')


FFSSS3 = Lorentz(name = 'FFSSS3',
                 spins = [2, 2, 1, 1, 1],
                 structure = 'FFCT3((P(-3,3)+P(-3,4)+P(-3,5))*(P(-3,3)+P(-3,4)+P(-3,5))) *FFCT3F1((P(-4,4)+P(-4,5))*(P(-4,4)+P(-4,5))) *(P(-1,4)*Gamma(-1,2,-2)*ProjM(-2,1) - P(-1,5)*Gamma(-1,2,-2)*ProjM(-2,1))')


FFSSS4 = Lorentz(name = 'FFSSS4',
                 spins = [2, 2, 1, 1, 1],
                 structure = '0.5 *FFCT3((P(-3,3)+P(-3,4)+P(-3,5))*(P(-3,3)+P(-3,4)+P(-3,5))) *( FFCT3F1((P(-4,3)+P(-4,5))*(P(-4,3)+P(-4,5)))*((P(-5,3)+P(-5,4)+P(-5,5))*(P(-5,3)-P(-5,5)))/((P(-6,3)+P(-6,4)+P(-6,5))*(P(-6,3)+P(-6,4)+P(-6,5))) + FFCT3F1((P(-7,4)+P(-7,5))*(P(-7,4)+P(-7,5)))*((P(-8,3)+P(-8,4)+P(-8,5))*(P(-8,4)-P(-8,5)))/((P(-9,3)+P(-9,4)+P(-9,5))*(P(-9,3)+P(-9,4)+P(-9,5))) ) *(P(-1,3)*Gamma(-1,2,-2)*ProjM(-2,1) + P(-1,4)*Gamma(-1,2,-2)*ProjM(-2,1) + P(-1,5)*Gamma(-1,2,-2)*ProjM(-2,1))')


FFSSS5 = Lorentz(name = 'FFSSS5',
                 spins = [2, 2, 1, 1, 1],
                 structure = 'FFCT3((P(-3,3)+P(-3,4)+P(-3,5))*(P(-3,3)+P(-3,4)+P(-3,5))) *FFCT3F0((P(-4,3)+P(-4,5))*(P(-4,3)+P(-4,5))) *(P(-1,3)*Gamma(-1,2,-2)*ProjM(-2,1) - P(-1,5)*Gamma(-1,2,-2)*ProjM(-2,1))')


FFSSS6 = Lorentz(name = 'FFSSS6',
                 spins = [2, 2, 1, 1, 1],
                 structure = 'FFCT3((P(-3,3)+P(-3,4)+P(-3,5))*(P(-3,3)+P(-3,4)+P(-3,5))) *FFCT3F0((P(-4,4)+P(-4,5))*(P(-4,4)+P(-4,5))) *(P(-1,4)*Gamma(-1,2,-2)*ProjM(-2,1) - P(-1,5)*Gamma(-1,2,-2)*ProjM(-2,1))')


FFSSS7 = Lorentz(name = 'FFSSS7',
                 spins = [2, 2, 1, 1, 1],
                 structure = '0.5 *FFCT3((P(-3,3)+P(-3,4)+P(-3,5))*(P(-3,3)+P(-3,4)+P(-3,5))) *( FFCT3F0((P(-4,3)+P(-4,5))*(P(-4,3)+P(-4,5)))*((P(-5,3)+P(-5,4)+P(-5,5))*(P(-5,3)-P(-5,5)))/((P(-6,3)+P(-6,4)+P(-6,5))*(P(-6,3)+P(-6,4)+P(-6,5))) + FFCT3F0((P(-7,4)+P(-7,5))*(P(-7,4)+P(-7,5)))*((P(-8,3)+P(-8,4)+P(-8,5))*(P(-8,4)-P(-8,5)))/((P(-9,3)+P(-9,4)+P(-9,5))*(P(-9,3)+P(-9,4)+P(-9,5))) ) *(P(-1,3)*Gamma(-1,2,-2)*ProjM(-2,1) + P(-1,4)*Gamma(-1,2,-2)*ProjM(-2,1) + P(-1,5)*Gamma(-1,2,-2)*ProjM(-2,1))')

