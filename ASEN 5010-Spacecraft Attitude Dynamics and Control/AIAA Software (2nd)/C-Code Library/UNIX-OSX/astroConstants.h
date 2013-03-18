/*
 *  astroConstants.h
 *  OrbitalMotion
 *
 *  Created by Hanspeter Schaub on 6/19/05.
 *
 *	Defines common astronomical constants using in the 
 *          Keplerian 2-body problem.  
 *
 */

#ifndef ASTROCONSTANTS
#define ASTROCONSTANTS


/*
  universial Gravitational Constant
  units are in km^3/s^2/kg
*/
#define G_UNIVERSIAL 6.67259e-20          

/*
  astronomical unit in units of kilometers
*/
#define AU          149597870.691     



/*
  common conversions
*/
#ifndef M_PI
#define M_PI	    3.141592653589793
#endif

#ifndef D2R
#define D2R         M_PI/180.
#endif
#ifndef R2D
#define R2D         180./M_PI
#endif



/*
  Gravitational Constants mu = G*m, where m is the planet of the attracting
  planet.  All units are km^3/s^2. 
  values are obtained from http://ssd.jpl.nasa.gov/astro_constants.html
*/
#define MU_SUN       1.32712440018e11  
#define MU_MERCURY   2.2032e4          
#define MU_VENUS     3.2485859e5       
#define MU_EARTH     3.986004418e5     
#define MU_MOON      4.9027988e3       
#define MU_MARS      4.28283e4         
#define MU_JUPITER   1.2671277e8       
#define MU_SATURN    3.79406e7         
#define MU_URANUS    5.79455e6         
#define MU_NEPTUNE   6.83653e6         
#define MU_PLUTO     983.                           




        


/*
  planet information for major solar system bodies. Units are in km.
  data taken from http://nssdc.gsfc.nasa.gov/planetary/planets.html
*/
/* Sun:  */
#define REQ_SUN      695000.           

/* Mercury: */
#define REQ_MERCURY    2439.7          
#define J2_MERCURY       60.0e-6
#define SMA_MERCURY       0.38709893*AU
#define I_MERCURY         7.00487*D2R
#define E_MERCURY         0.20563069

/* Venus:   */
#define REQ_VENUS      6051.8          
#define J2_VENUS          4.458e-6
#define SMA_VENUS         0.72333199*AU
#define I_VENUS           3.39471*D2R
#define E_VENUS           0.00677323

/* Earth:  */
#define REQ_EARTH      6378.14         
#define SMA_EARTH         1.00000011*AU
#define I_EARTH           0.00005*D2R
#define E_EARTH           0.01671022

/* Moon:  */
#define REQ_MOON       1737.4          
#define J2_MOON         202.7e-6
#define SMA_MOON          0.3844e6
#define E_MOON            0.0549

/* Mars:  */
#define REQ_MARS       3397.2          
#define J2_MARS        1960.45e-6
#define SMA_MARS          1.52366231*AU
#define I_MARS            1.85061*D2R
#define E_MARS            0.09341233

/* Jupiter:  */
#define REQ_JUPITER   71492.           
#define J2_JUPITER    14736.e-6
#define SMA_JUPITER       5.20336301*AU
#define I_JUPITER         1.30530*D2R
#define E_JUPITER         0.04839266

/* Saturn:  */
#define REQ_SATURN    60268.           
#define J2_SATURN     16298.e-6
#define SMA_SATURN        9.53707032*AU
#define I_SATURN          2.48446*D2R
#define E_SATURN          0.05415060

/* Uranus:   */
#define REQ_URANUS    25559.           
#define J2_URANUS      3343.43e-6
#define SMA_URANUS       19.19126393*AU
#define I_URANUS          0.76986*D2R
#define E_URANUS          0.04716771

/* Neptune:   */
#define REQ_NEPTUNE   24746.           
#define J2_NEPTUNE     3411.e-6
#define SMA_NEPTUNE      30.06896348*AU
#define I_NEPTUNE         1.76917*D2R
#define E_NEPTUNE         0.00858587

/* Pluto:  */
#define REQ_PLUTO      1137.           
#define SMA_PLUTO        39.48168677*AU
#define I_PLUTO          17.14175*D2R
#define E_PLUTO           0.24880766







/* 
  Zonal gravitational harmonics Ji for the Earth
*/
#define J2  1082.63e-6
#define J3    -2.52e-6
#define J4    -1.61e-6
#define J5    -0.15e-6
#define J6     0.57e-6






#endif      