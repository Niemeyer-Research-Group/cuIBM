/***************************************************************************//**
 * \file body.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Definition of the class \c body.
 */


#pragma once
#include <cusp/array1d.h>

/**
 * \class body
 * \brief Contains information about an immersed boundary in the flow.
 *
 * The number of boundary points and the reference x- and y- coordinates of the
 * boundary points are provided. The motion of the body is also described.
 * For uniform motion, the translational velocity needs to be specified.
 * For uniform rotation, the angular velocity has to be specified.
 * You can also provide options for oscillations in the x and y directions,
 * as well as pitching (rotational) oscillations.
 * When the body rotates, the center of rotation must be provided too.
 *
 */
class body
{
  public:
    int  numPoints; ///< number of boundary points
	
    cusp::array1d<double, cusp::host_memory>
    	X,         ///< reference x-coordinate of boundary points
		Y;         ///< reference y-coordinate of boundary points
         
	double X0[2];     ///< reference center of rotation
	     
	double Xc0[2],    ///< initial position of center of rotation
	     Theta0;    ///< initial angle of attack
    
    double Xc[2],     ///< actual center of rotation (x- and y-coordinates)
	     Theta,     ///< actual angle of attack (counterclockwise is positive)
	     vel[2],    ///< translational velocity (x- and y- components)
	     angVel;    ///< angular velocity (counterlockwise is positive)

	bool moving[2]; ///< flag to indicate if the body is moving (translating or rotating)

	double velocity[2], ///< uniform translational velocity (x- and y-components)
	     omega;       ///< uniform rotational velocity

	double xOscillation[3],     ///< amplitude, angular frequency and phase difference of oscillation in the x-direction
         yOscillation[3],     ///< amplitude, angular frequency and phase difference of oscillation in the y-direction
         pitchOscillation[3]; ///< amplitude, angular frequency and phase difference of pitch oscillation

	double	xCoefficient, ///< used in oscillating cylinder movement
			yCoefficient,
			uCoefficient,
			vCoefficient,
			xfrequency,
			yfrequency,
			xPhase,
			yPhase,
			uPhase,
			vPhase;
	// Note that the angular amplitude & phase difference are stored in radians, and angular frequncy is radians/time
	// But the inputs in bodies.yaml are provided in degrees, and degrees/time.
    
	/**
	 * \brief Updates the position of the center of rotation of the body,
	 *        as well as its the translational and rotational velocities.
	 *
	 * \param Time the time
	 */
	void update(double Time)
	{
		// if the body is translating
		if(moving[0])
		{
			Xc[0] = Xc0[0] + velocity[0]*Time + xOscillation[0]*sin(xOscillation[1]*Time + xOscillation[2]);
			Xc[1] = Xc0[1] + velocity[1]*Time + yOscillation[0]*sin(yOscillation[1]*Time + yOscillation[2]);
			vel[0]= velocity[0] + xOscillation[0]*xOscillation[1]*cos(xOscillation[1]*Time + xOscillation[2]);
			vel[1]= velocity[1] + yOscillation[0]*yOscillation[1]*cos(yOscillation[1]*Time + yOscillation[2]);
		}
		// if the body is rotating
		if(moving[1])
		{
			Theta = Theta0 + omega*Time + pitchOscillation[0]*sin(pitchOscillation[1]*Time + pitchOscillation[2]);
			angVel = omega + pitchOscillation[0]*pitchOscillation[1]*cos(pitchOscillation[1]*Time + pitchOscillation[2]);
		}
	}

	/**
	 * \brief Resets the position and motion of the body.
	 *
	 * The body center is now at (0,0) and is not moving anymore.
	 *
	 */
	void reset()
	{
		X0[0]  = X0[1]  = 0.0;
		Xc0[0] = Xc0[1] = 0.0;
		Theta0 = 0.0;
		Xc[0]  = Xc[1]  = 0.0;
		Theta  = 0.0;
		vel[0] = vel[1] = 0.0;
		angVel = 0.0;
		moving[0] = false;
		moving[1] = false;
		velocity[0] = 0.0;
		velocity[1] = 0.0;
		omega = 0.0;
		xOscillation[0] = xOscillation[1] = xOscillation[2] = 0.0;
		yOscillation[0] = yOscillation[1] = yOscillation[2] = 0.0;
		pitchOscillation[0] = pitchOscillation[1] = pitchOscillation[2] = 0.0;
	}
};
