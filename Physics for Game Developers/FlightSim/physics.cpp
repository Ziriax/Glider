//----------------------------------------------------------------------------------------------------//
/*
	PHYSICS FOR GAME DEVELOPERS
	
	CHAPTERS 7 & 15 EXAMPLE PROGRAM

	NAME:		FlightSim
	PURPOSE:	To demonstrate 3D rigid body flight simulation
	BY:			David Bourg
	DATE:		07/24/00
	COPYRIGHT:	Copyright 2000 by David Bourg
*/
//----------------------------------------------------------------------------------------------------//
#include "Physics.h"
#include <cstdio>
#include <iostream>

using namespace std;

//------------------------------------------------------------------------//
// Global variables
//------------------------------------------------------------------------//



void main()
{
	const float dt = 0.1;

	Glider glider;

	cout << "pos: " << glider.Airplane.vPosition << " rot: " << glider.Airplane.qOrientation << endl;

	glider.StepSimulation(dt);

	cout << "pos: " << glider.Airplane.vPosition << " rot: " << glider.Airplane.qOrientation << endl;

	printf("Hello\n");
	getchar();
}