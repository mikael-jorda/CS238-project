#include <iostream>
#include <string>
#include <thread>
#include <math.h>

// #include "taco_library/include/taco/optimizer/QP.h"
// #include "taco_library/include/taco/optimizer/QuadProgRT.h"

#include "model/ModelInterface.h"
#include "graphics/ChaiGraphics.h"
#include "simulation/Sai2Simulation.h"
#include <dynamics3d.h>

#include <fstream>

#include "force_sensor/ForceSensorSim.h"
#include "timer/LoopTimer.h"

#include <GLFW/glfw3.h> //must be loaded after loading opengl/glew as part of graphicsinterface

#include <signal.h>
bool fSimulationRunning = false;
void sighandler(int){fSimulationRunning = false;}

using namespace std;

const string world_file = "../resources/02-particle_filter/world.urdf";
const string robot_file = "../resources/02-particle_filter/4pbot_fixed.urdf";
const string robot_name = "4PBOT";

const string camera_name = "camera_fixed";

const string datafile_estimate = "../../02-particle_filter/data_logging/estimates.txt";
const string datafile_real = "../../02-particle_filter/data_logging/real.txt";
ofstream real_file, estimate_file;

Eigen::Vector3d robot_origin = Eigen::Vector3d(0.0, -1.5, 0.050001);

// filtering the force signals
// bool filter_forces = false;
bool filter_forces = true;
// filter eq : y[n] = (x[n-2] + 2x[n-1] + x[n])/gain + a1*y[n-2] + a2*y[n-1]
double coeff_controller_filter[2] = {-0.87521, 1.86689}; // assumes control freq = 1kHz, hard coded butterworth filter order 2 15Hz cutoff
double coeff_simulation_filter[2] = {-0.93553, 1.93338}; // assumes sim freq = 2kHz, hard coded butterworth filter order 2 15Hz cutoff
double controller_filter_gain = 4.806381793e+02;
double simulation_filter_gain = 1.861608837e+03;
Eigen::Vector2d Fepp, Fep, Feppf, Fepf, Fef;
Eigen::Vector2d Frpp, Frp, Frppf, Frpf, Frf;
Eigen::Vector2d Pepp, Pep, Peppf, Pepf, Pef;
Eigen::Vector2d Prpp, Prp, Prppf, Prpf, Prf;

bool contact_detected = false;

int sampleSize = 50;
int extraSamples = 0;
Eigen::MatrixXd particleSet = Eigen::MatrixXd::Zero(sampleSize+extraSamples,3);

// simulation loop
void control(Model::ModelInterface* robot, Simulation::Sai2Simulation* sim);
void simulation(Model::ModelInterface* robot, Simulation::Sai2Simulation* sim);

// initialize window manager
GLFWwindow* glfwInitialize();

// callback to print glfw errors
void glfwError(int error, const char* description);

// callback when a key is pressed
void keySelect(GLFWwindow* window, int key, int scancode, int action, int mods);

// callback when a mouse button is pressed
void mouseClick(GLFWwindow* window, int button, int action, int mods);

// flags for scene camera movement
bool fTransXp = false;
bool fTransXn = false;
bool fTransYp = false;
bool fTransYn = false;
bool fTransZp = false;
bool fTransZn = false;
bool fRotPanTilt = false;

/**
 * resample()
 * ----------------------------------------------------
 * multimodal resampling
 */
Eigen::VectorXd resample(Eigen::VectorXd normalizedWeights) 
{
	double numParticles = normalizedWeights.size();
	Eigen::VectorXd idx = Eigen::VectorXd::Zero(numParticles);
	Eigen::VectorXd cumSumWeights = Eigen::VectorXd::Zero(numParticles);
	int k = 0;

	cumSumWeights(0) = normalizedWeights(0); 
	for (int j=1; j<numParticles; j++)
	{
		cumSumWeights(j) = cumSumWeights(j-1) + normalizedWeights(j);
	}
	// cout << "cumulative sum: " << cumSumWeights << endl;

	for (int i=0; i<numParticles; i++)
	{
		// generate random number between 0 and 1
		double rand4resamp = ((double) rand() / (RAND_MAX));

		k = 0;
		while(cumSumWeights(k)<rand4resamp && k < numParticles-1)
		{
			k++;
		}

		idx(i) = k;
	}

	return idx;
}

/**
 * motionModel()
 * -------------------------------------------------------
 * redistribute particles according to normal distribution
 */
Eigen::MatrixXd motionModel(Eigen::MatrixXd particleSet) 
{
	double numParticles = particleSet.rows();
	Eigen::MatrixXd movedParticleSet = Eigen::MatrixXd::Zero(numParticles,3);
	Eigen::Vector3d particleLocation;
	float normalSample;

	for (int i=0; i<numParticles; i++)
	{
		// random device class instance, source of 'true' randomness for initializing random seed
	    random_device rd; 

	    // Mersenne twister PRNG, initialized with seed from random device instance
	    mt19937 gen(rd()); 

	    //particleLocation
	    particleLocation = particleSet.row(i);

	    // instance of class normal_distribution with specific mean and stddev
        normal_distribution<float> d(particleLocation(2), 0.005); 

        // get random number with normal distribution using gen as random source
        normalSample = d(gen); 

        // save location
        movedParticleSet.row(i) << 0, particleLocation(1), normalSample;
	}

	return movedParticleSet;
}

int main (int argc, char** argv) {
	cout << "Loading URDF world model file: " << world_file << endl;

	// set up signal handler
	signal(SIGABRT, &sighandler);
	signal(SIGTERM, &sighandler);
	signal(SIGINT, &sighandler);

	// load graphics scene
	auto graphics = new Graphics::ChaiGraphics(world_file, Graphics::urdf, false);
	Eigen::Vector3d camera_pos, camera_lookat, camera_vertical;
	graphics->getCameraPose(camera_name, camera_pos, camera_vertical, camera_lookat);

	vector<chai3d::cShapeSphere*> graphic_particles;
	for(int i=0; i < sampleSize + extraSamples; i++)
	{
		graphic_particles.push_back(new chai3d::cShapeSphere(0.05));
		graphic_particles[i]->m_material->setColorf(1.0,0.0,0.0);
		graphics->_world->addChild(graphic_particles[i]);
	}
	camera_lookat << 0.0, 0.0, 1.0;
	auto sphere_rpos = new chai3d::cShapeSphere(0.08);
	auto sphere_epos = new chai3d::cShapeSphere(0.08);
	auto line_rforce = new chai3d::cShapeLine();
	auto line_eforce = new chai3d::cShapeLine();
	graphics->_world->addChild(sphere_epos);
	graphics->_world->addChild(sphere_rpos);
	graphics->_world->addChild(line_eforce);
	graphics->_world->addChild(line_rforce);
	// change position
	// sphere->setLocalPos(Eigen::Vector3d(0.0, 0.0, 1.0));
	// color
	sphere_epos->m_material->setColorf(0.6, 0.5, 0.2);
	sphere_rpos->m_material->setColorf(0.2, 0.5, 0.3);

	line_rforce->m_colorPointA = chai3d::cColorf(0.2, 0.5, 0.3);
	line_rforce->m_colorPointB = chai3d::cColorf(0.2, 0.5, 0.3);
	line_eforce->m_colorPointA = chai3d::cColorf(0.6, 0.5, 0.2);
	line_eforce->m_colorPointB = chai3d::cColorf(0.6, 0.5, 0.2);
	line_rforce->setLineWidth(4.0);
	line_eforce->setLineWidth(4.0);

	sphere_rpos->setShowEnabled(false);
	sphere_epos->setShowEnabled(false);
	line_rforce->setShowEnabled(false);
	line_eforce->setShowEnabled(false);

	// load robots
	auto robot = new Model::ModelInterface(robot_file, Model::rbdl, Model::urdf, false);

	// load simulation world
	auto sim = new Simulation::Sai2Simulation(world_file, Simulation::urdf, false);
	sim->setCollisionRestitution(0);
	sim->setCoeffFrictionStatic(0.6);

	// set initial condition
	robot->_q << -75.0/180.0*M_PI,
				  25.0/180.0*M_PI,
				 -60.0/180.0*M_PI,
				  10.0/180.0*M_PI;
	sim->setJointPositions(robot_name, robot->_q);
	robot->updateModel();

	// initialize GLFW window
	GLFWwindow* window = glfwInitialize();

	double last_cursorx, last_cursory;

    // set callbacks
	glfwSetKeyCallback(window, keySelect);
	glfwSetMouseButtonCallback(window, mouseClick);

	// open files for data logging
	real_file.open(datafile_real);
	estimate_file.open(datafile_estimate);
	estimate_file << "time \t force \t position\n";
	real_file << "time \t force \t position\n";

	// initialize filter variables
	Fepp.setZero();
	Fep.setZero();
	Feppf.setZero();
	Fepf.setZero();
	Fef.setZero();
	Frpp.setZero();
	Frp.setZero();
	Frppf.setZero();
	Frpf.setZero();
	Frf.setZero();
	Pepp.setZero();
	Pep.setZero();
	Peppf.setZero();
	Pepf.setZero();
	Pef.setZero();
	Prpp.setZero();
	Prp.setZero();
	Prppf.setZero();
	Prpf.setZero();
	Prf.setZero();

	// start the simulation thread first
	fSimulationRunning = true;
	thread sim_thread(simulation, robot, sim);

	// next start the control thread
	thread ctrl_thread(control, robot, sim);
	
    // while window is open:
    while (!glfwWindowShouldClose(window)) {
		// update kinematic models
		// robot->updateModel();

		Eigen::Vector3d link_origin = Eigen::Vector3d::Zero();
		Eigen::Matrix3d Rlink = Eigen::Matrix3d::Identity();
		robot->position(link_origin, "link3", Eigen::Vector3d::Zero(3));
		robot->rotation(Rlink, "link3");

    	for(int i=0; i< sampleSize+extraSamples; i++)
    	{
    		Eigen::Vector3d tmp = particleSet.row(i);
			graphic_particles[i]->setLocalPos(robot_origin + link_origin + Rlink*tmp);
    	}

		// update contact position and force
		if(contact_detected)
		{
			sphere_rpos->setShowEnabled(true);
			sphere_epos->setShowEnabled(true);
			line_rforce->setShowEnabled(true);
			line_eforce->setShowEnabled(true);
			sphere_epos->setLocalPos(Eigen::Vector3d(0.02,Pef(0),Pef(1)));
			sphere_rpos->setLocalPos(Eigen::Vector3d(0.02,Prf(0),Prf(1)));
			line_eforce->m_pointA = chai3d::cVector3d(Eigen::Vector3d(0.02,Pef(0),Pef(1)));
			line_eforce->m_pointB = chai3d::cVector3d(Eigen::Vector3d(0.02,Pef(0)+0.05*Fef(0),Pef(1)+0.05*Fef(1)));
			line_rforce->m_pointA = chai3d::cVector3d(Eigen::Vector3d(0.02,Prf(0),Prf(1)));
			line_rforce->m_pointB = chai3d::cVector3d(Eigen::Vector3d(0.02,Prf(0)+0.05*Frf(0),Prf(1)+0.05*Frf(1)));
		}
		else
		{
			sphere_rpos->setShowEnabled(false);
			sphere_epos->setShowEnabled(false);
			line_rforce->setShowEnabled(false);
			line_eforce->setShowEnabled(false);
		}



		// update graphics. this automatically waits for the correct amount of time
		int width, height;
		glfwGetFramebufferSize(window, &width, &height);
		graphics->updateGraphics(robot_name, robot);
		graphics->render(camera_name, width, height);
		glfwSwapBuffers(window);
		glFinish();

	    // poll for events
	    glfwPollEvents();

		// move scene camera as required
    	// graphics->getCameraPose(camera_name, camera_pos, camera_vertical, camera_lookat);
    	Eigen::Vector3d cam_depth_axis;
    	cam_depth_axis = camera_lookat - camera_pos;
    	cam_depth_axis.normalize();
    	Eigen::Vector3d cam_up_axis;
    	// cam_up_axis = camera_vertical;
    	// cam_up_axis.normalize();
    	cam_up_axis << 0.0, 0.0, 1.0; //TODO: there might be a better way to do this
	    Eigen::Vector3d cam_roll_axis = (camera_lookat - camera_pos).cross(cam_up_axis);
    	cam_roll_axis.normalize();
    	Eigen::Vector3d cam_lookat_axis = camera_lookat;
    	cam_lookat_axis.normalize();
    	if (fTransXp) {
	    	camera_pos = camera_pos + 0.05*cam_roll_axis;
	    	camera_lookat = camera_lookat + 0.05*cam_roll_axis;
	    }
	    if (fTransXn) {
	    	camera_pos = camera_pos - 0.05*cam_roll_axis;
	    	camera_lookat = camera_lookat - 0.05*cam_roll_axis;
	    }
	    if (fTransYp) {
	    	// camera_pos = camera_pos + 0.05*cam_lookat_axis;
	    	camera_pos = camera_pos + 0.05*cam_up_axis;
	    	camera_lookat = camera_lookat + 0.05*cam_up_axis;
	    }
	    if (fTransYn) {
	    	// camera_pos = camera_pos - 0.05*cam_lookat_axis;
	    	camera_pos = camera_pos - 0.05*cam_up_axis;
	    	camera_lookat = camera_lookat - 0.05*cam_up_axis;
	    }
	    if (fTransZp) {
	    	camera_pos = camera_pos + 0.1*cam_depth_axis;
	    	camera_lookat = camera_lookat + 0.1*cam_depth_axis;
	    }	    
	    if (fTransZn) {
	    	camera_pos = camera_pos - 0.1*cam_depth_axis;
	    	camera_lookat = camera_lookat - 0.1*cam_depth_axis;
	    }
	    if (fRotPanTilt) {
	    	// get current cursor position
	    	double cursorx, cursory;
			glfwGetCursorPos(window, &cursorx, &cursory);
			//TODO: might need to re-scale from screen units to physical units
			double compass = 0.006*(cursorx - last_cursorx);
			double azimuth = 0.006*(cursory - last_cursory);
			double radius = (camera_pos - camera_lookat).norm();
			Eigen::Matrix3d m_tilt; m_tilt = Eigen::AngleAxisd(azimuth, -cam_roll_axis);
			camera_pos = camera_lookat + m_tilt*(camera_pos - camera_lookat);
			Eigen::Matrix3d m_pan; m_pan = Eigen::AngleAxisd(compass, -cam_up_axis);
			camera_pos = camera_lookat + m_pan*(camera_pos - camera_lookat);
	    }
	    graphics->setCameraPose(camera_name, camera_pos, cam_up_axis, camera_lookat);
	    glfwGetCursorPos(window, &last_cursorx, &last_cursory);
	}

	// stop simulation
	fSimulationRunning = false;
	sim_thread.join();
	ctrl_thread.join();

	// close files
	real_file.close();
	estimate_file.close();

    // destroy context
    glfwDestroyWindow(window);

    // terminate
    glfwTerminate();

	return 0;
}

//------------------------------------------------------------------------------
void control(Model::ModelInterface* robot, Simulation::Sai2Simulation* sim) {
	robot->updateModel();

	int dof = robot->dof();
	Eigen::VectorXd command_torques = Eigen::VectorXd::Zero(dof);
	Eigen::VectorXd gravity_compensation;
	Eigen::VectorXd nonlinear_terms;

	// joint task
	Eigen::VectorXd joint_task_torques = Eigen::VectorXd::Zero(dof);
	double kp_joint = 0.0;
	double kv_joint = 10.0;

	// position task
	const string link_name = "link4";
	const Eigen::Vector3d pos_in_link = Eigen::Vector3d(0.0, 0.0, 1.0);
	Eigen::MatrixXd J3d, Jv2d, Jw2d, Lambda, N, Jbar;
	Eigen::MatrixXd Jr = Eigen::MatrixXd::Zero(6, dof);
	J3d = Eigen::MatrixXd::Zero(6, dof);
	Jw2d = Eigen::MatrixXd::Zero(1, dof);
	Jv2d = Eigen::MatrixXd::Zero(2, dof);
	Lambda = Eigen::MatrixXd(2,2);
	Jbar = Eigen::MatrixXd(dof,2);
	N = Eigen::MatrixXd(dof,dof);

	Eigen::Vector2d pos_task_force;
	Eigen::VectorXd pos_task_torques;

	double kp_pos = 100.0;
	double kv_pos = 20.0;

	Eigen::Vector3d pos3d;
	Eigen::Vector2d x, xdot;
	Eigen::Vector3d x_init;
	robot->position(x_init, link_name, pos_in_link);
	Eigen::Vector2d xd = x_init.tail(2);

	// initialize momentum observer
	Eigen::VectorXd r = Eigen::VectorXd::Zero(dof); 
	Eigen::VectorXd r_int = Eigen::VectorXd::Zero(dof);
	Eigen::VectorXd p = Eigen::VectorXd::Zero(dof);
	Eigen::MatrixXd K0 = Eigen::MatrixXd::Zero(dof,dof);
	for(int i=0; i<dof; i++)
	{
		K0(i,i) = 100.0;
	}
	Eigen::VectorXd beta = Eigen::VectorXd::Zero(dof);

	// collision identification
	Eigen::VectorXd mu = Eigen::VectorXd::Zero(dof);
	double eps_mu = 0.7;
	Eigen::Vector3d Fi;
	int collision_index = 0;
	// Eigen::MatrixXd Ji = Eigen::MatrixXd(3,dof);
	Eigen::Vector2d estimated_Fc;
	Eigen::Vector2d estimated_rc;
	Eigen::Matrix3d Ri3d;
	// Eigen::Matrix2d Ri;
	Eigen::Vector3d pos_i3d;
	// Eigen::Vector2d pos_i;

	// initialize particle filter variables

	Eigen::MatrixXd particleSetInit = Eigen::MatrixXd::Zero(sampleSize+extraSamples,3);
	Eigen::MatrixXd particleForce = Eigen::MatrixXd::Zero(sampleSize+extraSamples,2);
	// Eigen::MatrixXd particleInfo = Eigen::MatrixXd::Zero(sampleSize+extraSamples,4); //contains y,z,weight,force
	Eigen::MatrixXd resampledParticleSet = Eigen::MatrixXd::Zero(sampleSize+extraSamples,3);
	Eigen::VectorXd weights = Eigen::VectorXd::Zero(sampleSize+extraSamples);
	Eigen::VectorXd resampledParticleIdx = Eigen::VectorXd::Zero(sampleSize+extraSamples);
	Eigen::MatrixXd Jparticle, A, B;
	Eigen::Vector3d Fopt3d;
	Eigen::Vector2d Fopt2d;

	// sample particles on link surface
	Eigen::Vector3d linkDepthVec;
	Eigen::VectorXd error = Eigen::VectorXd::Zero(4);
	linkDepthVec << 0,0.05,0;
	Eigen::Vector3d particleLocation;

	// initial uniform sampling
	for (int i=0; i<sampleSize; i++)
	{	
		int sign;
		if (i % 2 == 0)
		{
			sign = 1;
		}
		else
		{
			sign = -1;
		}
		particleLocation << sign*linkDepthVec + Eigen::Vector3d(0,0,(double)i/(sampleSize)); 
		// cout << "Particle : " << particleLocation << endl;
		particleSetInit.row(i) = particleLocation;
	}
	// particleSet = particleSetInit;

	// cout << particleSetInit << endl;

	// create a loop timer
	double dt = 0.001;
	double control_freq = 1/dt;
	LoopTimer timer;
	timer.setLoopFrequency(control_freq);   // 1 KHz
	double last_time = timer.elapsedTime(); //secs
	bool fTimerDidSleep = true;
	timer.initializeTimer(1000000); // 1 ms pause before starting loop

	unsigned long long controller_counter = 0;

	// Eigen::Vector3d sensed_force = Eigen::Vector3d::Zero();

	bool gpjs = true;
	// gpjs = false;

	while (fSimulationRunning) { //automatically set to false when simulation is quit
		fTimerDidSleep = timer.waitForNextLoop();

		// update time
		double curr_time = timer.elapsedTime();
		double loop_dt = curr_time - last_time;

		// read joint positions, velocities, update model
		sim->getJointPositions(robot_name, robot->_q);
		sim->getJointVelocities(robot_name, robot->_dq);
		robot->updateModel();
		robot->gravityVector(gravity_compensation);

		// update momentum residual
		p = robot->_M * robot->_dq;
		beta = gravity_compensation;   // negelecting coriolis for now
		r_int += (command_torques - beta + r)*dt;
		r = K0*(p - r_int);

		// collision identification and isolation
		mu.setZero();
		// pos_i.setZero();
		collision_index = 0;
		for(int i=0; i<dof; i++)
		{
			if(r(i) > eps_mu)
			{
				mu(i) = r(i);
				collision_index = i+1; //link where the collision happened
			}
		}
		switch(collision_index)
		{
			case 3 :
			contact_detected = true;
			// robot->J_0(J3d, "link3", Eigen::Vector3d::Zero());
			robot->rotation(Ri3d, "link3");
			robot->position(pos_i3d, "link3", Eigen::Vector3d::Zero());
			//string currentLink("link3");
			break;

			case 4 :
			contact_detected = true;
			// robot->J_0(J3d, "link4", Eigen::Vector3d::Zero());
			robot->rotation(Ri3d, "link4");
			robot->position(pos_i3d, "link4", Eigen::Vector3d::Zero());
			//string currentLink("link4");
			break;

			default :
			contact_detected = false;
			J3d.setZero();

			Fepp.setZero();
			Fep.setZero();
			Feppf.setZero();
			Fepf.setZero();
			Fef.setZero();
			Frpp.setZero();
			Frp.setZero();
			Frppf.setZero();
			Frpf.setZero();
			Frf.setZero();
			Pepp.setZero();
			Pep.setZero();
			Peppf.setZero();
			Pepf.setZero();
			Pef.setZero();
			Prpp.setZero();
			Prp.setZero();
			Prppf.setZero();
			Prpf.setZero();
			Prf.setZero();
		}
	
		// cout << "line 447 " << endl;

		// if(collision_index > 2 and collision_index <= dof)
		// {
		// 	// define matrices
		// 	Ri = Ri3d.block(1,1,2,2);
		// 	Ji = J3d.block(1,0,3,dof);
		// 	// Ji.block(0,0,2,dof) = J3d.block(1,0,2,dof);
		// 	// Ji.block(2,0,1,dof) = J3d.block(3,0,1,dof);

		// 	// find equivqlent force at the joint frame
		// 	Fi = (Ji.transpose()).colPivHouseholderQr().solve(mu);

		// 	// deduce applied force
		// 	estimated_Fc = Fi.head(2);
		// 	estimated_Fc = Ri.transpose()*estimated_Fc;

		// 	// find contact point via moment
		// 	double ry = 0;
		// 	if(estimated_Fc(0) > 0)
		// 	{
		// 		ry = -0.05;
		// 	}
		// 	else
		// 	{
		// 		ry = 0.05;
		// 	}
		// 	double rz = (ry*estimated_Fc(1) - Fi(2))/estimated_Fc(0);
		// 	estimated_rc << ry, rz;

		// 	// transform to base frame
		// 	pos_i = pos_i3d.tail(2);
		// 	estimated_rc = robot_origin.tail(2) + pos_i + Ri*estimated_rc;
		// 	estimated_Fc = Ri*estimated_Fc;
		// }



		// -------------------------------------------
		////////////////////////////// Particle Filter


		if(contact_detected)
		{

			// use initial uniform particle set or motion model for the next time step
			if (resampledParticleSet.isZero())
			{
				particleSet = particleSetInit;
			}
			else // move resampled particles according to normal distribution
			{
				particleSet = motionModel(resampledParticleSet);
			}

			// Add 20 extra random samples to the current particleSet
			for (int i=sampleSize+1; i<(sampleSize + extraSamples); i++)
			{	
				int sign;
				if (i % 2 == 0)
				{
					sign = 1;
				}
				else
				{
					sign = -1;
				}
				// particleLocation << Ri3d.transpose()*pos_i3d + sign*linkDepthVec + Eigen::Vector3d(0,0,i*1/(sampleSize/2));  // CHECK 
				particleLocation << sign*linkDepthVec + Eigen::Vector3d(0,0,(double)(i - sampleSize)/(extraSamples)); 
				particleSet.row(i) = particleLocation;
			}

			// compute importance weights
			for (int j=0; j<sampleSize + extraSamples; j++)
			{
				Eigen::Vector3d particleVec;
				particleVec = particleSet.row(j);

				if(collision_index == 3) //find the jacobian for the current particle location
				{
			 		robot->Jv(Jparticle, "link3", particleVec); 
				}
			 	else
			 	{
			 		robot->Jv(Jparticle, "link4", particleVec); 

			 	}

			 	//Minimize   FAF+BF
	            // A = Jparticle * Jparticle.transpose();
	            // B = 2 * r.transpose() * Jparticle.transpose();
	            // Fopt3d = -2 * A.inverse() * B.transpose();

				Fopt3d = (Jparticle.transpose()).colPivHouseholderQr().solve(r);
				// Fopt3d.setZero();
	            Fopt2d = Fopt3d.tail(2);
	            particleForce.row(j) = Fopt2d;

	            // Eigen::Vector3d Flocframe = Ri3d.transpose()*Fopt3d;
	            
	            //Subject to [0 -1 0]*F<=0
	            // if((Flocframe(1) > 0 && particleVec(1) < 0) || (Flocframe(1) < 0 && particleVec(1) > 0))
	            // {
	            	error = (r - Jparticle.transpose()*Fopt3d);  
	            	// error.setZero();  
				 	weights(j) = exp(-0.5 * error.norm()*error.norm()); 	
	            // }
	            // else
	            // {
		         	// weights(j) = 0; //discard particles that give an F in the wrong half plane   
	            // }	
			}

			// normalize weights
			weights = weights/ weights.sum();
			
			// resample
			resampledParticleIdx = resample(weights);

			for (int i=0; i<sampleSize + extraSamples; i++)
			{
				resampledParticleSet.row(i) = particleSet.row(resampledParticleIdx(i));
			}

			// select the average
			Eigen::Vector3d estimated_rc_3d = resampledParticleSet.colwise().mean();
			estimated_rc_3d = robot_origin + pos_i3d + Ri3d*estimated_rc_3d;
			estimated_rc = estimated_rc_3d.tail(2);
			// estimated_rc.setZero();
			estimated_Fc = particleForce.colwise().mean();
		}
		else
		{
			particleSet = particleSetInit;
			resampledParticleSet.setZero();
			estimated_Fc.setZero();
			estimated_rc.setZero();
		}

		// filter forces
		if(filter_forces)
		{
			Fef = (Fepp + 2*Fep + estimated_Fc)/controller_filter_gain + coeff_controller_filter[0]*Feppf + coeff_controller_filter[1]*Fepf;
			Fepp = Fep;
			Fep = estimated_Fc;
			Feppf = Fepf;
			Fepf = Fef;

			Pef = (Pepp + 2*Pep + estimated_rc)/controller_filter_gain + coeff_controller_filter[0]*Peppf + coeff_controller_filter[1]*Pepf;
			Pepp = Pep;
			Pep = estimated_rc;
			Peppf = Pepf;
			Pepf = Pef;
		}
		else
		{
			Fef = estimated_Fc;
			Pef = estimated_rc;
		}

		// -------------------------------------------
		// jacobian for control
		robot->J_0(J3d, link_name, pos_in_link);
		Jv2d = J3d.block(1,0,2,dof);
		robot->operationalSpaceMatrices(Lambda, Jbar, N, Jv2d);
		////////////////////////////// Compute joint torques
		double time = controller_counter/control_freq;


		//----------- pos task
		robot->position(pos3d, link_name, pos_in_link);
		x = pos3d.tail(2);
		xdot = Jv2d*robot->_dq;

		double mvt_freq = 0.3;
		xd(0) = x_init(1) + 0.15*sin(2*M_PI*mvt_freq* time);

		pos_task_force = Lambda*(-kp_pos*(x - xd) - kv_pos*xdot);
		if(!gpjs)
		{
			pos_task_force += Jbar.transpose() * gravity_compensation;
		}
		pos_task_torques = Jv2d.transpose() * pos_task_force;

		//----- Joint nullspace damping
		joint_task_torques = robot->_M*( - kv_joint*robot->_dq);

		//------ Final torques
		command_torques = pos_task_torques + N.transpose()*joint_task_torques;
		if(gpjs)
		{
			command_torques += gravity_compensation;
		}
		// command_torques.setZero();

		// -------------------------------------------
		sim->setJointTorques(robot_name, command_torques);
		

		// -------------------------------------------
		// if(controller_counter % 500 == 0)
		// {
		// 	cout << "Fc : " << Fef.transpose() << endl;
		// 	cout << "rc : " << estimated_rc.transpose() << endl;
		// 	cout << endl;
		// }
		// if(controller_counter == 3000)
		if(timer.elapsedTime() >= 1)
		{
			gpjs = false;
		}

		controller_counter++;

		// write estimates to file
		estimate_file << timer.elapsedTime() << '\t' << Fef.transpose() << '\t' << Pef.transpose() << endl;

		// -------------------------------------------
		// update last time
		last_time = curr_time;
	}

	double end_time = timer.elapsedTime();
    std::cout << "\n";
    std::cout << "Control Loop run time  : " << end_time << " seconds\n";
    std::cout << "Control Loop updates   : " << timer.elapsedCycles() << "\n";
    std::cout << "Control Loop frequency : " << timer.elapsedCycles()/end_time << "Hz\n";

}

//------------------------------------------------------------------------------
void simulation(Model::ModelInterface* robot, Simulation::Sai2Simulation* sim) {
	fSimulationRunning = true;

	// real forces and contact point
	std::string robot_name = "4PBOT";
	std::string link_name = "link3";
	// std::string link_name = "link4";
	std::vector<Eigen::Vector3d> point_list;
	std::vector<Eigen::Vector3d> force_list;
	Eigen::Vector3d sensed_force, sensed_moment;
	Eigen::Vector2d real_forces;
	Eigen::Vector3d contact_pos;
	Eigen::Vector2d real_contact_position;

	// create a timer
	LoopTimer timer;
	timer.initializeTimer();
	timer.setLoopFrequency(2000); 
	double last_time = timer.elapsedTime(); //secs
	bool fTimerDidSleep = true;

	unsigned long long simulation_counter = 0;

	while (fSimulationRunning) {
		fTimerDidSleep = timer.waitForNextLoop();

		// integrate forward
		double curr_time = timer.elapsedTime();
		double loop_dt = curr_time - last_time; 
		sim->integrate(loop_dt);

		// update forces and contact point
		point_list.clear();
		force_list.clear();
		sim->getContactList(
			point_list,
			force_list,
			robot_name,
			link_name);
		// zero out current forces
		sensed_force.setZero();
		contact_pos.setZero();
		// if list is empty, simply set forces to 0
		if(not point_list.empty()) 
		{
			// transform to sensor frame
			for (uint pt_ind=0; pt_ind < point_list.size(); ++pt_ind) 
			{
				sensed_force += force_list[pt_ind];
				contact_pos += point_list[pt_ind];
			}
			contact_pos /= point_list.size();
		}
		real_forces = sensed_force.tail(2);
		real_contact_position = contact_pos.tail(2);

		// filter forces
		if(filter_forces)
		{
			Frf = (Frpp + 2*Frp + real_forces)/simulation_filter_gain + coeff_simulation_filter[0]*Frppf + coeff_simulation_filter[1]*Frpf;
			Frpp = Frp;
			Frp = real_forces;
			Frppf = Frpf;
			Frpf = Frf;
		
			Prf = (Prpp + 2*Prp + real_contact_position)/simulation_filter_gain + coeff_simulation_filter[0]*Prppf + coeff_simulation_filter[1]*Prpf;
			Prpp = Prp;
			Prp = real_contact_position;
			Prppf = Prpf;
			Prpf = Prf;
		}
		else
		{
			Frf = real_forces;
			Prf = real_contact_position;
		}


		// if(simulation_counter % 1000 == 0)
		// {
		// 	cout << "real contact forces : " << Frf.transpose() << endl;
		// 	cout << "real contact pos : " << real_contact_position.transpose() << endl;
		// 	cout << endl;
		// }

		// if (!fTimerDidSleep) {
		// 	cout << "Warning: timer underflow! dt = 0.001;: " << loop_dt << "\n";
		// }

		// log real force and position values to file
		real_file << timer.elapsedTime() << '\t' << Frf.transpose() << '\t' << real_contact_position.transpose() << endl;

		//update last time
		last_time = curr_time;

		simulation_counter++;
	}

	double end_time = timer.elapsedTime();
    std::cout << "\n";
    std::cout << "Simulation Loop run time  : " << end_time << " seconds\n";
    std::cout << "Simulation Loop updates   : " << timer.elapsedCycles() << "\n";
    std::cout << "Simulation Loop frequency : " << timer.elapsedCycles()/end_time << "Hz\n";
}


//------------------------------------------------------------------------------
GLFWwindow* glfwInitialize() {
		/*------- Set up visualization -------*/
    // set up error callback
    glfwSetErrorCallback(glfwError);

    // initialize GLFW
    glfwInit();

    // retrieve resolution of computer display and position window accordingly
    GLFWmonitor* primary = glfwGetPrimaryMonitor();
    const GLFWvidmode* mode = glfwGetVideoMode(primary);

    // information about computer screen and GLUT display window
	int screenW = mode->width;
    int screenH = mode->height;
    int windowW = 0.8 * screenH;
    int windowH = 0.5 * screenH;
    int windowPosY = (screenH - windowH) / 2;
    int windowPosX = windowPosY;

    // create window and make it current
    glfwWindowHint(GLFW_VISIBLE, 0);
    GLFWwindow* window = glfwCreateWindow(windowW, windowH, "SAI2.0 - CS238", NULL, NULL);
	glfwSetWindowPos(window, windowPosX, windowPosY);
	glfwShowWindow(window);
    glfwMakeContextCurrent(window);
	glfwSwapInterval(1);

	return window;
}

//------------------------------------------------------------------------------

void glfwError(int error, const char* description) {
	cerr << "GLFW Error: " << description << endl;
	exit(1);
}

//------------------------------------------------------------------------------

void keySelect(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	bool set = (action != GLFW_RELEASE);
    switch(key) {
		case GLFW_KEY_ESCAPE:
			// exit application
			glfwSetWindowShouldClose(window,GL_TRUE);
			break;
		case GLFW_KEY_RIGHT:
			fTransXp = set;
			break;
		case GLFW_KEY_LEFT:
			fTransXn = set;
			break;
		case GLFW_KEY_UP:
			fTransYp = set;
			break;
		case GLFW_KEY_DOWN:
			fTransYn = set;
			break;
		case GLFW_KEY_A:
			fTransZp = set;
			break;
		case GLFW_KEY_Z:
			fTransZn = set;
			break;
		default:
			break;
    }
}

//------------------------------------------------------------------------------

void mouseClick(GLFWwindow* window, int button, int action, int mods) {
	bool set = (action != GLFW_RELEASE);
	//TODO: mouse interaction with robot
		switch (button) {
		// left click pans and tilts
		case GLFW_MOUSE_BUTTON_LEFT:
			fRotPanTilt = set;
			// NOTE: the code below is recommended but doesn't work well
			// if (fRotPanTilt) {
			// 	// lock cursor
			// 	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
			// } else {
			// 	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
			// }
			break;
		// if right click: don't handle. this is for menu selection
		case GLFW_MOUSE_BUTTON_RIGHT:
			//TODO: menu
			break;
		// if middle click: don't handle. doesn't work well on laptops
		case GLFW_MOUSE_BUTTON_MIDDLE:
			break;
		default:
			break;
	}
}