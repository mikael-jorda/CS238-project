#include <iostream>
#include <string>
#include <thread>
#include <math.h>

#include "model/ModelInterface.h"
#include "graphics/ChaiGraphics.h"
#include "simulation/Sai2Simulation.h"
#include <dynamics3d.h>

#include "force_sensor/ForceSensorSim.h"
#include "timer/LoopTimer.h"

#include <GLFW/glfw3.h> //must be loaded after loading opengl/glew as part of graphicsinterface

#include <signal.h>
bool fSimulationRunning = false;
void sighandler(int){fSimulationRunning = false;}

using namespace std;

const string world_file = "../resources/01-model_based_baseline/world.urdf";
const string robot_file = "../resources/01-model_based_baseline/4pbot_fixed.urdf";
const string robot_name = "4PBOT";

const string camera_name = "camera_fixed";

Eigen::Vector3d robot_origin = Eigen::Vector3d(0.0, -1.5, 0.050001);

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
	// load robots
	auto robot = new Model::ModelInterface(robot_file, Model::rbdl, Model::urdf, false);

	// load simulation world
	auto sim = new Simulation::Sai2Simulation(world_file, Simulation::urdf, false);
	sim->setCollisionRestitution(0);
	sim->setCoeffFrictionStatic(0);

	

	// set initial condition
	robot->_q << -50.0/180.0*M_PI,
				 -25.0/180.0*M_PI,
				 -35.0/180.0*M_PI,
				  10.0/180.0*M_PI;
	sim->setJointPositions(robot_name, robot->_q);
	robot->updateModel();

	// initialize GLFW window
	GLFWwindow* window = glfwInitialize();

	double last_cursorx, last_cursory;

    // set callbacks
	glfwSetKeyCallback(window, keySelect);
	glfwSetMouseButtonCallback(window, mouseClick);

	// start the simulation thread first
	fSimulationRunning = true;
	thread sim_thread(simulation, robot, sim);

	// next start the control thread
	thread ctrl_thread(control, robot, sim);
	
    // while window is open:
    while (!glfwWindowShouldClose(window)) {
		// update kinematic models
		// robot->updateModel();

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
	Eigen::MatrixXd J3d, Jv2d, Jw2d, J2d, Lambda, N, Jbar;
	J3d = Eigen::MatrixXd::Zero(6, dof);
	Jw2d = Eigen::MatrixXd::Zero(1, dof);
	Jv2d = Eigen::MatrixXd::Zero(2, dof);
	J2d = Eigen::MatrixXd::Zero(3, dof);
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
	Eigen::MatrixXd Ji = Eigen::MatrixXd(3,dof);
	Eigen::Vector2d estimated_Fc;
	Eigen::Vector2d estimated_rc;
	Eigen::Matrix3d Ri3d;
	Eigen::Matrix2d Ri;
	Eigen::Vector3d pos_i3d;
	Eigen::Vector2d pos_i;

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
		pos_i.setZero();
		collision_index = 0;
		for(int i=0; i<dof; i++)
		{
			if(r(i) > eps_mu)
			{
				mu(i) = r(i);
				collision_index = i+1;
			}
		}
		switch(collision_index)
		{
			// case 1 :
			// robot->J_0(J3d, "link1", Eigen::Vector3d::Zero());
			// break;

			// case 2 :
			// robot->J_0(J3d, "link2", Eigen::Vector3d::Zero());
			// break;

			case 3 :
			robot->J_0(J3d, "link3", Eigen::Vector3d::Zero());
			robot->rotation(Ri3d, "link3");
			robot->position(pos_i3d, "link3", Eigen::Vector3d::Zero());
			break;

			case 4 :
			robot->J_0(J3d, "link4", Eigen::Vector3d::Zero());
			robot->rotation(Ri3d, "link4");
			robot->position(pos_i3d, "link4", Eigen::Vector3d::Zero());
			break;

			default :
			J3d.setZero();
		}

		if(collision_index > 2 and collision_index <= dof)
		{
			// define matrices
			Ri = Ri3d.block(1,1,2,2);
			Ji = J3d.block(1,0,3,dof);
			// Ji.block(0,0,2,dof) = J3d.block(1,0,2,dof);
			// Ji.block(2,0,1,dof) = J3d.block(3,0,1,dof);

			// find equivqlent force at the joint frame
			Fi = (Ji.transpose()).colPivHouseholderQr().solve(mu);

			// deduce applied force
			estimated_Fc = Fi.head(2);
			estimated_Fc = Ri.transpose()*estimated_Fc;

			// find contact point via moment
			double ry = 0;
			if(estimated_Fc(0) > 0)
			{
				ry = -0.05;
			}
			else
			{
				ry = 0.05;
			}
			double rz = (ry*estimated_Fc(1) - Fi(2))/estimated_Fc(0);
			estimated_rc << ry, rz;

			// transform to base frame
			pos_i = pos_i3d.tail(2);
			estimated_rc = robot_origin.tail(2) + pos_i + Ri*estimated_rc;
			estimated_Fc = Ri*estimated_Fc;
		}

		// jacobian for control
		robot->J_0(J3d, link_name, pos_in_link);
		Jv2d = J3d.block(1,0,2,dof);
		robot->operationalSpaceMatrices(Lambda, Jbar, N, Jv2d);

		// -------------------------------------------
		////////////////////////////// Compute joint torques
		double time = controller_counter/control_freq;


		//----------- pos task
		robot->position(pos3d, link_name, pos_in_link);
		x = pos3d.tail(2);
		xdot = Jv2d*robot->_dq;

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
		if(controller_counter % 500 == 0)
		{
			// cout << "mu : " << mu.transpose() << endl;
			cout << "Fc : " << estimated_Fc.transpose() << endl;
			cout << "rc : " << estimated_rc.transpose() << endl;
			// cout << "Ji : " << Ji << endl;
			cout << endl;
		}
		if(controller_counter == 3000)
		{
			gpjs = false;
		}

		controller_counter++;

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

		if(simulation_counter % 1000 == 0)
		{
			cout << "real contact forces : " << real_forces.transpose() << endl;
			cout << "real contact pos : " << real_contact_position.transpose() << endl;
			cout << endl;
		}

		// if (!fTimerDidSleep) {
		// 	cout << "Warning: timer underflow! dt = 0.001;: " << loop_dt << "\n";
		// }

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
    GLFWwindow* window = glfwCreateWindow(windowW, windowH, "SAI2.0 - CS327a HW2", NULL, NULL);
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