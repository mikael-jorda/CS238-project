#!/usr/bin/env python

# read redis keys and dump them to a file
import redis, time, signal, sys
import os
import json

runloop = True
counter = 0

# handle ctrl-C and close the files
def signal_handler(signal, frame):
	global runloop
	runloop = False
	print('  Exiting data logger')

signal.signal(signal.SIGINT, signal_handler)

# data files
folder = 'data/'
if not os.path.exists(folder):
    os.makedirs(folder)

# date and time
header = time.strftime("%x").replace('/','-') + '_' + time.strftime("%X").replace(':','-')

file_commanded_torques = open(folder + '/' + header + '_fgc.txt','w')

file_commanded_torques.write('Command Torques\n')

# open redis server
r_server = redis.StrictRedis(host='localhost', port=6379, db=0)

# redis keys used in SAI2
SIM_TIMESTAMP_KEY = "sai2::bracing2d::simulation::timestamp";
JOINT_TORQUES_LOGGER_KEY = "sai2::bracing2d::actuators::fgc_logger";

# data logging frequency
logger_frequency = 1000.0  # Hz
logger_period = 1.0/logger_frequency
t_init = time.time()
t = t_init

print 'Start Logging Data\n'

while(runloop):
	t += logger_period

	commanded_torques = json.loads(r_server.get(JOINT_TORQUES_LOGGER_KEY))

	line = r_server.get(SIM_TIMESTAMP_KEY) + '\t' + " ".join([str(x) for x in commanded_torques]) + '\n'
	# line = '0\t' + " ".join([str(x) for x in force]) + '\n'
	# print line
	file_commanded_torques.write(line)

	counter = counter + 1

	time.sleep(max(0.0,t-time.time()))

elapsed_time = time.time() - t_init
print "Elapsed time : ", elapsed_time, " seconds"
print "Loop cycles  : ", counter
print "Frequency    : ", counter/elapsed_time, " Hz"

file_commanded_torques.close()